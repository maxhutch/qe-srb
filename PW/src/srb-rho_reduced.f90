#define W_TOL 0.00001

subroutine build_rho_reduced(states, betawfc, wg, wq, nspin, rho, becsum)
  use kinds, only : DP
  use srb_types, only : basis, nk_list
  use scf, only : scf_type
  use cell_base, only : omega, tpiba2
  use uspp, only : nkb
  use uspp_param, only : nh, nhm
  use ions_base, only : nat, ityp, nsp

  use gvecs, only : nls
  use gvect, only : ngm, g, nl
  use wavefunctions_module, only : psic
  use fft_base, only : dfftp, dffts
  use fft_interfaces, only : invfft, fwfft
  use symme, only : sym_rho 
  use scalapack_mod, only : scalapack_localindex
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm, me_pool, nproc_pool
  USE wvfct, only: ecutwfc_int => ecutwfc
  use buffers, only : get_buffer


  IMPLICIT NONE

  type(nk_list), intent(in) :: states
  type(nk_list), intent(in) :: betawfc
  REAL(DP), intent(in) :: wg(:,:)
  REAL(DP), intent(in) :: wq(:)
  integer, intent(in)   :: nspin
  COMPLEX(DP), intent(inout) :: rho(:,:,:)
  real(DP), intent(inout) :: becsum(:,:,:)

  ! locals
  integer :: npw, nbasis, nbnd, nk
  integer :: ibnd, ir
  complex(DP), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0)
  real(DP) :: w1
  complex(DP) :: wz
  real(DP), parameter :: k_gamma(3) = 0.d0
  complex(DP), allocatable :: tmp(:,:)
  complex(DP), allocatable :: betapsi(:,:)
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  real(DP), allocatable :: root(:)
  integer :: ioff, a, t, ijh, ih, jh, i, j, k, spin, q
  integer :: max_band
  complex(DP), pointer :: ptr(:,:) => NULL()
  real(DP) :: trace

  nbasis = size(states%host_ar, 1)
  nbnd   = size(states%host_ar, 2)
  nk     = states%nk / nspin

  allocate(root(nbasis))
  allocate(tmp(states%desc%nrl, states%desc%ncl))

  trace = 0.
  do k = 1+me_pool, nk, nproc_pool
   do spin = 1, nspin
    q = (k+(spin-1)*(nk+nproc_pool)-1)/nproc_pool + 1
    
    if (size(states%host_ar, 3) == 1) then
      ptr => states%host_ar(:,:,1)
      call get_buffer(ptr, nbasis*nbnd, states%file_unit, q)
    else
      allocate(ptr(nbasis,nbnd))
      ptr = states%host_ar(:,:,q)
    endif

  ! Form a charge density
  !
  ! ... here we sum for each k point the contribution
  ! ... of the wavefunctions to the charge

    
    max_band = 1
    do ibnd = 1, nbnd
      if (wg(ibnd, k+(spin-1)*nk) / wq(k+(spin-1)*nk) < W_TOL) exit
      max_band = ibnd
    enddo

    root(1:max_band) = sqrt(wg(1:max_band, k+(spin-1)*nk) / omega)
    forall(j = 1:max_band, i = 1:nbasis) ptr(i,j) = ptr(i,j) * root(j)

    call zherk('U', 'N', nbasis, max_band, one, &
               ptr(:,:), nbasis, one, &
               rho(:,:,spin), nbasis)

    do ibnd = max_band + 1, nbnd
      if (ibnd > size(wg,1)) exit 
      if (abs(wg(ibnd, k+(spin-1)*nk) / wq(k+(spin-1)*nk)) < W_TOL) cycle

      wz = cmplx(wg(ibnd, k+(spin-1)*nk) / omega, kind=DP)
      call ZHERK('U', 'N', nbasis, 1, wz, &
                 ptr(:, ibnd), nbasis, one, &
                 rho(:,:,spin), nbasis)
    enddo

   enddo
  enddo
  call mp_sum(rho, intra_pool_comm)
  if (size(states%host_ar, 3) == 1) deallocate(ptr)
  deallocate(root)

  ! ==================================================================
  ! Add the non-local part
  ! ==================================================================
  if (size(betawfc%host_ar,3) == 0) return 
  call start_clock('  addproj')

  kpoint: do k = 1+me_pool, nk, nproc_pool
   do  spin = 1,nspin 
    q = (k+(spin-1)*(nk+nproc_pool)-1)/nproc_pool + 1

  if (size(betawfc%host_ar, 3) == 1) then
    ptr => betawfc%host_ar(:,:,1)
    call get_buffer(ptr, nkb*nbnd, betawfc%file_unit, q)
  else
    ptr => betawfc%host_ar(:,:,q)
  endif

  band: DO ibnd = 1, nbnd
    !
    if (abs(wg(ibnd, k+(spin-1)*nk) / wq(k+(spin-1)*nk)) < W_TOL) cycle
    w1 = wg(ibnd, k+(spin-1)*nk) !/ omega
    ioff = 0
    type: do t = 1, nsp 
    atom: do a = 1, nat
      if (ityp(a) /= t) cycle
      ijh = 0
      proj1: do ih = 1, nh(t)
        ijh = ijh + 1
        becsum(ijh, a, spin) = becsum(ijh, a, spin) + &
          w1 * DBLE(ptr(ioff + ih, ibnd) * CONJG(ptr(ioff + ih, ibnd)))
        proj2: do jh = (ih + 1), nh(t)
          ijh = ijh + 1
          becsum(ijh, a, spin) = becsum(ijh, a, spin) + &
            2.d0 * w1 * DBLE(ptr(ioff + jh, ibnd) * CONJG(ptr(ioff + ih, ibnd)))
        enddo proj2
      enddo proj1
      ioff = ioff + nh(t)
    enddo atom
    enddo type
  enddo band
   enddo
  enddo kpoint
  call mp_sum(becsum, intra_pool_comm)

  call stop_clock('  addproj')

  nullify(ptr)

  return

end subroutine build_rho_reduced
