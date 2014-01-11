#define W_TOL 0.00001

subroutine build_rho_reduced(states, betawfc, wg, wq, nspin, rho, becsum)
  use kinds, only : DP
  use srb_types, only : basis, nk_list
  use srb_matrix, only : dmat, block_outer, col_scal, copy_dmat, print_dmat
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
  type(dmat), intent(inout) :: rho(:)
  real(DP), intent(inout) :: becsum(:,:,:)

  ! locals
  integer :: npw, nbasis, nbnd, nk
  integer :: ibnd, ibnd_g, ir
  complex(DP), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0)
  real(DP) :: w1
  complex(DP) :: wz
  real(DP), parameter :: k_gamma(3) = 0.d0
  type(dmat) :: tmp
  complex(DP), allocatable :: betapsi(:,:)
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  real(DP), allocatable :: root(:)
  integer :: ioff, a, t, ijh, ih, jh, i, j, k, spin, q, ptr
  integer :: max_band
  real(DP) :: trace
  integer, external :: indxl2g

  nbasis = size(states%host_ar(1)%dat, 1)
  nbnd   = size(states%host_ar(1)%dat, 2)
  nk     = states%nk / nspin

  allocate(root(nbasis))
  call copy_dmat(tmp, states%host_ar(1))

  trace = 0.
  do k = 1+me_pool, nk, nproc_pool
   do spin = 1, nspin
    q = (k+(spin-1)*(nk+nproc_pool)-1)/nproc_pool + 1
    

  ! Form a charge density
  !
  ! ... here we sum for each k point the contribution
  ! ... of the wavefunctions to the charge
    if (size(states%host_ar) == 1) then
      ptr = 1
      call get_buffer(states%host_ar(1)%dat, nbasis*nbnd, states%file_unit, q)
    else
      ptr = q
    endif
#if 1
    call col_scal(states%host_ar(ptr), tmp, wg(:,k+(spin-1)*nk)/omega)
    call block_outer(nbasis, nbnd, &
                     one, tmp%dat, nbasis, &
                          states%host_ar(ptr)%dat, nbasis, &
                     one, rho(spin))
    trace = trace + sum(wg(1:5,k+(spin-1)*nk)/omega)

#else
    max_band = 1
    do ibnd = 1, nbnd
      if (wg(ibnd, k+(spin-1)*nk) / wq(k+(spin-1)*nk) < W_TOL) exit
      max_band = ibnd
    enddo

    root(1:max_band) = sqrt(wg(1:max_band, k+(spin-1)*nk) / omega)
    forall(j = 1:max_band, i = 1:nbasis) tmp%dat(i,j) = states%host_ar(ptr)%dat(i,j) * root(j)

    call zherk('U', 'N', nbasis, max_band, one, &
               tmp%dat, nbasis, one, &
               rho(spin)%dat, nbasis)

    do ibnd = max_band + 1, nbnd
      if (ibnd > size(wg,1)) exit 
      if (abs(wg(ibnd, k+(spin-1)*nk) / wq(k+(spin-1)*nk)) < W_TOL) cycle

      wz = cmplx(wg(ibnd, k+(spin-1)*nk) / omega, kind=DP)
      call ZHERK('U', 'N', nbasis, 1, wz, &
                 tmp%dat(:, ibnd), nbasis, one, &
                 rho(spin)%dat, nbasis)
    enddo
#endif
   enddo
  enddo
  write(*,*) "Trace(rho) = ", trace
  do spin = 1, nspin
     call mp_sum(rho(spin)%dat, intra_pool_comm)
  enddo

  deallocate(root)

  ! ==================================================================
  ! Add the non-local part
  ! ==================================================================
  if (size(betawfc%host_ar) == 0) return 
  call start_clock('  addproj')

  kpoint: do k = 1+me_pool, nk, nproc_pool
   do  spin = 1,nspin 
    q = (k+(spin-1)*(nk+nproc_pool)-1)/nproc_pool + 1

  if (size(betawfc%host_ar) == 1) then
    ptr = 1 
    call get_buffer(betawfc%host_ar(1)%dat, nkb*nbnd, betawfc%file_unit, q)
  else
    ptr = q
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
          w1 * DBLE(betawfc%host_ar(ptr)%dat(ioff + ih, ibnd) * CONJG(betawfc%host_ar(ptr)%dat(ioff + ih, ibnd)))
        proj2: do jh = (ih + 1), nh(t)
          ijh = ijh + 1
          becsum(ijh, a, spin) = becsum(ijh, a, spin) + &
            2.d0 * w1 * DBLE(betawfc%host_ar(ptr)%dat(ioff + jh, ibnd) * CONJG(betawfc%host_ar(ptr)%dat(ioff + ih, ibnd)))
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


  return

end subroutine build_rho_reduced
