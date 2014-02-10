#define W_TOL 0.00001
!#define AVOID_FFT

subroutine build_rho(states, betawfc, wg, wq, opt_basis, nspin, rho, becsum)
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
  use mp_global, only : intra_pool_comm
  USE wvfct, only: ecutwfc_int => ecutwfc
  use buffers, only : get_buffer

  IMPLICIT NONE
  type(nk_list), intent(in) :: states
  type(nk_list), intent(in) :: betawfc
  REAL(DP), intent(in) :: wg(:,:)
  REAL(DP), intent(in) :: wq(:)
  type(basis), intent(in) :: opt_basis
  integer, intent(in)   :: nspin
  type(scf_type), intent(inout) :: rho
  real(DP), intent(inout) :: becsum(:,:,:)

  ! locals
  integer :: npw, nbasis, nbnd, nk
  integer :: ibnd, ir, max_band
  complex(DP), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0)
  real(DP) :: w1
  real(DP), parameter :: k_gamma(3) = 0.d0
  complex(DP), allocatable :: tmp(:,:)
  complex(DP), pointer :: ptr(:,:) => NULL()
  complex(DP), allocatable :: betapsi(:,:)
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  integer :: ioff, a, t, ijh, ih, jh, k, spin

  npw   = size(opt_basis%elements, 1)
  nbnd  = size(states%host_ar(1)%dat, 2)
  nk    = states%nk / nspin
  allocate(igk(ngm), g2kin(ngm))
  call gk_sort(k_gamma, ngm, g, ecutwfc_int/tpiba2, npw, igk, g2kin)
  nbasis = opt_basis%length

  do k = 1, nk * nspin
    spin = (k-1) / nk + 1
  ! Form a charge density
  !
  ! ... here we sum for each k point the contribution
  ! ... of the wavefunctions to the charge
  max_band = 1
  DO ibnd = 1, nbnd
      !
      if (ibnd > size(wg,1)) exit 
      if (abs(wg(ibnd,k) / wq(k)) < W_TOL) cycle
      max_band = ibnd
  enddo

  if (size(states%host_ar) == 1) then
    ptr => states%host_ar(1)%dat
    call get_buffer(ptr, nbasis*nbnd, states%file_unit, k)
  else
    ptr => states%host_ar(k)%dat
  endif

#ifdef AVOID_FFT
  ! We have <r|B>, so we can just basis transform to get <r|psi>
  allocate(tmp(dffts%nnr, max_band))
  call start_clock( '  gemm')
  call ZGEMM('N', 'N', dffts%nnr, max_band, nbasis, one, &
              opt_basis%elem_rs, dffts%nnr, &
              ptr, nbasis, zero, &
              tmp, dffts%nnr)
  call stop_clock( '  gemm')
  ! Accumulate |<r|psi>|^2 with weights
  do ibnd = 1, nbnd
      if (abs(wg(ibnd,k) / wq(k)) < W_TOL) cycle
      w1 = wg(ibnd,k) / omega
      do ir = 1, dffts%nnr
        rho%of_r(ir, spin) = rho%of_r(ir, spin) + w1 * (DBLE(tmp(ir, ibnd))**2 + AIMAG(tmp(ir, ibnd))**2) 
      enddo
  enddo
#else
  ! FFTs are awesome beacuse they do basis transforms quickly.  Let's use them
  allocate(tmp(npw, max_band))

  ! Transform <B|psi> to <G|psi>
  call start_clock( '  gemm')
  call ZGEMM('N', 'N', npw, max_band, nbasis, one, &
              opt_basis%elements, npw, &
              ptr, nbasis, zero, &
              tmp, npw)
  call stop_clock( '  gemm')

  do ibnd = 1, nbnd
      if (abs(wg(ibnd,k) / wq(k)) < W_TOL) cycle
      w1 = wg(ibnd,k) / omega
      psic(:) = ( 0.D0, 0.D0 )

      ! Transform <G|psi> to <r|psi>
      call start_clock('  fft')
      psic(nls(igk(1:npw))) = tmp(1:npw, ibnd)
      CALL invfft ('Wave', psic, dffts)

      ! Accumulate |<r|psi>|^2 with weights
      call stop_clock('  fft')
      do ir = 1, dffts%nnr
        rho%of_r(ir, spin) = rho%of_r(ir, spin) + w1 * (DBLE(psic(ir))**2 + AIMAG(psic(ir))**2) 
      enddo
  enddo
#endif

  deallocate(tmp)
  enddo 

  deallocate(igk, g2kin)

  ! ==================================================================
  ! Add the non-local part
  ! ==================================================================
  if (size(betawfc%host_ar) == 0) return 

  call start_clock('  addproj')

  kpoint: do k = 1, nk * nspin
    spin = (k-1) / nk + 1

  if (size(betawfc%host_ar) == 1) then
    ptr => betawfc%host_ar(1)%dat
    call get_buffer(ptr, nkb*nbnd, betawfc%file_unit, k)
  else
    ptr => betawfc%host_ar(k)%dat
  endif

  band: DO ibnd = 1, nbnd
    !
    if (abs(wg(ibnd,k) / wq(k)) < W_TOL) cycle
    w1 = wg(ibnd,k) !/ omega
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

  enddo kpoint
  call stop_clock('  addproj')

  nullify(ptr)

  return

end subroutine build_rho
