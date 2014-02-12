subroutine transform_rho(rho_srb, opt_basis, rho)
  use kinds, only : DP
  use srb_types, only : basis
  use srb_matrix, only : dmat, copy_dmat, diag, setup_dmat, pool_scope, serial_scope
  use srb_matrix, only : distribute
  use scf, only : scf_type
  use cell_base, only : omega, tpiba2
  use uspp, only : nkb
  use uspp_param, only : nh, nhm
  use ions_base, only : nat, ityp, nsp

  use srb, only : decomp_size
  use buffers, only : get_buffer, save_buffer, open_buffer
  use gvecs, only : nls
  use gvect, only : ngm, g, nl
  use wavefunctions_module, only : psic
  use fft_base, only : dfftp, dffts
  use fft_interfaces, only : invfft, fwfft
  use symme, only : sym_rho 
  use scalapack_mod, only : scalapack_localindex
  use mp, only : mp_sum
  use mp_global, only : me_image, inter_pot_comm, my_pot_id
  use input_parameters, only : dens_tol
  USE wvfct, only: ecutwfc_int => ecutwfc

  IMPLICIT NONE

  type(dmat), intent(inout) :: rho_srb(:)
  type(basis), intent(inout) :: opt_basis
  type(scf_type), intent(inout) :: rho

  ! locals
  integer :: npw, nbasis, nspin
  integer :: i, ibnd, ir, max_band
  complex(DP), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0)
  integer :: spin, offset, num
  complex(DP), allocatable :: tmp(:,:)

  real(DP), parameter :: k_gamma(3) = 0.d0
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)

  real(DP), allocatable :: S(:)
  type(dmat) :: sv, rho_big
  complex(DP), allocatable :: work(:)
  real(DP), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  integer :: lwork, lrwork, liwork
  real(DP) :: trace, trace2
  integer, save :: funit = -128
  logical :: info

  npw   = size(opt_basis%elements, 1)
  nbasis = opt_basis%length
  nspin = size(rho_srb)

  allocate(igk(ngm), g2kin(ngm))
  call gk_sort(k_gamma, ngm, g, ecutwfc_int/tpiba2, npw, igk, g2kin)

  allocate(S(nbasis))
  call setup_dmat(rho_big, nbasis, nbasis, scope_in = pool_scope)
  call setup_dmat(sv, nbasis, nbasis, scope_in = serial_scope)


  do spin = 1, nspin

  ! find the rank-1 decomposition
  trace = 0.d0
  do i = 1, nbasis
  trace = trace + abs(rho_srb(spin)%dat(i,i))
  enddo
  call mp_sum(rho_srb(spin)%dat, inter_pot_comm)
  call distribute(rho_srb(spin), rho_big, my_pot_id)

  call start_clock('  svd')
  call diag(rho_big, S, sv)
  call stop_clock('  svd')

  S = abs(S)
  trace = sum(S)
  max_band = 1
  do ibnd = 1, nbasis
      if (S(nbasis+1-ibnd)*(nbasis+1.-ibnd)/trace < dens_tol) exit
      max_band = ibnd
  enddo
  if (me_image == 0) write(*,*) max_band, " of ", nbasis

  ! transform the representative wave-functions to <G|
  allocate(tmp(npw, max_band))
  call start_clock('  gemm')
  call ZGEMM('N', 'N', npw, max_band, nbasis, one, &
              opt_basis%elements, npw, &
              sv%dat(:,nbasis+1-max_band), nbasis, zero, &
              tmp, npw)
  call stop_clock('  gemm')

  ! accumulate the left singular vectors
  do ibnd = 1, max_band
    psic(:) = ( 0.D0, 0.D0 )

!    call save_buffer(rho_srb(:,nbasis-max_band + ibnd, spin), nbasis, funit+1, ibnd)
    ! Transform <G|u> to <r|u>
    call start_clock('  fft')
    psic(nls(igk(1:npw))) = tmp(1:npw, ibnd)
    CALL invfft ('Wave', psic, dffts)
    call stop_clock('  fft')

!    call save_buffer(psic, dffts%nnr, funit, ibnd)

    ! Accumulate |<r|u>|^2 with weights
    call start_clock('  acc')
    do ir = 1, dffts%nnr
      rho%of_r(ir, spin) = rho%of_r(ir, spin) + S(nbasis-max_band+ibnd) * (DBLE(psic(ir))**2 + AIMAG(psic(ir))**2)
    enddo
    call stop_clock('  acc')
  enddo
  decomp_size = max_band

  deallocate(tmp)
  enddo 
  deallocate(S, igk, g2kin)

  return

end subroutine transform_rho
