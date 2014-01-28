
subroutine backload(qpoints, srb, states)
  use kinds, only : DP
  use srb_types, only : basis, nk_list
  use kpoint, only : k_list
  use buffers, only : get_buffer, save_buffer
  USE klist, ONLY : xk, nks
  use io_files, only : iunwfc, nwordwfc
  use lsda_mod, only : nspin, isk
  USE gvect, only : g, ngm, ig_l2g
  USE wvfct, only: ecutwfc
  use cell_base, only : bg, tpiba2, at
  use mp, only : mp_max, mp_sum
  use mp_global, only : me_pool, nproc_pool, root_pool, intra_pool_comm
  use mp_global, only : npot, my_pot_id
  use mp_wave, only : mergewf, splitwf

  implicit none

  type(k_list), intent(in) :: qpoints
  TYPE(basis), INTENT(IN) :: srb
  type(nk_list), intent(inout)   :: states

  integer :: npw_orig, npwx_orig, ngk_gamma, npwx_gamma, igwx
  integer, allocatable :: ngk_orig(:), igk_orig(:), igk_l2g_orig(:), igk_gamma(:), igk_l2g_gamma(:)
  real(DP), allocatable :: g2kin(:)

  complex(DP), allocatable :: evc(:,:), evc_tmp(:)
  complex(DP), pointer :: ptr(:,:)
  integer :: this_k, k, q, ibnd, s
  integer :: npw, nbnd, nbasis
  COMPLEX(DP), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0)

  npw = size(srb%elements, 1)
  nbasis = size(srb%elements, 2)
  nbnd = states%nbnd

  allocate(ngk_orig(nks))
  call n_plane_waves(ecutwfc, tpiba2, nks, xk, g, ngm, npwx_orig, ngk_orig)
  allocate(igk_orig(npwx_orig), g2kin(npwx_orig), igk_l2g_orig(npwx_orig))
  allocate(igk_gamma(npwx_orig), igk_l2g_gamma(npwx_orig))
  call gk_sort (xk(:,1), ngm, g, ecutwfc/tpiba2, ngk_gamma, igk_gamma, g2kin)
  write(*,*) npwx_orig, npw, ngk_orig(1), ngk_gamma
  igk_l2g_gamma(1:ngk_gamma) = ig_l2g(igk_gamma(1:ngk_gamma))
  deallocate(igk_gamma)

  allocate(evc(npwx_orig, states%nbnd))
  do q = 1, qpoints%nred
   do s = 1, nspin
    this_k = 0
    do k = 1, nks
      if (all(abs(qpoints%xr(:,q) - matmul(transpose(at), xk(:,k)) ) < 0.00001) .and. isk(k) == s) then
        this_k = k
        exit
      endif
    enddo
    if (this_k == 0) cycle
    write(*,*) "backloading ", q, " to ", this_k
    ! load <B|psi>
    if (MOD(q-1, npot) == my_pot_id) then
      if (size(states%host_ar) == 1) then
        ptr => states%host_ar(1)%dat
        call get_buffer(ptr, nbasis*nbnd, states%file_unit, (q+(s-1)*(qpoints%nred+npot)-1)/npot + 1)
      else
        ptr => states%host_ar((q+(s-1)*(qpoints%nred+npot)-1)/npot + 1)%dat
      endif
    else
      if (size(states%host_ar) == 1) then
        ptr => states%host_ar(1)%dat
      else
        allocate(ptr(nbasis,nbnd))
      endif
      ptr = cmplx(0.d0, kind=DP)
    endif
    call mp_sum(ptr, intra_pool_comm)

    ! transform to <G|psi> at Gamma
    evc = cmplx(0.d0, kind=DP)
    call zgemm('N', 'N', npw, nbnd, nbasis, one, &
               srb%elements, npw, &
               ptr, nbasis, zero, &
               evc, npwx_orig)

    ! re-order to proper k-point
    call gk_sort (xk(:,this_k), ngm, g, ecutwfc/tpiba2, npw_orig, igk_orig, g2kin)
    igk_l2g_orig(1:npw_orig) = ig_l2g(igk_orig(1:npw_orig))
    igwx = maxval(igk_l2g_orig)
    call mp_max(igwx, intra_pool_comm)
    allocate(evc_tmp(igwx))
    do ibnd = 1, nbnd
      evc_tmp = cmplx(0.d0, kind=DP)
      call mergewf(evc(1:ngk_gamma, ibnd), evc_tmp, ngk_gamma, igk_l2g_gamma, &
                 me_pool, nproc_pool, root_pool, intra_pool_comm)
      evc(:, ibnd) = cmplx(0.d0, kind=DP)
      call splitwf(evc(1:npw_orig, ibnd), evc_tmp, npw_orig,  igk_l2g_orig,  &
                 me_pool, nproc_pool, root_pool, intra_pool_comm)
    enddo
    deallocate(evc_tmp)

    ! resave them
    CALL save_buffer ( evc, nwordwfc, iunwfc, this_k)
   enddo
  enddo
  deallocate(evc)
  deallocate(igk_orig)
  deallocate(g2kin)
  deallocate(igk_l2g_orig)
  deallocate(igk_l2g_gamma)
  deallocate(ngk_orig)

end subroutine backload

