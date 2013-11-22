subroutine build_projs(opt_basis, xq, nq, pp)
  USE kinds, ONLY : DP
  USE srb_types, ONLY : basis, pseudop

  USE uspp, only : nkb
  use uspp_param, only : nh, nhm
  USE gvect, ONLY : ngm, g
  USE cell_base, ONLY : tpiba2, bg
  use wvfct, only : npwx_int => npwx, ecutwfc
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm, me_pool, nproc_pool
  use buffers, only : open_buffer, save_buffer, get_buffer, close_buffer
  IMPLICIT NONE

  TYPE(basis), intent(in) :: opt_basis
  REAL(DP),    intent(in) :: xq(:,:)
  integer,     intent(in) :: nq
  type(pseudop), intent(inout) :: pp

  ! locals
  integer q, q2, left, npw, npwx_tmp
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  real(DP) :: k_gamma(3)
  complex(DP), allocatable :: vkb(:,:)
  complex(DP), allocatable :: projs_l(:,:)
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP)
  complex(DP), parameter :: one = cmplx(1.d0, kind=DP)
  real(DP) :: qcart(3), qtmp(3)
  logical ::  info

  integer,save :: old_basis_length

  ! projectors are q-dependent [:-(], so need to generate gvector info
  npw = size(opt_basis%elements, 1)
  k_gamma = 0.d0

  call start_clock('  proj_init')
  allocate(igk(ngm), g2kin(ngm))
  call gk_sort(k_gamma(:), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
  ! construct <G|\beta>

  ! if the save file hasn't been opened, open one
  if (pp%projs_unit < 0) then
    pp%projs_unit = - pp%projs_unit
    call open_buffer(pp%projs_unit, 'projs', opt_basis%length * nkb, 1, info)
    old_basis_length = opt_basis%length
  endif

  ! if we have a longer buffer - adjust size of file dynamically
  if (opt_basis%length > old_basis_length) then
    call close_buffer(pp%projs_unit,'delete')
    call open_buffer(pp%projs_unit, 'projs', opt_basis%length * nkb, 1, info)
    old_basis_length = opt_basis%length
  endif


#define BLOCK 8
  allocate(vkb(npw, nkb*BLOCK))
  allocate(projs_l(opt_basis%length, nkb*BLOCK))

    npwx_tmp = npwx_int
    npwx_int = npw

  do q = 1, nq, BLOCK
   left = MIN(nq - q + 1, BLOCK)
   do q2 = 0, left-1
    qtmp = xq(:,q+q2) - floor(xq(:,q+q2)) 
    qcart = matmul( bg, qtmp ) 
  
    call start_clock('  proj_init')
    call init_us_2_srb(npw, igk, qcart, -1, vkb(:,1+q2*nkb:q2*nkb+nkb)) 
    call stop_clock('  proj_init')
   enddo

    ! construct <B|\beta> = <B|G><G|\beta>
    call start_clock('  proj_gemm')
    call ZGEMM('C', 'N', opt_basis%length, nkb*left, npw, one, &
               opt_basis%elements, npw, & 
               vkb, npw, zero, &
               projs_l, opt_basis%length) 
    call stop_clock('  proj_gemm')

    call start_clock('  proj_comm')
    call mp_sum(projs_l(:,1:left*nkb), intra_pool_comm) 
    call stop_clock('  proj_comm')


    call start_clock('  proj_save')
    ! save the projs
   do q2 = 0, left-1
    if (MOD(q+q2-1, nproc_pool) == me_pool) then
      call save_buffer(projs_l(:,1+q2*nkb:q2*nkb+nkb), opt_basis%length * nkb, pp%projs_unit, (q+q2-1)/nproc_pool+1)
    endif
   enddo
    call stop_clock('  proj_save')

  enddo

  deallocate(vkb, projs_l)
  deallocate(igk, g2kin)
  npwx_int = npwx_tmp

  return

end subroutine build_projs

