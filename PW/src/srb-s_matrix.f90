subroutine build_s_matrix(pp, q, Hk)
  USE kinds, ONLY : DP

  use srb_types, only : pseudop, kproblem
  use srb_matrix, only : block_outer, add_diag
  USE uspp, only : nkb, qq
  use uspp_param, only : nh, nhm
  use ions_base, only : nat, ityp, nsp
  use scalapack_mod, only : scalapack_localindex
  use buffers, only : open_buffer, save_buffer, get_buffer, close_buffer
  use mp_global, only : nproc_pot

  IMPLICIT NONE

  type(pseudop),  intent(inout) :: pp
  integer,        intent(in)    :: q
  type(kproblem), intent(inout)   :: Hk
  ! locals
  complex(DP), allocatable :: S_half(:,:)
  real(DP), allocatable :: rwork(:)
  integer :: a, t, ioff
  integer :: i, j, i_l, j_l
  logical :: islocal
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP)
  complex(DP), parameter :: one = cmplx(1.d0, kind=DP)
  integer :: nbasis
  logical :: info

  integer,save :: old_size_s_matrix

  nbasis = size(pp%projs, 1)

  ! if the save file hasn't been opened, open one
  if (pp%s_unit < 0) then
    pp%s_unit = - pp%s_unit
    call open_buffer(pp%s_unit, 's_matrix', size(Hk%S), 1, info)
    old_size_s_matrix=size(Hk%S)
  endif

  ! if the buffer increases in size - dynamically change it
  if (size(Hk%S) > old_size_s_matrix) then
    call close_buffer(pp%s_unit,'delete')
    call open_buffer(pp%s_unit, 's_matrix', size(Hk%S), 1, info)
    old_size_s_matrix=size(Hk%S)
  endif

  ! if saved, just reload it
  if (q > 0) then
    call get_buffer(Hk%S, size(Hk%S), pp%s_unit, q)
    return
  endif

  ! Make <B|S|B>
  allocate(S_half(nbasis, pp%nkb_l), rwork(2*nhm*nbasis))
  Hk%S = cmplx(0.d0, kind=DP)
  ioff = 1 !index offset
  do t = 1, nsp
   do a = 1, nat, nproc_pot
    if (ityp(a) /= t) cycle
    ! Do the left side of the transformation
    call zlacrm(nbasis, nh(t), &
                pp%projs(:,ioff:ioff+nh(t)), nbasis, &
                qq(:,:,t), nhm, &
                S_half(:, ioff:ioff+nh(t)), nbasis, &
                rwork)
    ioff = ioff + nh(t)
   enddo
  enddo
  ! Do the right side of the transformation, summing into Hk%S
  call block_outer(nbasis, pp%nkb_l, &
                   one,  S_half, nbasis, &
                         pp%projs, nbasis, &
                   zero, Hk%S, Hk%desc) 
  deallocate(S_half, rwork)

  call add_diag(nbasis, one, Hk%S, Hk%desc)

  call save_buffer(Hk%S, size(Hk%S), pp%s_unit, -q)

  return

end subroutine build_s_matrix

