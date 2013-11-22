subroutine build_s_matrix(pp, q, S_matrix)
  USE kinds, ONLY : DP

  use srb_types, only : pseudop
  USE uspp, only : nkb, qq
  use uspp_param, only : nh, nhm
  use ions_base, only : nat, ityp, nsp
  use scalapack_mod, only : scalapack_localindex
  use buffers, only : open_buffer, save_buffer, get_buffer, close_buffer

  IMPLICIT NONE

  type(pseudop), intent(inout) :: pp
  integer,       intent(in)    :: q
  complex(DP),   intent(out)   :: S_matrix(:,:)
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

#define __SSDIAG
#ifndef __SSDIAG
  complex(DP), allocatable :: S_tmp(:,:)
  nbasis = size(pp%projs, 1)
  allocate(S_tmp(nbasis, nbasis))
  S_tmp = cmplx(0.d0, kind=DP)
#endif

  nbasis = size(pp%projs, 1)

  ! if the save file hasn't been opened, open one
  if (pp%s_unit < 0) then
    pp%s_unit = - pp%s_unit
    call open_buffer(pp%s_unit, 's_matrix', size(S_matrix), 1, info)
    old_size_s_matrix=size(S_matrix)
  endif

  ! if the buffer increases in size - dynamically change it
  if (size(S_matrix) > old_size_s_matrix) then
    call close_buffer(pp%s_unit,'delete')
    call open_buffer(pp%s_unit, 's_matrix', size(S_matrix), 1, info)
    old_size_s_matrix=size(S_matrix)
  endif

  if (q > 0) then
    call get_buffer(S_matrix, size(S_matrix), pp%s_unit, q)
    return
  endif

  allocate(S_half(nbasis, nkb), rwork(2*nhm*nbasis))


  ! Make <B|S|B>
  S_matrix = cmplx(0.d0, kind=DP)
  ioff = 1 !index offset
  do t = 1, nsp
   do a = 1, nat
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
  ! Do the right side of the transformation, summing into S_matrix
#ifdef __SSDIAG
  call ZGEMM('N', 'C', nbasis, nbasis, nkb, one, &
               S_half(:,:), nbasis, &
               pp%projs(:,:), nbasis, zero, &
               S_matrix(:,:), nbasis)
#else
  call ZGEMM('N', 'C', nbasis, nbasis, nh(t), one, &
               S_half(:,:), nbasis, &
               pp%projs(:,:), nbasis, zero, &
               S_tmp(:,:), nbasis)
#endif
  deallocate(S_half, rwork)

#ifdef __SSDIAG
  ! add identity to S
  forall(i = 1:nbasis) S_matrix(i,i) = S_matrix(i,i) + one
#else
  do i = 1, nbasis; do j = 1, nbasis
    call scalapack_localindex(i, j, i_l, j_l, islocal)
    if (islocal) then
      if (i == j) then
        S_matrix(i_l, j_l) = S_tmp(i, j) + one
      else
        S_matrix(i_l, j_l) = S_tmp(i, j)
      endif
    endif
  enddo; enddo
  deallocate(S_tmp)
#endif

  call save_buffer(S_matrix, size(S_matrix), pp%s_unit, -q)

  return

end subroutine build_s_matrix

