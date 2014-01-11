!
! This routine takes the bloch states and projectors from the qpoint loop and stores
! them for later use in build_dens
!
#define __SSDIAG
subroutine store_states(wfc, pp, k, states, betawfc)
  use kinds, only : DP
  use uspp_param, only : nh
  use srb_types, only : pseudop, nk_list
  use srb_matrix, only : print_dmat, dmat
  use scalapack_mod, only : scalapack_localindex
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm
  use buffers, only : open_buffer, save_buffer, close_buffer

  IMPLICIT NONE

  type(dmat), intent(in)    :: wfc
  type(pseudop), intent(in)    :: pp
  integer, intent(in)        :: k
  type(nk_list), intent(inout)   :: states
  type(nk_list), intent(inout)   :: betawfc

  integer :: nbasis, nbnd, nkb
  integer :: i, j, i_l, j_l
  integer :: t
  logical :: islocal
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
  COMPLEX(DP), pointer :: tmp(:,:) => NULL(), tmp2(:,:) => NULL()
  logical stat
  integer :: ptr

  integer,save :: old_states_length

  nbasis = size(states%host_ar(1)%dat,  1)
  nbnd   = size(states%host_ar(1)%dat,  2)

  ! check to see if we're caching to file
  if (size(states%host_ar) == 1) then
    if (states%file_unit < 0) then
      states%file_unit = 4826
      call open_buffer(states%file_unit, 'states', nbasis*nbnd, 1, stat)
      old_states_length = nbasis*nbnd
    endif 
    if (nbasis*nbnd > old_states_length) then
      call close_buffer(states%file_unit,'delete')
      call open_buffer(states%file_unit, 'states', nbasis*nbnd, 1, stat)
      old_states_length = nbasis*nbnd
    endif
    call save_buffer(wfc%dat, nbasis*nbnd, states%file_unit, k)
  else
    states%host_ar(k)%dat = wfc%dat   
  endif

!  if (abs(sum(dble(conjg(tmp(:,1)) * tmp(:,1))) - 1) > 1.d-6) &
!  write(*,*) "NORM: ", sum(dble(conjg(tmp(:,1)) * tmp(:,1)))

  if (pp%us) then

  ! check to see if we're caching to file
  if (size(betawfc%host_ar) == 1) then
    ptr = 1
  else
    ptr = k
  endif

  ! produce <\beta|\psi>
  do t = 1, size(pp%na)
    call print_dmat(pp%projs(t))
    call print_dmat(states%host_ar(1))
    call print_dmat(betawfc%host_ar(1))
    call pZGEMM('C', 'N', pp%na(t)*nh(t), states%host_ar(1)%desc(4), nbasis, &
               one,  pp%projs(t)%dat, 1, 1, pp%projs(t)%desc, &
                     wfc%dat, 1, 1, states%host_ar(1)%desc, &
               zero, betawfc%host_ar(ptr)%dat, pp%na_off(pp%nt_off(t)), 1, betawfc%host_ar(1)%desc)
  enddo

  if (size(betawfc%host_ar) == 1) then
    if (betawfc%file_unit < 0) then
      betawfc%file_unit = 4827
      call open_buffer(betawfc%file_unit, 'bstates', size(betawfc%host_ar(1)%dat), 1, stat)
    endif  
    call save_buffer(betawfc%host_ar(1)%dat, size(betawfc%host_ar(1)%dat), betawfc%file_unit, k)
  endif

  endif


  return

end subroutine store_states
