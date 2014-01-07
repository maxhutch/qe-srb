!
! This routine takes the bloch states and projectors from the qpoint loop and stores
! them for later use in build_dens
!
#define __SSDIAG
subroutine store_states(wfc, projs, k, states, betawfc)
  use kinds, only : DP
  use srb_types, only : pseudop, nk_list
  use scalapack_mod, only : scalapack_localindex
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm
  use buffers, only : open_buffer, save_buffer, close_buffer

  IMPLICIT NONE

  complex(DP), intent(in)    :: wfc(:,:)
  type(pseudop), intent(in)    :: projs
  integer, intent(in)        :: k
  type(nk_list), intent(inout)   :: states
  type(nk_list), intent(inout)   :: betawfc

  integer :: nbasis, nbnd, nkb
  integer :: i, j, i_l, j_l
  logical :: islocal
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
  COMPLEX(DP), pointer :: tmp(:,:) => NULL(), tmp2(:,:) => NULL()
  logical stat

  integer,save :: old_states_length

  nbasis = size(states%host_ar,  1)
  nkb    = size(betawfc%host_ar, 1)
  nbnd   = states%nbnd

  ! check to see if we're caching to file
  if (size(states%host_ar, 3) == 1) then
    tmp => states%host_ar(:,:,1)
  else
    tmp => states%host_ar(:,:,k)
  endif 

  tmp = wfc(:,1:nbnd)

  if (size(states%host_ar, 3) == 1) then
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
    call save_buffer(tmp, nbasis*nbnd, states%file_unit, k)
  endif

!  if (abs(sum(dble(conjg(tmp(:,1)) * tmp(:,1))) - 1) > 1.d-6) &
!  write(*,*) "NORM: ", sum(dble(conjg(tmp(:,1)) * tmp(:,1)))

  if (projs%us) then

  ! check to see if we're caching to file
  if (size(betawfc%host_ar, 3) == 1) then
    tmp2 => betawfc%host_ar(:,:,1)
  else
    tmp2 => betawfc%host_ar(:,:,k)
  endif

  ! produce <\beta|\psi>
  call ZGEMM('C', 'N', nkb, nbnd, nbasis, one, &
             projs%projs, nbasis, &
             tmp, nbasis, zero, &
             tmp2, nkb)

  if (size(betawfc%host_ar, 3) == 1) then
    if (betawfc%file_unit < 0) then
      betawfc%file_unit = 4827
      call open_buffer(betawfc%file_unit, 'bstates', nkb*nbnd, 1, stat)
    endif  
    call save_buffer(tmp2, nkb*nbnd, betawfc%file_unit, k)
  endif

  endif

  nullify(tmp)
  nullify(tmp2)

  return

end subroutine store_states
