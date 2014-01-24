subroutine load_projs(q, pp)
  USE kinds, ONLY : DP
  USE srb_types, ONLY : pseudop

  use mp_global, only : npot
  use buffers, only : get_buffer
  IMPLICIT NONE

  integer,     intent(in) :: q
  type(pseudop), intent(inout) :: pp
  integer :: t
  
  do t = 1, pp%ntyp
    if (size(pp%projs(t)%dat) > 0) &
      call get_buffer(pp%projs(t)%dat, size(pp%projs(t)%dat), pp%p_unit(t), (q-1)/npot+1)
  enddo

  return

end subroutine load_projs

