subroutine load_projs(q, pp)
  USE kinds, ONLY : DP
  USE srb_types, ONLY : pseudop

  use mp_global, only : npot
  use buffers, only : get_buffer
  IMPLICIT NONE

  integer,     intent(in) :: q
  type(pseudop), intent(inout) :: pp
  integer :: t
  
  call start_clock('  proj_load')
#if 1
  do t = 1, pp%ntyp
    call get_buffer(pp%projs(t)%dat, size(pp%projs(t)%dat), pp%p_unit(t), (q-1)/npot+1)
  enddo
#else
  call get_buffer(pp%projs, size(pp%projs), pp%projs_unit, q)
#endif
  call stop_clock('  proj_load')

  return

end subroutine load_projs

