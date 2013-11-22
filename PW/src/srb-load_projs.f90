subroutine load_projs(q, pp)
  USE kinds, ONLY : DP
  USE srb_types, ONLY : pseudop

  use mp_global, only : nproc_pool 
  use buffers, only : get_buffer
  IMPLICIT NONE

  integer,     intent(in) :: q
  type(pseudop), intent(inout) :: pp


  call start_clock('  proj_load')
#if 1
  call get_buffer(pp%projs, size(pp%projs), pp%projs_unit, (q-1)/nproc_pool+1)
#else
  call get_buffer(pp%projs, size(pp%projs), pp%projs_unit, q)
#endif
  call stop_clock('  proj_load')

  return

end subroutine load_projs

