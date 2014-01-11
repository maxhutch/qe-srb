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
  integer :: a, t, ioff, joff
  integer :: i, j, i_l, j_l
  logical :: islocal
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP)
  complex(DP), parameter :: one = cmplx(1.d0, kind=DP)
  integer :: nbasis
  logical :: info, prow, pcol
  complex(DP) :: trace

  integer,save :: old_size_s_matrix

  nbasis = size(pp%projs(1)%dat, 1)

  ! if the save file hasn't been opened, open one
  if (pp%s_unit < 0) then
    pp%s_unit = - pp%s_unit
    call open_buffer(pp%s_unit, 's_matrix', size(Hk%S%dat), 1, info)
    old_size_s_matrix=size(Hk%S%dat)
  endif

  ! if the buffer increases in size - dynamically change it
  if (size(Hk%S%dat) > old_size_s_matrix) then
    call close_buffer(pp%s_unit,'delete')
    call open_buffer(pp%s_unit, 's_matrix', size(Hk%S%dat), 1, info)
    old_size_s_matrix=size(Hk%S%dat)
  endif

  ! if saved, just reload it
  if (q > 0) then
    call get_buffer(Hk%S%dat, size(Hk%S%dat), pp%s_unit, q)
    return
  endif

  ! Make <B|S|B>
  Hk%S%dat = 0.d0
  do t = 1, pp%ntyp
    allocate(S_half(nbasis, size(pp%projs(t)%dat,2)))
    allocate(rwork(2*nhm*nbasis))
    ioff = 1 !index offset
    do a = 1, pp%na(t)
      call infog2l(1+(a-1)*nh(t), 1, &
                   pp%projs(t)%desc, pp%projs(t)%nprow, pp%projs(t)%npcol, &
                                     pp%projs(t)%myrow, pp%projs(t)%mycol, &
                   i_l, j_l, prow, pcol)
      if (prow == pp%projs(t)%myrow .and. pcol == pp%projs(t)%mycol) then
        ! Do the left side of the transformation
        call zlacrm(nbasis, nh(t), &
                    pp%projs(t)%dat(:,ioff:ioff+nh(t)), nbasis, &
                    qq(:,:,t), nhm, &
                    S_half(:, ioff:ioff+nh(t)), nbasis, &
                    rwork)
        ioff = ioff + nh(t)
      endif
    enddo

    ! Do the right side of the transformation, summing into S_matrix
    call block_outer(nbasis, size(pp%projs(t)%dat,2), &
                     one,  S_half, nbasis, &
                           pp%projs(t)%dat, nbasis, &
                     one, Hk%S)
    deallocate(S_half)
    deallocate(rwork)
  enddo

  call add_diag(Hk%S, one)

#if 0
  trace = 0.d0
  do a = 1, nbasis
    trace = trace + Hk%S%dat(a,a)
  enddo
  write(*,*) "Trace[S] = ", trace
#endif

  call save_buffer(Hk%S%dat, size(Hk%S%dat), pp%s_unit, -q)

  return

end subroutine build_s_matrix

