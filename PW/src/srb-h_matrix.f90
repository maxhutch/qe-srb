!
!> Build an matrix representation of the Hamiltonian using coefficents and a qpoint
!
SUBROUTINE build_h_matrix(ham, qpoint, pp, spin, Hk, q)
  USE kinds,   ONLY : DP
  USE srb_types, ONLY : basis, ham_expansion, pseudop, kproblem
  use srb_matrix, only : add_diag, block_outer, print_dmat
  use mp_global, only : nproc_pot
  use constants, only : rytoev
  use uspp, only : deeq, nkb
  use uspp_param, only : nh, nhm
  use ions_base, only : nat, ityp, nsp
  use cell_base, only : tpiba, bg
  use scalapack_mod, only : scalapack_localindex
  use buffers, only : open_buffer, save_buffer, get_buffer, close_buffer

  IMPLICIT NONE

  ! arguments
  TYPE(ham_expansion),      INTENT(in)  :: ham
  REAL(DP),                 INTENT(in)  :: qpoint(3)
  TYPE(pseudop), INTENT(INOUT)  :: pp
  integer, intent(in)  :: spin
  type(kproblem), INTENT(INOUT) :: Hk
  integer :: q

  ! locals
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
  REAL(DP) :: tqvec(3), sqvec(3)
  COMPLEX(DP) :: zqvec(3), zqvec2
  COMPLEX(DP), allocatable :: V_half(:,:)
  real(DP), allocatable :: rwork(:)
  integer :: i, j, i_l, j_l, a, t, ioff, prow, pcol
  integer :: nbasis
  integer, save :: old_size_h_matrix
  logical :: islocal, info
  complex(DP) :: trace 

  ! ============================================================
  ! Allocations and setup 
  ! ============================================================
  nbasis = ham%con(1)%desc(3)
  ! ===========================================================================
  ! Compute non-local component  V^{NL} = \sum <\beta|D|\beta>
  ! ===========================================================================

  ! if no projectors, start with 0
  if (.not. allocated(pp%projs)) then
    Hk%H%dat = cmplx(0.d0, kind=DP)
  ! if NCPP, we might be able to load V^nl
  else if (.not. pp%us .and. q > 0) then
    call get_buffer(Hk%H%dat, size(Hk%H%dat), pp%h_unit, q)
  ! if USPP or V^nl has changed, compute it
  else 
    do t = 1, pp%ntyp
      allocate(V_half(nbasis, size(pp%projs(t)%dat,2)))
      allocate(rwork(2*nhm*nbasis))
      ioff = 1 !index offset
      do a = 1, pp%na(t)
        call infog2l(1, 1+(a-1)*nh(t), &
                     pp%projs(t)%desc, pp%projs(t)%nprow, pp%projs(t)%npcol, &
                                       pp%projs(t)%myrow, pp%projs(t)%mycol, &
                     i_l, j_l, prow, pcol)
        if (prow == pp%projs(t)%myrow .and. pcol == pp%projs(t)%mycol) then
          ! Do the left side of the transformation
          call zlacrm(nbasis, nh(t), &
                      pp%projs(t)%dat(:,ioff:ioff+nh(t)), nbasis, &
                      deeq(:,:,a+pp%nt_off(t)-1,spin), nhm, &
                      V_half(:, ioff:ioff+nh(t)), nbasis, &
                      rwork)
          ioff = ioff + nh(t)
        endif
      enddo

      ! Do the right side of the transformation, summing into S_matrix
      call block_outer(nbasis, size(pp%projs(t)%dat,2), &
                       one,  V_half, nbasis, &
                             pp%projs(t)%dat, nbasis, &
                       zero, Hk%H)
      deallocate(V_half)
      deallocate(rwork)
    enddo

    ! if NCPP, save it
    if (.not. pp%us) then
      ! open a buffer if we don't have one
      if (pp%h_unit < 0) then
        pp%h_unit = - pp%h_unit
        call open_buffer(pp%h_unit, 'h_matrix', size(Hk%H%dat), 1, info)
        old_size_h_matrix=size(Hk%H%dat)
      endif
      ! if the buffer increases in size - dynamically change it
      if (size(Hk%H%dat) > old_size_h_matrix) then
        call close_buffer(pp%h_unit,'delete')
        call open_buffer(pp%h_unit, 'h_matrix', size(Hk%H%dat), 1, info)
        old_size_h_matrix=size(Hk%H%dat)
      endif
      !save it
      call save_buffer(Hk%H%dat, size(Hk%H%dat), pp%h_unit, -q)
     endif
  endif

  ! ===========================================================================
  ! Transform qvec for adding kinetic energy
  ! ===========================================================================
  sqvec = qpoint - floor( qpoint ) ! bring to [0,1)^3
  tqvec = matmul(bg, sqvec) * tpiba ! redefine consistent cartesian coordinates
  zqvec = cmplx(tqvec, kind=DP) ! caste up to complex
  zqvec2 = cmplx(dot_product(tqvec, tqvec), kind = DP ) !likewise with the square

  ! ===========================================================================
  ! Add kinetic energy and local potential 
  ! ===========================================================================
  Hk%H%dat = Hk%H%dat + ham%con(spin)%dat
  call zgemv('T', 3, size(ham%con(1)%dat), &
             cmplx(2.d0, kind=DP),  ham%lin, 3, &
                   zqvec,   1, &
             one, Hk%H%dat, 1)
  call add_diag(Hk%H, zqvec2)


!  else
#if 0
  ioff = 1
  do t = 1, nsp
   do a = 1, nat
    if (ityp(a) /= t) cycle
    do i = 1, nh(t)
      call zherk('U', 'N', nbasis, 1, cmplx(deeq(i,i,a,spin)), &
                 pp%projs(t)%dat(:,ioff+i-1), nbasis, one, &
                 Hk%H%dat, nbasis)
    enddo
    ioff = ioff + nh(t)
   enddo
  enddo
#endif
!  endif

END SUBROUTINE build_h_matrix

