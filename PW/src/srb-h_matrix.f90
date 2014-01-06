!
!> Build an matrix representation of the Hamiltonian using coefficents and a qpoint
!
SUBROUTINE build_h_matrix(ham, qpoint, pp, spin, Hk)
  USE kinds,   ONLY : DP
  USE srb_types, ONLY : basis, ham_expansion, pseudop, kproblem
  use srb_matrix, only : add_diag, block_outer
  use mp_global, only : nproc_pot
  use constants, only : rytoev
  use uspp, only : deeq, nkb
  use uspp_param, only : nh, nhm
  use ions_base, only : nat, ityp, nsp
  use cell_base, only : tpiba, bg
  use scalapack_mod, only : scalapack_localindex
  use buffers, only : open_buffer, save_buffer, get_buffer

  IMPLICIT NONE

  ! arguments
  TYPE(ham_expansion),      INTENT(in)  :: ham
  REAL(DP),                 INTENT(in)  :: qpoint(3)
  TYPE(pseudop), INTENT(INOUT)  :: pp
  integer, intent(in)  :: spin
  type(kproblem), INTENT(INOUT) :: Hk

  ! locals
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
  REAL(DP) :: tqvec(3), sqvec(3)
  COMPLEX(DP) :: zqvec(3), zqvec2
  COMPLEX(DP), allocatable :: V_half(:,:)
  real(DP), allocatable :: rwork(:)
  integer :: i, j, i_l, j_l, a, t, ioff
  logical :: islocal, info


  ! ============================================================
  ! Allocations and setup 
  ! ============================================================
  ! ===========================================================================
  ! Transform qvec
  ! ===========================================================================
  sqvec = qpoint - floor( qpoint ) ! bring to [0,1)^3
  tqvec = matmul(bg, sqvec) * tpiba ! redefine consistent cartesian coordinates
  zqvec = cmplx(tqvec, kind=DP) ! caste up to complex
  zqvec2 = cmplx(dot_product(tqvec, tqvec), kind = DP ) !likewise with the square

  Hk%H = ham%con(:,:,spin)
  call zgemv('T', 3, size(ham%con(:,:,1)), &
             cmplx(2.d0, kind=DP),  ham%lin, 3, &
                   zqvec,   1, &
             one, Hk%H, 1)
  call add_diag(ham%length, zqvec2, Hk%H, Hk%desc)

  ! ===========================================================================
  ! Add non-local component  V^{NL} = \sum <\beta|D|\beta>
  ! ===========================================================================
  if (size(pp%projs) == 1) then
    return
  endif

  if (pp%us) then
    allocate(V_half(ham%length, nkb))
    allocate(rwork(2*nhm*ham%length))
    ioff = 1 !index offset
    do t = 1, nsp
     do a = 1, nat, nproc_pot
       if (ityp(a) /= t) cycle
       ! Do the left side of the transformation
       call zlacrm(ham%length, nh(t), &
                   pp%projs(:,ioff:ioff+nh(t)), ham%length, &
                   deeq(:,:,a,spin), nhm, &
                   V_half(:, ioff:ioff+nh(t)), ham%length, &
                   rwork)
       ioff = ioff + nh(t)
     enddo
    enddo

    ! Do the right side of the transformation, summing into S_matrix
    call block_outer(ham%length, pp%nkb_l, &
                     one,  V_half, ham%length, &
                           pp%projs, ham%length, &
                     one, Hk%H, Hk%desc)
  deallocate(V_half)
  deallocate(rwork)

  else

  ioff = 1
  do t = 1, nsp
   do a = 1, nat
    if (ityp(a) /= t) cycle
    do i = 1, nh(t)
      call zherk('U', 'N', ham%length, 1, cmplx(1./deeq(i,i,a,spin)), &
                 pp%projs(:,ioff+i-1), ham%length, one, &
                 Hk%H, ham%length)
    enddo
    ioff = ioff + nh(t)
   enddo
  enddo
  endif


END SUBROUTINE build_h_matrix

