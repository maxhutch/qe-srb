!
!> Build an matrix representation of the Hamiltonian using coefficents and a qpoint
!
SUBROUTINE build_h_matrix(ham, qpoint, pp, spin, ham_matrix)
  USE kinds,   ONLY : DP
  USE srb_types, ONLY : basis, ham_expansion, pseudop
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
  COMPLEX(DP), INTENT(OUT) :: ham_matrix(:,:)

  ! locals
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
  REAL(DP) :: tqvec(3), sqvec(3)
  COMPLEX(DP) :: zqvec(3), zqvec2
  COMPLEX(DP), allocatable :: V_nl(:,:), V_half(:,:)
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
#define __SSDIAG
#ifdef __SSDIAG
  forall(i = 1 : ham%length, j = 1 : ham%length)
    ! constant term
    ham_matrix(i, j) = ham%con(i,j, spin)
    ! linear term
    ham_matrix(i, j) = ham_matrix(i, j) + 2.d0 * sum(ham%lin(:,i,j) * zqvec)
!    ham_matrix(i, j) = 2.d0 * sum(ham%lin(:,i,j) * zqvec)
  end forall
  ! quadratic term
  forall(i = 1 : ham%length) ham_matrix(i, i) = ham_matrix(i, i) + zqvec2
#else
  do j = 1, ham%length
    do i = 1, ham%length
      call scalapack_localindex(i, j, i_l, j_l, islocal)
      if (islocal) then
        ! constant term
        ham_matrix(i_l, j_l) = ham%con(i,j, spin)
        ! linear term
        ham_matrix(i_l, j_l) = ham_matrix(i_l, j_l) + 2.d0 * sum(ham%lin(:,i,j) * zqvec)
        ! quadratic term
        if (i == j) ham_matrix(i_l, j_l) = ham_matrix(i_l, j_l) + zqvec2
      end if 
    enddo
  enddo
#endif

  ! ===========================================================================
  ! Add non-local component  V^{NL} = \sum <\beta|D|\beta>
  ! ===========================================================================
  if (size(pp%projs) == 1) then
    return
  endif

  if (pp%us) then
  allocate(V_nl(ham%length, ham%length), V_half(ham%length, nkb))
  allocate(rwork(2*nhm*ham%length))
  ioff = 1 !index offset
  do t = 1, nsp
   do a = 1, nat
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
  call ZGEMM('N', 'C', ham%length, ham%length, nkb, one, &
               V_half(:, :), ham%length, &
               pp%projs(:, :), ham%length, zero, &
               V_nl(:,:), ham%length)
#ifdef __SSDIAG
  ham_matrix = ham_matrix + V_nl
#else
  do i = 1, ham%length; do j = 1, ham%length
    call scalapack_localindex(i, j, i_l, j_l, islocal)
    if (islocal) ham_matrix(i_l, j_l) = ham_matrix(i_l, j_l) + V_nl(i,j)
  enddo; enddo
#endif
  deallocate(V_nl, V_half)
  deallocate(rwork)
  else
  ioff = 1
  do t = 1, nsp
   do a = 1, nat
    if (ityp(a) /= t) cycle
    do i = 1, nh(t)
      call zherk('U', 'N', ham%length, 1, cmplx(1./deeq(i,i,a,spin)), &
                 pp%projs(:,ioff+i-1), ham%length, one, &
                 ham_matrix, ham%length)
    enddo
    ioff = ioff + nh(t)
   enddo
  enddo
  endif


END SUBROUTINE build_h_matrix

