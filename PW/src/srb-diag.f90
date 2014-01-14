!
! Author: Max Hutchinson
! 6/12
!
#define DLEN 9
#define ETHR_FAC 1.D-2

recursive SUBROUTINE diagonalize (Hk, evals, evecs, num_opt, meth_opt, P, Pinv, btype_opt)
  ! Diagonalizes a Hermitian matrix, providing eigenvalues and eigenvectors
  USE kinds, ONLY: DP
  use srb_types, only : kproblem
  use srb_matrix, only : dmat, print_dmat
  USE constants,        ONLY : pi

  USE scalapack_mod, only : scalapack_diag, nprow, npcol, ctx_sq, desc_sq
  use control_flags, only : ethr
  USE mp_global,        ONLY : intra_bgrp_comm, intra_pool_comm
  USE mp,               ONLY : mp_sum, mp_barrier

  IMPLICIT NONE

  ! arguments
  type(kproblem), intent(INOUT) :: Hk
  type(dmat), intent(INOUT) :: evecs
  REAL(DP),    intent(inOUT)   :: evals(:)
  integer,     intent(in), optional :: num_opt !>!< number of eigenvectors to compute
  integer,     intent(in), optional :: meth_opt
  complex(DP), intent(in), optional :: P(:,:)
  complex(DP), intent(in), optional :: Pinv(:,:)
  integer, intent(in), optional :: btype_opt(:)
 
  COMPLEX(DP), allocatable :: z(:,:), work(:)
  COMPLEX(DP), allocatable :: psi(:,:), hpsi(:,:), spsi(:,:)
  COMPLEX(DP), allocatable :: sc(:,:), hc(:,:), vc(:,:)
  real(DP), allocatable :: ew(:)
  logical, allocatable :: conv(:)
  COMPLEX(DP), ALLOCATABLE :: cg(:), scg(:), g0(:), lagrange(:)
  real(DP), allocatable :: rwork(:), gap(:)
  integer, allocatable :: iwork(:), ifail(:), iclustr(:)
  integer :: lwork, lrwork, liwork, ne_out, nv_out, ierr
  real(DP) :: abstol, ethr_empty
  integer :: num, num2, n, np, nb1, nbn, m, iter, moved, i, j, meth
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
  REAL(DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
                              cg0, e0, es(2)
  REAL(DP)                 :: theta, cost, sint, cos2t, sin2t
  real(DP) :: avg_iter
  integer :: notconv, nbase
  integer :: kter, dav_iter, maxter = 100
  integer, allocatable :: btype(:)

  real(dp),external :: DLAMCH, dznrm2
  complex(DP), external :: ZDOTC

  n = Hk%H%desc(3)
  if (present(num_opt)) then
    num = num_opt
  else
    num = n
  end if
  if (.not. present(meth_opt) .or. ethr < 1.D-11) then
    meth = 0
  else
    meth = 0 !meth_opt
  endif
  allocate(btype(num))
  if (present(btype_opt)) then
    btype(1:num) = btype_opt(1:num)
  else
    btype = 1
  endif 

  ethr_empty = MAX( ( ethr * 5.D0 ), 1.D-5 )

  select case (meth) 
  case ( 0 ) ! Direct eigen-solve
  if (Hk%generalized) then
    if (num > 0) then
      allocate(ifail(n))
      allocate(iclustr(2*Hk%H%nprow*Hk%H%npcol), gap(Hk%H%nprow*Hk%H%npcol))
      allocate(work(1), rwork(1), iwork(1))
      allocate(z(size(Hk%H%dat,1),size(Hk%H%dat,2)))
      !call print_dmat(Hk%H)
      !call print_dmat(Hk%S)
      !call print_dmat(evecs)
      abstol = ethr
      call pzhegvx(1, 'V', 'I', 'U', n, &
                   Hk%H%dat, 1, 1, Hk%H%desc, &
                   Hk%S%dat, 1, 1, Hk%S%desc, &
                   0, 0, 1, num, &
                   abstol, ne_out, nv_out, evals, &
                   -1.d0, &
                   z, 1, 1, Hk%H%desc, &
                   work, -1, rwork, -1, iwork, -1, &
                   ifail, iclustr, gap, ierr)
      lwork = work(1); deallocate(work); allocate(work(lwork))
      lrwork = rwork(1); deallocate(rwork); allocate(rwork(lrwork))
      liwork = iwork(1); deallocate(iwork); allocate(iwork(liwork))
      call pzhegvx(1, 'V', 'I', 'U', n, &
                   Hk%H%dat, 1, 1, Hk%H%desc, &
                   Hk%S%dat, 1, 1, Hk%S%desc, &
                   0, 0, 1, num, &
                   abstol, ne_out, nv_out, evals, &
                   -1.d0, &
                   z, 1, 1, Hk%H%desc, &
                   work, lwork, rwork, lrwork, iwork, liwork, &
                   ifail, iclustr, gap, ierr)
      if (ierr /= 0) write(*,*) "zhegvx error: ", ierr
      call pzgemr2d(n, num, z, 1, 1, Hk%H%desc, evecs%dat, 1, 1, evecs%desc, Hk%H%desc(2))
      write(*,*) "Worked once?"
      deallocate(z)
      deallocate(work, rwork, iwork)
      deallocate(ifail, iclustr, gap)
    else
      lwork = 2*n
      allocate(work(lwork), rwork(3*n))
      call ZHEGV(1, 'V', 'U', n, &
                 Hk%H, n, &
                 Hk%S, n, &
                 evals, work, lwork, rwork, ierr)
      if (ierr /= 0) write(*,*) "zhegv error: ", ierr
      evecs%dat = Hk%H%dat(:,1:num)
      deallocate(work)
    endif
  else
    if (num == 0) then
      lwork = 2*n
      allocate(work(lwork), rwork(3*n))
      CALL ZHEEV('N', 'U', n, Hk%H%dat, n, evals, work, lwork, rwork, ierr) 
      evecs%dat = Hk%H%dat(:,1:num)
      deallocate(work, rwork)
    else
      allocate(ifail(n))
      allocate(z(size(Hk%H%dat,1),size(Hk%H%dat,2)))
      allocate(iclustr(2*Hk%H%nprow*Hk%H%npcol), gap(Hk%H%nprow*Hk%H%npcol))
      allocate(work(1), rwork(1), iwork(1))
      abstol = ethr
      CALL pzheevx('V', 'I', 'U', n, &
                  Hk%H%dat, 1, 1, Hk%H%desc, &
                  0, 0, 1, num, &
                  abstol, ne_out, nv_out, evals, &
                  -1.d0, &
                  z, 1, 1, Hk%H%desc, &
                  work, -1, rwork, -1, iwork, -1, &
                  ifail, iclustr, gap, ierr) 
      lwork = work(1); deallocate(work); allocate(work(lwork))
      lrwork = rwork(1); deallocate(rwork); allocate(rwork(lrwork))
      liwork = iwork(1); deallocate(iwork); allocate(iwork(liwork))
      CALL pzheevx('V', 'I', 'U', n, &
                  Hk%H%dat, 1, 1, Hk%H%desc, &
                  0, 0, 1, num, &
                  abstol, ne_out, nv_out, evals, &
                  -1.d0, &
                  z, 1, 1, Hk%H%desc, &
                  work, lwork, rwork, lrwork, iwork, liwork, &
                  ifail, iclustr, gap, ierr) 
      call pzgemr2d(n, num, z, 1, 1, Hk%H%desc, evecs%dat, 1, 1, evecs%desc, Hk%H%desc(2))
      deallocate(z)
      deallocate(ifail, iclustr, gap)
      deallocate(work, rwork, iwork)
    end if
  end if
  return
#if 0
  case ( 2 ) ! Block Davidson
  num2 = n
  allocate(psi(n,num2), hpsi(n,num2), spsi(n,num2))
  allocate(sc(num2,num2), hc(num2,num2), vc(num2,num2))
  allocate(ew(num2), conv(num))
  allocate(z(n,n))

  !
  avg_iter = 0.D0
  notconv  = num
  nbase = num
  conv = .FALSE.
  !
  IF ( present(smatrix) ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,1:num) = evecs%dat(:,1:num)
  !
  ! ... hpsi contains h times the basis vectors
  !
  call zhemm('L','U', n, num, &
             one, matrix, n, &
                  psi, n, &
             zero, hpsi, n)
  !
  ! ... spsi contains s times the basis vectors
  !
  if ( present(smatrix) ) call zhemm('L','U', n, num, &
             one, smatrix, n, &
                  psi, n, &
             zero, spsi, n)
  !
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  !
  CALL ZGEMM( 'C', 'N', nbase, nbase, n, &
             one, psi, n, &
                  hpsi, n, &
             zero, hc, num2 )
  !
  CALL mp_sum( hc( :, 1:nbase ), intra_bgrp_comm )
  !
  IF ( present(smatrix) ) THEN
     !
     CALL ZGEMM( 'C', 'N', nbase, nbase, n, &
                one, psi, n, &
                     spsi, n, &
                zero, sc, num2 )
     !     
  ELSE
     !
     CALL ZGEMM( 'C', 'N', nbase, nbase, n, &
                one, psi, n, &
                     psi, n, &
                zero, sc, num2 )
     !
  END IF
  !
  CALL mp_sum( sc( :, 1:nbase ), intra_bgrp_comm )

  ! ... diagonalize the reduced hamiltonian
  CALL cdiaghg( nbase, num, hc, sc, num2, ew, vc )
  !
  evals(1:num) = ew(1:num)
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     if (nbase * 2 >= n) then
       write(*,*) "Falling back to exact diagonalization"
       if (present(smatrix)) then
       call diagonalize (matrix, evals, evecs, num, smatrix, meth_opt = 0)
       else
       call diagonalize (matrix, evals, evecs, num, meth_opt = 0)
       endif
       deallocate(psi, hpsi, spsi)
       deallocate(sc, hc, vc)
       deallocate(ew, conv)
       return
     endif
     !
     dav_iter = kter
     !
     np = 0
     !
     DO i = 1, num
        !
        IF ( .NOT. conv(i) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) vc(:,np) = vc(:,i)
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = evals(i)
           !
        END IF
        !
     END DO
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     IF ( present(smatrix) ) THEN
        !
        CALL ZGEMM( 'N', 'N', n, notconv, nbase, &
                   one, spsi, n, & 
                        vc, num2, &
                   zero, psi(1,nb1), n )
        !     
     ELSE
        !
        CALL ZGEMM( 'N', 'N', n, notconv, nbase, &
                   one, psi, n, & 
                        vc, num2, &
                   zero, psi(1,nb1), n )
        !
     END IF
     !
     DO np = 1, notconv
        !
        call zdscal(n, - ew(nbase+np), psi(:,nbase+np), 1)
        !
     END DO
     !
     CALL ZGEMM( 'N', 'N', n, notconv, nbase, &
                one, hpsi, n, &
                     vc, num2, &
                one, psi(1,nb1), n )
     !
     !
     ! ... approximate inverse iteration
     ! ... use hpsi as scratch space
     spsi(:,nb1:nb1+notconv-1) = psi(:,nb1:nb1+notconv-1) 
     do j = 1, 0
     call zhemm('L','U', n, notconv, &
                one, Pinv, n, &
                     spsi(1,nb1), n, &
                zero, hpsi(1,nb1), n)
#if 1
     do i = 0, notconv-1
       if (present(smatrix)) then
       call zhemv('U', n, &
                  cmplx(ew(nb1+i), kind=DP), smatrix, n, &
                       hpsi(1,nb1+i), 1, &
                  zero, spsi(1,nb1+i), 1)
       else
         spsi(:,nb1+i) = hpsi(:,nb1+i)
         !call zaxpy(n, cmplx(ew(nb1+i), kind=DP), hpsi(1,nb1+i), 1, psi(1,nb1+i), 1)
       endif
       call zaxpy(n, one, psi(1,nb1+i), 1, spsi(1,nb1+i), 1) 
     enddo
     enddo

     do i = 0, notconv-1
     theta = (ew(nb1+i) - evals(1))/(evals(num) - evals(1))
     z = P * (1 - theta) + Pinv * theta
     call zhemm('L','U', n, 1, &
                one, z, n, &
                     spsi(1,nb1+i), n, &
                zero, hpsi(1,nb1+i), n)
     enddo
#endif
     psi(:,nb1:nb1+notconv-1) = hpsi(:,nb1:nb1+notconv-1) 
     !write(*,*) psi(:,nb1+notconv-1)
     !CALL g_psi( npwx, npw, notconv, npol, psi(1,1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notconv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notconv
     !
     DO i = 1, notconv
        !
        nbn = nbase + i
        !
!        IF ( npol == 1 ) THEN
           !
           ew(n) = dznrm2(n, psi(1,nbn), 1)
           !
!        ELSE
           !
!           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 ) + &
!                   ddot( 2*npw, psi(1,2,nbn), 1, psi(1,2,nbn), 1 )
           !
!        END IF
        !
     END DO
     !
     CALL mp_sum( ew( 1:notconv ), intra_bgrp_comm )
     !
     DO i = 1, notconv
        !
        call zdscal(n, 1.D0/SQRT(ew(n)), psi(:,nbase+i), 1)
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     !
     call zhemm('L','U', n, notconv, &
                one, matrix, n, &
                     psi(1,nb1), n, &
                zero, hpsi(1,nb1), n)
     if ( present(smatrix) )  call zhemm('L','U', n, notconv, &
                one, smatrix, n, &
                     psi(1,nb1), n, &
                zero, spsi(1,nb1), n)
     !
     ! ... update the reduced hamiltonian
     !
     CALL ZGEMM( 'C', 'N', nbase+notconv, notconv, n, &
                one, psi, n, &
                hpsi(1,nb1), n, &
                zero, hc(1,nb1), num2 )
     !
     CALL mp_sum( hc( :, nb1:nb1+notconv-1 ), intra_bgrp_comm )
     !
     IF ( present(smatrix) ) THEN
        !
        CALL ZGEMM( 'C', 'N', nbase+notconv, notconv, n, &
                   one, psi, n, &
                   spsi(1,nb1), n, &
                   zero, sc(1,nb1), num2 )
        !     
     ELSE
        !
        CALL ZGEMM( 'C', 'N', nbase+notconv, notconv, n, &
                   one, psi, n, &
                   psi(1,nb1), n, &
                   zero, sc(1,nb1), num2 )
        !
     END IF
     !
     CALL mp_sum( sc( :, nb1:nb1+notconv-1 ), intra_bgrp_comm )
     !
     nbase = nbase + notconv
     !
     DO i = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real 
        !
        hc(i,i) = CMPLX( REAL( hc(i,i) ), 0.D0 ,kind=DP)
        sc(i,i) = CMPLX( REAL( sc(i,i) ), 0.D0 ,kind=DP)
        !
        DO m = i + 1, nbase
           !
           hc(m,i) = CONJG( hc(i,m) )
           sc(m,i) = CONJG( sc(i,m) )
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL cdiaghg( nbase, num, hc, sc, num2, ew, vc )
     !
     ! ... test for convergence
     !
!     write(*,*) btype
     WHERE( btype(1:num) == 1 )
        !
        conv(1:num) = ( ( ABS( ew(1:num) - evals(1:num) ) < ethr * ETHR_FAC) )
        !
     ELSEWHERE
        !
        conv(1:num) = ( ( ABS( ew(1:num) - evals(1:num) ) < ethr_empty * ETHR_FAC ) )
        !
     END WHERE
     !
     notconv = COUNT( .NOT. conv(:) )
     !
     evals(1:num) = ew(1:num)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notconv == 0 .OR. &
          nbase+notconv > num2 .OR. dav_iter == maxter ) THEN
        !
        CALL ZGEMM( 'N', 'N', n, num, nbase, &
                    one, psi, n, &
                         vc, num2, &
                    zero, evecs%dat, n )
        !
        IF ( notconv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notconv
           !
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:num) = evecs(:,1:num)
        !
        IF ( present(smatrix) ) THEN
           !
           CALL ZGEMM( 'N', 'N', n, num, nbase, &
                       one, spsi, n, &
                            vc, num2, &
                       zero, psi(1,num+1), n )
           !
           spsi(:,1:num) = psi(:,num+1:num+num)
           !
        END IF
        !
        CALL ZGEMM( 'N', 'N', n, num, nbase, &
                    one, hpsi, n, &
                         vc, num2, &
                    zero, psi(1,num+1), n )
        !
        hpsi(:,1:num) = psi(:,num+1:num+num)
        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = num
        !
        hc(:,1:nbase) = ZERO
        sc(:,1:nbase) = ZERO
        vc(:,1:nbase) = ZERO
        !
        DO i = 1, nbase
           !
!           hc(n,n) = REAL( evals(n) )
           hc(i,i) = CMPLX( evals(n), 0.0_DP ,kind=DP)
           !
           sc(i,i) = ONE
           vc(i,i) = ONE
           !
        END DO
        !
        !
     END IF
     !
  END DO iterate
  write(*,*) "It", kter, notconv, nbase

  deallocate(psi, hpsi, spsi)
  deallocate(sc, hc, vc)
  deallocate(ew, conv)
  deallocate(btype)
  deallocate(z)
  return
#endif
  end select 





#if 0
  !
  ! ... every eigenfunction is calculated separately
  !
  DO m = 1, num
     !
     IF ( .true. ) then !btype(m) == 1 ) THEN
        !
        ethr_m = ethr * 1.D-2
        !
     ELSE
        !
        ethr_m = ethr * 1.D-2
        !
     END IF
     !
     spsi     = zero
     scg      = zero
     hpsi     = zero
     g        = zero
     cg       = zero
     g0       = zero
     ppsi     = zero
     lagrange = zero
     !
     ! ... calculate S|psi>
     !
     CALL zhemv('U', n, one, &
                smatrix, n, &
                evecs(:,m), 1, zero,  &
                spsi, 1)
     !
     ! ... orthogonalize starting eigenfunction to those already calculated
     !
     CALL ZGEMV('C', n, m, &
                one,  evecs, n, &
                      spsi, 1, &
                zero, lagrange, 1)
     !
     CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )
     !
     psi_norm = DBLE( lagrange(m) )
     !
     DO j = 1, m - 1
        !
        evecs(:,m)  = evecs(:,m) - lagrange(j) * evecs(:,j)
        !
        psi_norm = psi_norm - &
                   ( DBLE( lagrange(j) )**2 + AIMAG( lagrange(j) )**2 )
        !
     END DO
     !
     psi_norm = SQRT( psi_norm )
     !
     evecs(:,m) = evecs(:,m) / psi_norm
     !
     ! ... calculate starting gradient (|hpsi> = H|psi>) ...
     !
     CALL zhemv('U', n, one, &
                matrix, n, &
                evecs(:,m), 1, zero,  &
                hpsi, 1)
     CALL zhemv('U', n, one, &
                smatrix, n, &
                evecs(:,m), 1, zero,  &
                spsi, 1)
     !
     ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
     !
     ! ... NB:  ddot(2*npw,a,1,b,1) = REAL( zdotc(npw,a,1,b,1) )
     !
     evals(m) = real(zdotc(n, evecs(:,m), 1, hpsi, 1)) 
     !
     CALL mp_sum( evals(m), intra_bgrp_comm )
     !
     ! ... start iteration for this band
     !
     iterate: DO iter = 1, 1000!maxter
        !
        ! ... calculate  P (PHP)|y>
        ! ... ( P = preconditioning matrix, assumed diagonal )
        !
        if (present(P) .and. present(Pinv)) then
        call zhemm('l', 'U', n, 2, one, &
                   Pinv, n, &
                   hspsi, n, zero, &
                   gppsi, n)
        call zhemv('U', n, one, &
                   Pinv, n, &
                   spsi, 1, zero, &
                   ppsi, 1)
        else
          g(:)    = hpsi(:) 
          ppsi(:) = spsi(:) 
        endif
        !
        ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
        !
        es(1) = real(zdotc(n, spsi, 1, g, 1)) 
        es(2) = real(zdotc(n, spsi, 1, ppsi, 1)) 
        !
        CALL mp_sum( es , intra_bgrp_comm )
        !
        es(1) = es(1) / es(2)
        !
        g(:) = g(:) - es(1) * ppsi(:)
        !
        ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y> ensures that 
        ! ... <g| S P^2|y> = 0
        ! ... orthogonalize to lowest eigenfunctions (already calculated)
        !
        ! ... scg is used as workspace
        !
        CALL zhemv('U', n, one, &
                   smatrix, n, &
                   g, 1, zero,  &
                   scg, 1)
        !
        CALL ZGEMV( 'C', n, ( m - 1 ), one, evecs, &
                    n, scg, 1, zero, lagrange, 1  )
        !
        CALL mp_sum( lagrange( 1:m-1 ), intra_bgrp_comm )
        !
        DO j = 1, ( m - 1 )
           !
           g(:)   = g(:)   - lagrange(j) * evecs(:,j)
           scg(:) = scg(:) - lagrange(j) * evecs(:,j)
           !
        END DO
        !
        IF ( iter /= 1 ) THEN
           !
           ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
           !
           gg1 = real(zdotc(n, g, 1, g0, 1)) 
           !
           CALL mp_sum( gg1, intra_bgrp_comm )
           !
        END IF
        !
        ! ... gg is <g(n+1)|S|g(n+1)>
        !
!        g0(:) = scg(:)
       !
!        g0(:) = g0(:) !* precondition(:)
        if (present(P) .and. present(Pinv)) then
        call zhemv('U', n, one, &
                   P, n, &
                   scg, 1, zero, &
                   g0, 1)
        else
          g0 = scg
        endif
        !
        gg = real(zdotc( n, g, 1, g0, 1)) 
        !
        CALL mp_sum( gg, intra_bgrp_comm )
        !
        IF ( iter == 1 ) THEN
           !
           ! ... starting iteration, the conjugate gradient |cg> = |g>
           !
           gg0 = gg
           !
           cg(:) = g(:)
           !
        ELSE
           !
           ! ... |cg(n+1)> = |g(n+1)> + gamma(n) * |cg(n)>
           !
           ! ... Polak-Ribiere formula :
           !
           gamma = ( gg - gg1 ) / gg0
           gg0   = gg
           !
           cg(:) = cg(:) * gamma
           cg(:) = g + cg(:)
           !
           ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)> 
           ! ... is not 0. In fact :
           ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
           !
           psi_norm = gamma * cg0 * sint
           !
           cg(:) = cg(:) - psi_norm * evecs(:,m)
           !
        END IF
        !
        ! ... |cg> contains now the conjugate gradient
        !
        ! ... |scg> is S|cg>
        !
        CALL zhemv('U', n, one, &
                   matrix, n, &
                   cg, 1, zero,  &
                   ppsi, 1)
        CALL zhemv('U', n, one, &
                   smatrix, n, &
                   cg, 1, zero,  &
                   scg, 1)
        !
        cg0 = real(zdotc(n, cg, 1, scg, 1)) 
        !
        CALL mp_sum(  cg0 , intra_bgrp_comm )
        !
        cg0 = SQRT( cg0 )
        !
        ! ... |ppsi> contains now HP|cg>
        ! ... minimize <y(t)|PHP|y(t)> , where :
        ! ...                         |y(t)> = cos(t)|y> + sin(t)/cg0 |cg>
        ! ... Note that  <y|P^2S|y> = 1, <y|P^2S|cg> = 0 ,
        ! ...           <cg|P^2S|cg> = cg0^2
        ! ... so that the result is correctly normalized :
        ! ...                           <y(t)|P^2S|y(t)> = 1
        !
        a0 = 2.D0 * real(zdotc(n, evecs(:,m), 1, ppsi, 1)) / cg0
        !
        CALL mp_sum(  a0 , intra_bgrp_comm )
        !
        b0 = real(zdotc(n, cg, 1, ppsi, 1)) / cg0**2
        !
        CALL mp_sum(  b0 , intra_bgrp_comm )
        !
        e0 = evals(m)
        !
        theta = 0.5D0 * ATAN( a0 / ( e0 - b0 ) )
        !
        cost = COS( theta )
        sint = SIN( theta )
        !
        cos2t = cost*cost - sint*sint
        sin2t = 2.D0*cost*sint
        !
        es(1) = 0.5D0 * (   ( e0 - b0 ) * cos2t + a0 * sin2t + e0 + b0 )
        es(2) = 0.5D0 * ( - ( e0 - b0 ) * cos2t - a0 * sin2t + e0 + b0 )
        !
        ! ... there are two possible solutions, choose the minimum
        !
        IF ( es(2) < es(1) ) THEN
           !
           theta = theta + 0.5D0 * pi
           !
           cost = COS( theta )
           sint = SIN( theta )
           !
        END IF
        !
        ! ... new estimate of the eigenvalue
        !
        evals(m) = MIN( es(1), es(2) )
        !
        ! ... upgrade |psi>
        !
        evecs(:,m) = cost * evecs(:,m) + sint / cg0 * cg(:)
        !
        ! ... here one could test convergence on the energy
        !
        IF ( ABS( evals(m) - e0 ) < ethr_m ) EXIT iterate
        !
        ! ... upgrade H|psi> and S|psi>
        !
        spsi(:) = cost * spsi(:) + sint / cg0 * scg(:)
        !
        hpsi(:) = cost * hpsi(:) + sint / cg0 * ppsi(:)
        !
     END DO iterate
     !
     !IF ( iter >= maxter ) notconv = notconv + 1
     !
     !write(*,*) iter
     avg_iter = avg_iter + iter + 1
     !
     ! ... reorder evals if they are not in the right order
     ! ... ( this CAN and WILL happen in not-so-special cases )
     !
     IF ( m > 1 .AND. .true. ) then !reorder ) THEN
        !
        IF ( evals(m) - evals(m-1) < - 2.D0 * ethr_m ) THEN
           !
           ! ... if the last calculated eigenvalue is not the largest...
           !
           DO i = m - 2, 1, - 1
              !
              IF ( evals(m) - evals(i) > 2.D0 * ethr_m ) EXIT
              !
           END DO
           !
           i = i + 1
           !
           moved = moved + 1
           !
           ! ... last calculated eigenvalue should be in the 
           ! ... i-th position: reorder
           !
           e0 = evals(m)
           !
           ppsi(:) = evecs(:,m)
           !
           DO j = m, i + 1, - 1
              !
              evals(j) = evals(j-1)
              !
              evecs(:,j) = evecs(:,j-1)
              !
           END DO
           !
           evals(i) = e0
           !
           evecs(:,i) = ppsi(:)
           !
           ! ... this procedure should be good if only a few inversions occur,
           ! ... extremely inefficient if eigenvectors are often in bad order
           ! ... ( but this should not happen )
           !
        END IF
        !
     END IF
     !
  END DO
  !
  avg_iter = avg_iter / DBLE( num )
  write(*,*) "avg_iter:", avg_iter
  !
  DEALLOCATE( lagrange )
  DEALLOCATE( gppsi )
  DEALLOCATE( g0 )
  DEALLOCATE( cg )
  DEALLOCATE( hspsi )
  DEALLOCATE( scg )
#endif

END SUBROUTINE diagonalize


