!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this module contains routines for fft with an custom selected cutoff

MODULE fft_custom_gwl

  USE kinds, ONLY: DP
  USE parallel_include
  
  USE fft_types, ONLY: fft_dlay_descriptor
  
  IMPLICIT NONE

  TYPE fft_cus
  
        ! ... data structure containing all information
        ! ... about fft data distribution for a given
        ! ... potential grid, and its wave functions sub-grid.

     TYPE ( fft_dlay_descriptor ) :: dfftt ! descriptor for custom grim
  
     REAL(kind=DP) :: ecutt!custom cutoff in rydberg
     REAL(kind=DP) :: dual_t!dual facor
     REAL(kind=DP) :: gcutmt
     INTEGER :: nr1t,nr2t,nr3t
     INTEGER :: nrx1t,nrx2t,nrx3t
     INTEGER :: nrxxt
     INTEGER :: ngmt,ngmt_l,ngmt_g
     INTEGER, DIMENSION(:), POINTER :: nlt,nltm
     REAL(kind=DP), DIMENSION(:), POINTER :: ggt
     REAL(kind=DP), DIMENSION(:,:),POINTER :: gt
     INTEGER, DIMENSION(:), POINTER :: ig_l2gt 
     INTEGER :: gstart_t
     INTEGER,  DIMENSION(:), POINTER :: ig1t,ig2t,ig3t
     INTEGER :: nlgt
     INTEGER :: npwt,npwxt


!we redifine the cell for arbitrary cell possibility

     REAL(DP) :: alat_t = 0.0_DP
     REAl(DP) :: omega_t = 0.0_DP
     REAL(DP) :: tpiba_t  = 0.0_DP, tpiba2_t = 0.0_DP
     REAL(DP) :: at_t(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )
     REAL(DP) :: bg_t(3,3) = RESHAPE( (/ 0.0_DP /), (/ 3, 3 /), (/ 0.0_DP /) )

  
     
  END TYPE fft_cus


!=-------------------------------------------------------

!=----------------------------------------------------------------------=!
CONTAINS
!=----------------------------------------------------------------------=!
!
  
  SUBROUTINE set_custom_grid(fc)

  !-----------------------------------------------------------------------
  !     This routine computes the dimensions of the minimum FFT grid
  !     compatible with the input cut-off
  !
  !     NB: The values of nr1, nr2, nr3 are computed only if they are not
  !     given as input parameters. Input values are kept otherwise.
  !
  USE io_global,  ONLY : stdout
  use fft_scalar, only: allowed
  implicit none

  type(fft_cus) :: fc

  integer, parameter :: nmax = 5000
  ! an unreasonably big number for a FFT grid
  !
  ! the values of nr1, nr2, nr3 are computed only if they are not given
  ! as input parameters
  !
  fc%gcutmt = fc%dual_t*fc%ecutt / fc%tpiba2_t
  fc%nr1t=0
  fc%nr2t=0
  fc%nr3t=0

  if (fc%nr1t == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     fc%nr1t = int (2 * sqrt (fc%gcutmt) * sqrt (fc%at_t (1, 1) **2 + fc%at_t (2, 1) &
          **2 + fc%at_t (3, 1) **2) ) + 1
10   continue
     if (fc%nr1t > nmax) &
          call errore ('set_fft_dim', 'nr1 is unreasonably large', fc%nr1t)
     if (allowed (fc%nr1t) ) goto 15
     fc%nr1t = fc%nr1t + 1
     goto 10
  else
     if (.not.allowed (fc%nr1t) ) call errore ('set_fft_dim', &
          'input nr1t value not allowed', 1)
  endif
15 continue
   !
  if (fc%nr2t == 0) then
     !
     ! estimate nr1 and check if it is an allowed value for FFT
     !
     fc%nr2t = int (2 * sqrt (fc%gcutmt) * sqrt (fc%at_t (1, 2) **2 + fc%at_t (2, 2) &
          **2 + fc%at_t (3, 2) **2) ) + 1
20   continue
     if (fc%nr2t > nmax) &
          call errore ('set_fft_dim', 'nr2t is unreasonably large', fc%nr2t)
     if (allowed (fc%nr2t) ) goto 25
     fc%nr2t = fc%nr2t + 1
     goto 20
  else
     if (.not.allowed (fc%nr2t) ) call errore ('set_fft_dim', &
          'input nr2t value not allowed', 2)
  endif
25 continue
  !
  if (fc%nr3t == 0) then
     !
     ! estimate nr3 and check if it is an allowed value for FFT
     !
     fc%nr3t = int (2 * sqrt (fc%gcutmt) * sqrt (fc%at_t (1, 3) **2 + fc%at_t (2, 3) &
          **2 + fc%at_t (3, 3) **2) ) + 1
30   continue
     if (fc%nr3t > nmax) &
          call errore ('set_fft_dim', 'nr3 is unreasonably large', fc%nr3t)
     if (allowed (fc%nr3t) ) goto 35
     fc%nr3t = fc%nr3t + 1
     goto 30
  else
     if (.not.allowed (fc%nr3t) ) call errore ('set_fft_dim', &
          'input nr3t value not allowed', 3)
  endif
35 continue
  !
  !    here we compute nr3s if it is not in input
  !
  
  !
  
    return
  END SUBROUTINE set_custom_grid


  SUBROUTINE data_structure_custom(fc)
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the fft arrays
  ! (both the smooth and the hard mesh)
  ! In the parallel case, it distributes columns to processes, too
  !
  USE io_global,  ONLY : stdout
  USE kinds,      ONLY : DP
  USE klist,      ONLY : xk, nks
!  USE gvect,      ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
!                         ngm, ngm_l, ngm_g, gcutm, ecutwfc
!  USE gsmooth,    ONLY : nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs, &
!                         ngms, ngms_l, ngms_g, gcutms
  USE mp,         ONLY : mp_sum, mp_max,mp_barrier
  USE mp_global,  ONLY : intra_pool_comm, nproc_pool, me_pool, my_image_id, &
                         inter_pool_comm,root_pool
  USE mp_world,   ONLY : world_comm, nproc
  USE stick_base
  USE fft_scalar, ONLY : good_fft_dimension
  USE fft_types,  ONLY : fft_dlay_allocate, fft_dlay_set, fft_dlay_scalar
  !
  !
  IMPLICIT NONE

  TYPE(fft_cus) :: fc

  INTEGER :: n1, n2, n3, i1, i2, i3
  ! counters on G space
  !

  REAL(DP) :: amod
  ! modulus of G vectors

  INTEGER, ALLOCATABLE :: stw(:,:)
  ! sticks maps

  INTEGER :: ub(3), lb(3)
  ! upper and lower bounds for maps

  REAL(DP) :: gkcut
  ! cut-off for the wavefunctions

  INTEGER  :: ncplane, nxx
  INTEGER  :: ncplanes, nxxs

#ifdef __PARA
  INTEGER, ALLOCATABLE :: st(:,:), sts(:,:)
  ! sticks maps

  INTEGER, ALLOCATABLE :: ngc (:), ngcs (:), ngkc (:)
  INTEGER  ::  ncp (nproc), nct, nkcp (nproc), ncts, ncps(nproc)
  INTEGER  ::  ngp (nproc), ngps(nproc), ngkp (nproc), ncp_(nproc),&
       i, j, jj, idum

  !      nxx              !  local fft data dim
  !      ncplane,        &!  number of columns in a plane
  !      nct,            &!  total number of non-zero columns
  !      ncp(nproc),     &!  number of (density) columns per proc


  LOGICAL :: tk = .true.
  ! map type: true for full space sticks map, false for half space sticks map
  INTEGER, ALLOCATABLE :: in1(:), in2(:), idx(:)
  ! sticks coordinates

  !
  !  Subroutine body
  !

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT1.0'
  call flush_unit(stdout)
  
  
  tk = .false.

  !
  ! set the values of fft arrays
  !
  fc%nrx1t  = good_fft_dimension (fc%nr1t)
  fc%nrx2t  = fc%nr2t          ! nrx2 is there just for compatibility
  fc%nrx3t  = good_fft_dimension (fc%nr3t)


  ! compute number of points per plane
  ncplane  = fc%nrx1t * fc%nrx2t
  ncplanes = fc%nrx1t * fc%nrx2t
  call mp_barrier( world_comm )
  write(stdout,*) 'ATT1.1'
  call flush_unit(stdout)

  !
  ! check the number of plane per process
  !
  IF ( fc%nr3t < nproc_pool ) &
    CALL infomsg ('data_structure', 'some processors have no planes ')


  !
  ! compute gkcut calling an internal procedure
  !
  CALL calculate_gkcut()
  !
  !     Now compute for each point of the big plane how many column have
  !     non zero vectors on the smooth and thick mesh
  !
  call mp_barrier( world_comm )
  write(stdout,*) 'ATT1.2'
  call flush_unit(stdout)

  n1 = fc%nr1t + 1
  n2 = fc%nr2t + 1
  n3 = fc%nr3t + 1

  ub =  (/  n1,  n2,  n3 /)
  lb =  (/ -n1, -n2, -n3 /)

  ALLOCATE( stw ( lb(1) : ub(1), lb(2) : ub(2) ) )
  ALLOCATE( st  ( lb(1) : ub(1), lb(2) : ub(2) ) )
  ALLOCATE( sts ( lb(1) : ub(1), lb(2) : ub(2) ) )

 !
! ...     Fill in the stick maps, for given g-space base (b1,b2,b3)
! ...     and cut-offs
! ...     The value of the element (i,j) of the map ( st ) is equal to the
! ...     number of G-vector belonging to the (i,j) stick.
!
  call mp_barrier( world_comm )
  write(stdout,*) 'ATT1.3'
  call flush_unit(stdout)

  CALL sticks_maps( tk, ub, lb, fc%bg_t(:,1), fc%bg_t(:,2), fc%bg_t(:,3), fc%gcutmt, gkcut, fc%gcutmt, st, &
       &stw, sts ,me_pool,nproc_pool,intra_pool_comm)
  write(stdout,*) 'ATT1.3.1'
  call flush_unit(stdout)

  
  nct  = count( st  > 0 )
  ncts = count( sts > 0 )

  IF ( nct > ncplane )    &
     &    CALL errore('data_structure','too many sticks',1)

  IF ( ncts > ncplanes )  &
     &    CALL errore('data_structure','too many sticks',2)

  IF ( nct  == 0 ) &
     &    CALL errore('data_structure','number of sticks 0', 1)

  IF ( ncts == 0 ) &
     &    CALL errore('data_structure','number smooth sticks 0', 1)

    
  !
  !   local pointers deallocated at the end
  !
  ALLOCATE( in1( nct ), in2( nct ) )
  ALLOCATE( ngc( nct ), ngcs( nct ), ngkc( nct ) )
  ALLOCATE( idx( nct ) )

!
! ...     initialize the sticks indexes array ist
! ...     nct counts columns containing G-vectors for the dense grid
! ...     ncts counts columns contaning G-vectors for the smooth grid
!
  write(stdout,*) 'ATT1.5'
  call flush_unit(stdout)


  CALL sticks_countg( tk, ub, lb, st, stw, sts, in1, in2, ngc, ngkc, ngcs )

  CALL sticks_sort( ngc, ngkc, ngcs, nct, idx ,nproc_pool)

  CALL sticks_dist( tk, ub, lb, idx, in1, in2, ngc, ngkc, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngkp, ngps, st, stw, sts ,nproc_pool)

  CALL sticks_pairup( tk, ub, lb, idx, in1, in2, ngc, ngkc, ngcs, nct, &
          ncp, nkcp, ncps, ngp, ngkp, ngps, st, stw, sts ,nproc_pool)

  write(stdout,*) 'ATT1.6'
  call flush_unit(stdout)

  !  set the total number of G vectors

  IF( tk ) THEN
    fc%ngmt  = ngp ( me_pool + 1 )
 ELSE
    IF( st( 0, 0 ) == ( me_pool + 1 ) ) THEN
       fc%ngmt  = ngp ( me_pool + 1 ) / 2 + 1
    ELSE
       fc%ngmt  = ngp ( me_pool + 1 ) / 2
    ENDIF
  ENDIF

  CALL fft_dlay_allocate( fc%dfftt, me_pool,root_pool,nproc_pool,intra_pool_comm ,0, fc%nrx1t,  fc%nrx2t  )
 

  !  here set the fft data layout structures for dense and smooth mesh,
  !  according to stick distribution

  CALL fft_dlay_set( fc%dfftt, &
       tk, nct, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, ub, lb, idx, in1(:), in2(:), ncp, nkcp, ngp, ngkp, st, stw)


  !  if tk = .FALSE. only half reciprocal space is considered, then we
  !  need to correct the number of sticks

  IF( .not. tk ) THEN
    nct  = nct*2  - 1
    ncts = ncts*2 - 1
  ENDIF

  !
  ! set the number of plane per process
  !

  ! npp ( 1 : nproc_pool ) = dfftp%npp ( 1 : nproc_pool )
  ! npps( 1 : nproc_pool ) = dffts%npp ( 1 : nproc_pool )

  IF ( fc%dfftt%nnp /= ncplane )    &
     &    CALL errore('data_structure','inconsistent plane dimension on dense grid', abs(fc%dfftt%nnp-ncplane) )

  
  WRITE( stdout, '(/5x,"Planes per process (custom) : nr3t =", &
       &        i4," npp = ",i4," ncplane =",i6)') fc%nr3t, fc%dfftt%npp (me_pool + 1) , ncplane



  WRITE( stdout,*)
  WRITE( stdout,'(5X,                                                     &
    & "Proc/  planes cols     G    "/5X)')
  DO i=1,nproc_pool
    WRITE( stdout,'(5x,i4,1x,i5,i7,i9)') i, fc%dfftt%npp(i), ncp(i), ngp(i)
  ENDDO
  IF ( nproc_pool > 1 ) THEN
      WRITE( stdout,'(5x,"tot",2x,i5,i7,i9)') &
      sum(fc%dfftt%npp(1:nproc_pool)), sum(ncp(1:nproc_pool)), sum(ngp(1:nproc_pool))
  ENDIF
  WRITE( stdout,*)

  DEALLOCATE( stw, st, sts, in1, in2, idx, ngc, ngcs, ngkc )

  !
  !   ncp0 = starting column for each processor
  !

  ! ncp0( 1:nproc_pool )  = dfftp%iss( 1:nproc_pool )
  ! ncp0s( 1:nproc_pool ) = dffts%iss( 1:nproc_pool )

  !
  !  array ipc and ipcl ( ipc contain the number of the
  !                       column for that processor or zero if the
  !                       column do not belong to the processor,
  !                       ipcl contains the point in the plane for
  !                       each column)
  !
  !  ipc ( 1:ncplane )    = >  dfftp%isind( 1:ncplane )
  !  icpl( 1:nct )        = >  dfftp%ismap( 1:nct )

  !  ipcs ( 1:ncplanes )  = >  dffts%isind( 1:ncplanes )
  !  icpls( 1:ncts )      = >  dffts%ismap( 1:ncts )

  fc%nrxxt  = fc%dfftt%nnr

  !
  ! nxx is just a copy
  !
  nxx   = fc%nrxxt
  nxxs  = fc%nrxxt


#else

  fc%nrx1t = good_fft_dimension (fc%nr1t)

  !
  !     nrx2 and nrx3 are there just for compatibility
  !
  fc%nrx2t = fc%nr2t
  fc%nrx3t = fc%nr3t

  fc%nrxxt = fc%nrx1t * fc%nrx2t * fc%nrx3t
  

  ! nxx is just a copy
  !
  nxx   = fc%nrxxt
  nxxs  = fc%nrxxt

  CALL fft_dlay_allocate( fc%dfftt, me_pool,root_pool,nproc_pool, intra_pool_comm,0,max(fc%nrx1t, fc%nrx3t),  fc%nrx2t  )
 

  CALL calculate_gkcut()

  !
  !     compute the number of g necessary to the calculation
  !
  n1 = fc%nr1t + 1
  n2 = fc%nr2t + 1
  n3 = fc%nr3t + 1

  fc%ngmt = 0
 

  ub =  (/  n1,  n2,  n3 /)
  lb =  (/ -n1, -n2, -n3 /)
!
  ALLOCATE( stw ( lb(2):ub(2), lb(3):ub(3) ) )
  stw = 0

  DO i1 = - n1, n1
     !
     ! Gamma-only: exclude space with x<0
     !
     IF (i1 < 0) GOTO 10
     !
     DO i2 = - n2, n2
        !
        ! Gamma-only: exclude plane with x=0, y<0
        !
        IF(i1 == 0.and. i2 < 0) GOTO 20
        !
        DO i3 = - n3, n3
           !
           ! Gamma-only: exclude line with x=0, y=0, z<0
           !
           IF(i1 == 0 .and. i2 == 0 .and. i3 < 0) GOTO 30
           !
           amod = (i1 * fc%bg_t (1, 1) + i2 * fc%bg_t (1, 2) + i3 * fc%bg_t (1, 3) ) **2 + &
                  (i1 * fc%bg_t (2, 1) + i2 * fc%bg_t (2, 2) + i3 * fc%bg_t (2, 3) ) **2 + &
                  (i1 * fc%bg_t (3, 1) + i2 * fc%bg_t (3, 2) + i3 * fc%bg_t (3, 3) ) **2
           IF (amod <= fc%gcutmt)  fc%ngmt  = fc%ngmt  + 1
           IF (amod <= fc%gcutmt ) THEN
              stw( i2, i3 ) = 1
              stw( -i2, -i3 ) = 1
           ENDIF
30         CONTINUE
        ENDDO
20      CONTINUE
     ENDDO
10   CONTINUE
  ENDDO

  CALL fft_dlay_scalar( fc%dfftt, ub, lb, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, stw )


  DEALLOCATE( stw )

#endif

  IF( nxx < fc%dfftt%nnr ) &
     CALL errore( ' data_structure ', ' inconsistent value for nxx ', abs( nxx - fc%dfftt%nnr ) )

  !
  !     compute the global number of g, i.e. the sum over all processors
  !     within a pool
  !
  fc%ngmt_l  = fc%ngmt
  fc%ngmt_g  = fc%ngmt  ; CALL mp_sum( fc%ngmt_g , intra_pool_comm )


!  IF( use_task_groups ) THEN
  IF(.false.) THEN !ATTENZIONE
    !
    !  Initialize task groups.
    !  Note that this call modify dffts adding task group data.
    !
    CALL task_groups_init( fc%dfftt )
    !
  ENDIF


CONTAINS

  SUBROUTINE calculate_gkcut()

    USE io_global, ONLY : stdout

    INTEGER :: kpoint

   
    IF (nks == 0) THEN
       !
       ! if k-points are automatically generated (which happens later)
       ! use max(bg)/2 as an estimate of the largest k-point
       !
       gkcut = 0.5d0 * max ( &
          sqrt (sum(fc%bg_t (1:3, 1)**2) ), &
          sqrt (sum(fc%bg_t (1:3, 2)**2) ), &
          sqrt (sum(fc%bg_t (1:3, 3)**2) ) )
    ELSE
       gkcut = 0.0d0
       DO kpoint = 1, nks
          gkcut = max (gkcut, sqrt ( sum(xk (1:3, kpoint)**2) ) )
       ENDDO
    ENDIF
    gkcut = (sqrt (fc%ecutt) / fc%tpiba_t + gkcut)**2
    !
    ! ... find maximum value among all the processors
    !
    CALL mp_max (gkcut, inter_pool_comm )

  END SUBROUTINE calculate_gkcut


END SUBROUTINE data_structure_custom




SUBROUTINE initialize_fft_custom(fc)
!this subroutines initialize all the fft stuff for the custom defined grid

  USE kinds,              ONLY : DP
  USE gvect,              ONLY : g,mill
  USE cell_base,          ONLY : at, bg,tpiba2,tpiba,omega,alat
  USE io_global,          ONLY : stdout
  use control_flags, ONLY : gamma_only
  USE mp, ONLY : mp_barrier
  USE mp_world, ONLY : world_comm

  implicit none

  TYPE (fft_cus) :: fc

  INTEGER :: ng,n1t,n2t,n3t

 
  fc%at_t(1:3,1:3)=at(1:3,1:3)
  fc%bg_t(1:3,1:3)=bg(1:3,1:3)
  fc%alat_t=alat
  fc%omega_t=omega
  fc%tpiba_t=tpiba
  fc%tpiba2_t=tpiba2

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT1'
  call flush_unit(stdout)
  call  set_custom_grid(fc)

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT2'
  call flush_unit(stdout)

  call data_structure_custom(fc)

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT3'
  call flush_unit(stdout)
  
  allocate(fc%nlt(fc%ngmt))
  allocate(fc%nltm(fc%ngmt))
  write(stdout,*) 'ATT4'
  call flush_unit(stdout)
  
  call ggent(fc)
  
  return
END SUBROUTINE initialize_fft_custom


SUBROUTINE initialize_fft_custom_cell(fc)
!this subroutines initialize all the fft stuff for the custom defined grid
!for an arbitratry cell

  USE kinds,              ONLY : DP
  USE gvect,              ONLY : g,mill
  USE cell_base,          ONLY : at, bg,tpiba2,tpiba,omega,alat
  USE io_global,          ONLY : stdout
  use control_flags, ONLY : gamma_only
  USE mp, ONLY : mp_barrier
  USE mp_world, ONLY : world_comm

  implicit none

  TYPE (fft_cus) :: fc

  INTEGER :: ng,n1t,n2t,n3t

!the following must be provided from input
 
  !fc%at_t(1:3,1:3)=at(1:3,1:3)
  !fc%bg_t(1:3,1:3)=bg(1:3,1:3)
  !fc%alat_t=alat
  !fc%omega_t=omega
  !fc%tpiba_t=tpiba
  !fc%tpiba2_t=tpiba2

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT1'
  call flush_unit(stdout)
  call  set_custom_grid(fc)

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT2'
  call flush_unit(stdout)

  call data_structure_custom(fc)

  call mp_barrier( world_comm )
  write(stdout,*) 'ATT3'
  call flush_unit(stdout)
  
  allocate(fc%nlt(fc%ngmt))
  allocate(fc%nltm(fc%ngmt))
  write(stdout,*) 'ATT4'
  call flush_unit(stdout)
  
  call ggent(fc)
  
  return
END SUBROUTINE initialize_fft_custom_cell



SUBROUTINE ggent(fc)
  USE kinds,              ONLY : DP
  USE control_flags, ONLY : gamma_only
  USE constants,          ONLY : eps8

 implicit none

 TYPE(fft_cus) :: fc

   !
   REAL(DP) ::  t (3), tt, swap
   !
   INTEGER :: ngmx, n1, n2, n3, n1s, n2s, n3s
   !
   REAL(DP), ALLOCATABLE :: g2sort_g(:)
   ! array containing all g vectors, on all processors: replicated data
   INTEGER, ALLOCATABLE :: mill_g(:,:), mill_unsorted(:,:)
   ! array containing all g vectors generators, on all processors:
   !     replicated data
   INTEGER, ALLOCATABLE :: igsrt(:)
   !

#ifdef __PARA
   INTEGER :: m1, m2, mc
   !
#endif
   INTEGER :: i, j, k, ipol, ng, igl, iswap, indsw


   allocate(fc%gt(3,fc%ngmt),fc%ggt(fc%ngmt))


   ALLOCATE( fc%ig_l2gt( fc%ngmt_l ) )
   ALLOCATE( mill_g( 3, fc%ngmt_g ),mill_unsorted( 3, fc%ngmt_g ) )
   ALLOCATE( igsrt( fc%ngmt_g ) )
   ALLOCATE( g2sort_g( fc%ngmt_g ) )
   ALLOCATE(fc%ig1t(fc%ngmt),fc%ig2t(fc%ngmt),fc%ig3t(fc%ngmt))
   
   g2sort_g(:) = 1.0d20
   !
   ! save present value of ngm in ngmx variable
   !
   ngmx = fc%ngmt
   !
   fc%ngmt = 0
   iloop: DO i = -fc%nr1t-1, fc%nr1t+1
      !
      ! gamma-only: exclude space with x < 0
      !
      IF ( gamma_only .and. i < 0) CYCLE iloop
      jloop: DO j = -fc%nr2t-1, fc%nr2t+1
         !
         ! gamma-only: exclude plane with x = 0, y < 0
         !
         IF ( gamma_only .and. i == 0 .and. j < 0) CYCLE jloop
         kloop: DO k = -fc%nr3t-1, fc%nr3t+1
            !
            ! gamma-only: exclude line with x = 0, y = 0, z < 0
            !
            IF ( gamma_only .and. i == 0 .and. j == 0 .and. k < 0) CYCLE kloop
            t(:) = i * fc%bg_t (:,1) + j * fc%bg_t (:,2) + k * fc%bg_t (:,3)
            tt = sum(t(:)**2)
            IF (tt <= fc%gcutmt) THEN
               fc%ngmt = fc%ngmt + 1
               IF (fc%ngmt > fc%ngmt_g) CALL errore ('ggent', 'too many g-vectors', fc%ngmt)
               mill_unsorted( :, fc%ngmt ) = (/ i,j,k /)
               IF ( tt > eps8 ) THEN
                  g2sort_g(fc%ngmt) = tt
               ELSE
                  g2sort_g(fc%ngmt) = 0.d0
               ENDIF
            ENDIF
         ENDDO kloop
      ENDDO jloop
   ENDDO iloop

   IF (fc%ngmt  /= fc%ngmt_g ) &
         CALL errore ('ggen', 'g-vectors missing !', abs(fc%ngmt - fc%ngmt_g))
   igsrt(1) = 0
   CALL hpsort_eps( fc%ngmt_g, g2sort_g, igsrt, eps8 )
   mill_g(1,:) = mill_unsorted(1,igsrt(:))
   mill_g(2,:) = mill_unsorted(2,igsrt(:))
   mill_g(3,:) = mill_unsorted(3,igsrt(:))
   DEALLOCATE( g2sort_g, igsrt, mill_unsorted )
   fc%ngmt = 0
  
   ngloop: DO ng = 1, fc%ngmt_g
      i = mill_g(1, ng)
      j = mill_g(2, ng)
      k = mill_g(3, ng)

#ifdef __PARA
      m1 = mod (i, fc%nr1t) + 1
      IF (m1 < 1) m1 = m1 + fc%nr1t
      m2 = mod (j, fc%nr2t) + 1
      IF (m2 < 1) m2 = m2 + fc%nr2t
      mc = m1 + (m2 - 1) * fc%nrx1t
      IF ( fc%dfftt%isind ( mc ) == 0) CYCLE ngloop
#endif

      fc%ngmt = fc%ngmt + 1

      !  Here map local and global g index !!!
      fc%ig_l2gt( fc%ngmt ) = ng

      fc%gt (1:3, fc%ngmt) = i * fc%bg_t (:, 1) + j * fc%bg_t (:, 2) + k * fc%bg_t (:, 3)
      fc%ggt (fc%ngmt) = sum(fc%gt (1:3, fc%ngmt)**2)

      IF (fc%ngmt > ngmx) CALL errore ('ggen', 'too many g-vectors', fc%ngmt)
   ENDDO ngloop
   IF (fc%ngmt /= ngmx) &
      CALL errore ('ggent', 'g-vectors missing !', abs(fc%ngmt - ngmx))

   !
   !  here to initialize berry_phase
   !  CALL berry_setup(ngm, ngm_g, nr1, nr2, nr3, mill_g)
   !
   !     determine first nonzero g vector
   !
   IF (fc%ggt(1).le.eps8) THEN
      fc%gstart_t=2
   ELSE
      fc%gstart_t=1
   ENDIF
   !
   !     Now set nl and nls with the correct fft correspondence
   !
   DO ng = 1, fc%ngmt
      n1 = nint (sum(fc%gt (:, ng) * fc%at_t (:, 1))) + 1
      fc%ig1t (ng) = n1 - 1
      IF (n1<1) n1 = n1 + fc%nr1t
      

      n2 = nint (sum(fc%gt (:, ng) * fc%at_t (:, 2))) + 1
      fc%ig2t (ng) = n2 - 1
      IF (n2<1) n2 = n2 + fc%nr2t
      

      n3 = nint (sum(fc%gt (:, ng) * fc%at_t (:, 3))) + 1
      fc%ig3t (ng) = n3 - 1
      IF (n3<1) n3 = n3 + fc%nr3t
      

      IF (n1>fc%nr1t .or. n2>fc%nr2t .or. n3>fc%nr3t) &
         CALL errore('ggent','Mesh too small?',ng)

#if defined (__PARA) && !defined (__USE_3D_FFT)
      fc%nlt (ng) = n3 + ( fc%dfftt%isind (n1 + (n2 - 1) * fc%nrx1t) - 1) * fc%nrx3t
#else
      fc%nlt (ng) = n1 + (n2 - 1) * fc%nrx1t + (n3 - 1) * fc%nrx1t * fc%nrx2t

#endif
   ENDDO
    !
   DEALLOCATE( mill_g )
   !
   ! calculate number of G shells: ngl

    IF ( gamma_only) THEN

      DO ng = 1, fc%ngmt
      n1 = -fc%ig1t (ng) + 1
      IF (n1 < 1) n1 = n1 + fc%nr1t
      

      n2 = -fc%ig2t (ng) + 1
      IF (n2 < 1) n2 = n2 + fc%nr2t
      

      n3 = -fc%ig3t (ng) + 1
      IF (n3 < 1) n3 = n3 + fc%nr3t
      

      IF (n1>fc%nr1t .or. n2>fc%nr2t .or. n3>fc%nr3t) THEN
         CALL errore('ggent meno','Mesh too small?',ng)
      ENDIF

#if defined (__PARA) && !defined (__USE_3D_FFT)
      fc%nltm(ng) = n3 + (fc%dfftt%isind (n1 + (n2 - 1) * fc%nrx1t) - 1) * fc%nrx3t
     
#else
      fc%nltm(ng) = n1 + (n2 - 1) * fc%nrx1t + (n3 - 1) * fc%nrx1t * fc%nrx2t
      
#endif
   ENDDO


   ENDIF
 
!set npwt,npwxt

   if(gamma_only) then
      fc%npwt=0
      fc%npwxt=0
       do ng = 1, fc%ngmt
          tt = (fc%gt (1, ng) ) **2 + (fc%gt (2, ng) ) **2 + (fc%gt (3, ng) ) **2
          if (tt <= fc%ecutt / fc%tpiba2_t) then
             !
           ! here if |k+G|^2 <= Ecut increase the number of G inside the sphere
           !
             fc%npwt = fc%npwt + 1
          endif
       enddo
       fc%npwxt=fc%npwt
   endif

  return

END SUBROUTINE ggent
        
SUBROUTINE deallocate_fft_custom(fc)
!this subroutine deallocates all the fft custom stuff
  USE fft_types, ONLY : fft_dlay_deallocate

  implicit none

  TYPE(fft_cus) :: fc

  deallocate(fc%nlt,fc%nltm)
  call fft_dlay_deallocate(fc%dfftt)
  deallocate(fc%ig_l2gt,fc%ggt,fc%gt)
  deallocate(fc%ig1t,fc%ig2t,fc%ig3t)

  return

END SUBROUTINE deallocate_fft_custom
  


SUBROUTINE cft3t( fc, f, n1, n2, n3, nx1, nx2, nx3, sign )
  !----------------------------------------------------------------------------
  !
  ! ... sign = +-1 : parallel 3d fft for rho and for the potential
  ! ... sign = +-2 : parallel 3d fft for wavefunctions
  !
  ! ... sign = +   : G-space to R-space, output = \sum_G f(G)exp(+iG*R)
  ! ...              fft along z using pencils        (cft_1z)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along y (using planes) and x (cft_2xy)
  ! ... sign = -   : R-space to G-space, output = \int_R f(R)exp(-iG*R)/Omega
  ! ...              fft along x and y(using planes)  (cft_2xy)
  ! ...              transpose across nodes           (fft_scatter)
  ! ...                 and reorder
  ! ...              fft along z using pencils        (cft_1z)
  !
  ! ...  The array "planes" signals whether a fft is needed along y :
  ! ...    planes(i)=0 : column f(i,*,*) empty , don't do fft along y
  ! ...    planes(i)=1 : column f(i,*,*) filled, fft along y needed
  ! ...  "empty" = no active components are present in f(i,*,*)
  ! ...            after (sign>0) or before (sign<0) the fft on z direction
  !
  ! ...  Note that if sign=+/-1 (fft on rho and pot.) all fft's are needed
  ! ...  and all planes(i) are set to 1
  !

  USE kinds,        ONLY : DP
  USE fft_parallel, ONLY : tg_cft3s
  USE fft_scalar,   ONLY : cfft3ds, cfft3d  !  common scalar fft driver
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  TYPE(fft_cus) :: fc
  INTEGER,     INTENT(IN)    :: n1, n2, n3, nx1, nx2, nx3, sign

#if defined (__PARA) && !defined(__USE_3D_FFT)
!
  COMPLEX(DP), INTENT(INOUT) :: f( fc%dfftt%nnr )
  !
  ! ... call the general purpose parallel driver
  !
  call start_clock('cft3t')

  CALL tg_cft3s( f, fc%dfftt, sign )
  !write(stdout,*) fc%dfftt%nr1,fc%dfftt%nr1x,fc%dfftt%nr2,fc%dfftt%nr2x,fc%dfftt%nr3,fc%dfftt%nr3x
  !
#else
  !
  ! ... serial case
  !
  COMPLEX(DP), INTENT(INOUT) :: f(nx1*nx2*nx3)
  !
  !
 
  !
  ! ... sign = +-1 : complete 3d fft (for rho and for the potential)
  !
  IF ( sign == 1 ) THEN
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, 1 )
     !
  ELSE IF ( sign == -1 ) THEN
     !
     CALL cfft3d( f, n1, n2, n3, nx1, nx2, nx3, -1 )
     !
     ! ... sign = +-2 : if available, call the "short" fft (for psi's)
     !
  ELSE IF ( sign == 2 ) THEN
     !
     CALL cfft3ds( f, n1, n2, n3, nx1, nx2, nx3, 1, fc%dfftt%isind, fc%dfftt%iplw )
     !
  ELSE IF ( sign == -2 ) THEN
     !
     CALL cfft3ds( f, n1, n2, n3, nx1, nx2, nx3, -1, fc%dfftt%isind, fc%dfftt%iplw )
     !
  ELSE
     !
     CALL errore( 'cft3t', 'what should i do?', 1 )
     !
  END IF
  !

  !
#endif
  !
   call stop_clock('cft3t')
  RETURN
  !
END SUBROUTINE cft3t




END MODULE fft_custom_gwl
