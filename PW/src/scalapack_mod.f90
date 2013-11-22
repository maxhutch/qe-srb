! ----------------------------------------------------------------------
  module scalapack_mod
! ----------------------------------------------------------------------

! An interface to scalapack routines
! designed to make distributed matrix inversion easier

  use kinds, only : dp
  use io_global, only : stdout
  implicit none
  public

! parameters
  integer,parameter :: DLEN_ = 9
  integer,parameter :: DTYPE_ = 1, &
                       CTXT_  = 2, &
                       M_     = 3, &
                       N_     = 4, &
                       MB_    = 5, &
                       NB_    = 6, &
                       RSRC_  = 7, &
                       CSRC_  = 8, &
                       LLD_   = 9 
! useful variables
  integer :: nprocs, nprow, npcol
  integer :: ctx_sq, ctx_re
  integer :: iam, myrow, mycol

  integer :: nblock
  integer :: desc_sq(DLEN_), desc_re(DLEN_)
  
  integer,external :: numroc, ilcm, iceil

  contains


! ----------------------------------------------------------------------
  subroutine scalapack_init
! ----------------------------------------------------------------------

  call blacs_pinfo( iam, nprocs )

  call scalapack_defgrid( nprocs, nprow, npcol )

  call blacs_get( -1, 0, ctx_sq )
  call blacs_gridinit( ctx_sq, 'c', nprow, npcol )
  call blacs_gridinfo( ctx_sq, nprow, npcol, myrow, mycol )
!  write(*,*) "sq: ", nprow, npcol, myrow, mycol

!  call blacs_get( -1, 0, ctx_re )
!  call blacs_gridinit( ctx_re, 'c', nprocs, 1 )
!  call blacs_gridinfo( ctx_re, nprow, npcol, myrow, mycol )
!  write(*,*) "re: ", nprow, npcol, myrow, mycol

  end subroutine scalapack_init

! ----------------------------------------------------------------------
  subroutine scalapack_distrib( nr, nc, nr_l, nc_l )
! ----------------------------------------------------------------------
  integer,intent(in) :: nr, nc
  integer,intent(out) :: nr_l, nc_l
  integer :: info

  call scalapack_blocksize( nblock, 32, nr, nc, nprow, npcol )
  ! padding constants
  nr_l = numroc( nr, nblock, myrow, 0, nprow )
  nc_l = numroc( nc, nblock, mycol, 0, npcol )
!  write(*,*) nr_l, nc_l
!  write(200+iam,*) ' nr = ', nr
!  write(200+iam,*) ' nc = ', nc
!  write(200+iam,*) ' nblock = ', nblock
!  write(200+iam,*) ' nr_l = ', nr_l
!  write(200+iam,*) ' nc_l = ', nc_l
  call descinit( desc_sq, nr, nc, nblock, nblock, 0, 0, ctx_sq, max(1,nr_l), info )
!  write(200+iam,*) ' desc_sq = ', desc_sq
  end subroutine scalapack_distrib

! ----------------------------------------------------------------------
  subroutine scalapack_distrib_re( nr, nc, nr_i )
! ----------------------------------------------------------------------
  integer,intent(in) :: nr, nr_i, nc
  integer :: info
  integer :: nr_l, nc_l
!  nr_l = numroc( nr, nr_i, myrow, 0, nprocs )
!  nc_l = numroc( nc, nc, mycol, 0, 1 )
!  write(*,*) nr_l, nc_l
  ! padding constants
  call descinit( desc_re, nr, nc, nr_i, nc, 0, 0, ctx_re, nr_l, info )
  end subroutine scalapack_distrib_re

! ----------------------------------------------------------------------
  subroutine scalapack_localindex( i, j, li, lj, islocal )
! ----------------------------------------------------------------------
  integer,intent(in) :: i, j
  integer,intent(out) :: li, lj 
  logical,intent(out) :: islocal
  integer :: prow, pcol

  call infog2l( i, j, desc_sq, nprow, npcol, myrow, mycol, li, lj, prow, pcol )
  islocal = .false.
  if( myrow==prow .AND. mycol==pcol ) islocal = .true.
  
  end subroutine scalapack_localindex


! ----------------------------------------------------------------------
  subroutine scalapack_invert( n, a_l )
! ----------------------------------------------------------------------

  integer :: n
  complex(dp) :: a_l(*)
  integer,allocatable :: ipiv(:), iwork(:)
  complex(dp),allocatable :: work(:)
  integer :: lcm, nr_l, nc_l
  integer,save :: lwork, liwork
  integer :: info
  integer,save :: first_time=0

  ! workspace dimensions
  nr_l = numroc( n, nblock, myrow, 0, nprow )
  nc_l = numroc( n, nblock, mycol, 0, npcol )
!  lwork = nr_l * nblock
!  if( nprow==npcol ) then
!    liwork = nc_l + nblock
!  else
!    lcm = ilcm( nprow, npcol )
!    liwork = nc_l + max( iceil( iceil(nr_l,nblock) , (lcm/nprow) ), nblock )
!  endif

  call blacs_barrier( ctx_sq, 'All' )

  allocate( ipiv(nr_l+nblock) )
  call pzgetrf( n, n, a_l, 1, 1, desc_sq, ipiv, info )

  ! workspace query
  if( first_time /= n )  then
    lwork=-1
    liwork = -1
    allocate( work(1), iwork(1) )
    call pzgetri( n, a_l, 1, 1, desc_sq, ipiv, &
                  work, lwork, iwork, liwork, info )
    lwork = work(1)
    liwork = iwork(1)
    deallocate( work, iwork )
    first_time = n
  endif

  allocate( work(lwork), iwork(liwork), stat=info )
  call errore('scalapack_invert','unable to allocate workspace',abs(info))

  call pzgetri( n, a_l, 1, 1, desc_sq, ipiv, &
                work, lwork, iwork, liwork, info )

  deallocate( ipiv, work, iwork )

  end subroutine scalapack_invert


! ----------------------------------------------------------------------
  subroutine scalapack_diag( n, a_l, eig, z_l )
! ----------------------------------------------------------------------

  integer :: n
  complex(dp) :: a_l(*), z_l(*)
  real(dp) :: eig(n)
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer,allocatable :: iwork(:)
  integer,save :: lwork, lrwork, liwork
  integer :: info
  integer,save :: first_time=0
  integer :: descz(DLEN_)

  integer :: nb, np0, nq0, np, nq

  integer :: ia, ja
  real(dp) :: trace
  complex(dp),external :: pzlatra
  complex(dp) :: tmp

  call blacs_barrier( ctx_sq, 'All' )

  ! assume same descriptor for z as a
  descz = desc_sq

  ! workspace query
  if( first_time /= n )  then
!    write(*,*) 'first_time'
    lwork=-1
    lrwork = -1
    liwork = -1
    allocate( work(1), rwork(1), iwork(1) )
    call pzheevd( 'V', 'U', n, a_l, 1, 1, desc_sq, eig, z_l, 1, 1, descz, &
                   work, lwork, rwork, lrwork, iwork, liwork, info )
    lwork = work(1)
    lrwork = rwork(1)
    liwork = iwork(1)
    deallocate( work, rwork, iwork )
!    first_time = n
!    write(*,*) 'first_time', lwork, lrwork, liwork
!    write(*,*) ' n = ', n
!    nb = desc_sq( MB_ )
!    !np0 = numroc( n, nb, 0, 0, nprow )
!    np0 = numroc( max(n,nb,2), nb, 0, 0, nprow )
!    nq0 = numroc( max(n,nb,2), nb, 0, 0, npcol )
!    !lwork = (np0+nq0+nb)*nb+3*n+n**2
!    !lrwork = 2*n + 2*n-2
!    lwork = n + (np0+nq0+nb)*nb
!    np = numroc( n, nb, myrow, iarow, nprow )
!    nq = numroc( n, nb, mycol, iacol, npcol )
!    lrwork = 1+9*n+3*np*nq
!    liwork = 7*n+8*npcol+2
!    write(*,*) ' lwork >= ', lwork
!    write(*,*) ' lrwork >= ', lrwork
  endif

  allocate( work(lwork), rwork(lrwork), iwork(liwork), stat=info )
  call errore('scalapack_diag','unable to allocate workspace',abs(info))

  !write(*,*) ' pzheev '
  !call pzheev( 'V', 'U', n, a_l, 1, 1, desc_sq, eig, z_l, 1, 1, descz, &
  !              work, lwork, rwork, lrwork, info )
  !write(*,*) ' pzheevd '
  call pzheevd( 'V', 'U', n, a_l, 1, 1, desc_sq, eig, z_l, 1, 1, descz, &
                work, lwork, rwork, lrwork, iwork, liwork, info )
  if (info /= 0) write(*,*) "pzheevd error: ", info
  !write(*,*) ' done '
  deallocate( work, rwork, iwork )

  end subroutine scalapack_diag

! ----------------------------------------------------------------------
  subroutine scalapack_svd(m, n, a_l, sv, v_l )
! ----------------------------------------------------------------------

  integer :: m, n
  complex(dp) :: a_l(*), v_l(*)
  real(dp) :: sv(*)
  complex(dp),allocatable :: work(:)
  real(dp),allocatable :: rwork(:)
  integer,save :: lwork, lrwork
  integer :: info
  integer,save :: first_time=0
  integer :: descv(DLEN_)

  call blacs_barrier( ctx_sq, 'All' )

  ! assume same descriptor for z as a
  if (n >= m) then
    descv = desc_sq
  else
    ! this could be better
    descv = desc_sq
  endif


  ! workspace query
  if( first_time /= n )  then
    lwork=-1
    allocate( work(1), rwork(1))
    call pzgesvd('V','N',m,n,a_l,1,1,desc_sq,sv,v_l,1,1,descv, &
                 v_l,1,1,descv,work,lwork,rwork,info)
    if (info /= 0) write(*,*) "pzgesvd error: ", info
    lwork = work(1)
    lrwork = rwork(1)
    deallocate(work, rwork)
  endif
  allocate( work(lwork), rwork(lrwork), stat=info )
  call errore('scalapack_svd','unable to allocate workspace',abs(info))
  call pzgesvd('V','N',m,n,a_l,1,1,desc_sq,sv,v_l,1,1,descv, &
                 v_l,1,1,descv,work,lwork,rwork,info)
  if (info /= 0) write(*,*) "pzgesvd error: ", info
  deallocate( work, rwork )

  end subroutine scalapack_svd

! ----------------------------------------------------------------------
  subroutine scalapack_defgrid( nprocs, nprow, npcol )
! ----------------------------------------------------------------------

  integer,intent(in) :: nprocs
  integer,intent(out) :: nprow, npcol
  integer :: sqside

  ! try to make a square grid
  sqside = sqrt( dble(nprocs) )
  nprow = max( sqside, 1 )
  npcol = nprocs / nprow
  do while( nprocs - nprow*npcol > 0 ) 
    nprow=nprow+1
    npcol = nprocs / nprow
  enddo
  
  end subroutine scalapack_defgrid


! ----------------------------------------------------------------------
  subroutine scalapack_blocksize( nb, nbpref, nr, nc, nprow, npcol )
! ----------------------------------------------------------------------
  integer,intent(out) :: nb
  integer,intent(in) :: nbpref, nr, nc, nprow, npcol

  nb = min (nr / nprow, nc / npcol)
  if (nbpref.gt.0) then
     nb = min (nb, nbpref)
  endif
  nb = min (nb, nr, nc)
  if (min(nr,nc).lt.10) nb = 1

  end subroutine scalapack_blocksize


  end module scalapack_mod
