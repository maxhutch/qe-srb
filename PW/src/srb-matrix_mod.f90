module srb_matrix
  use kinds, only : DP
  type dmat
    integer :: desc(9)
    complex(DP), allocatable :: dat(:,:)
    integer :: myrow
    integer :: mycol
    integer :: nprow
    integer :: npcol
  end type dmat

  type mydesc
    integer :: desc(9)
    integer :: myrow
    integer :: mycol
    integer :: nprow
    integer :: npcol
    integer :: nrl
    integer :: ncl
  end type mydesc

  integer, parameter :: nscope = 3

  integer :: serial_scope = 1
  integer :: pot_scope = 2
  integer :: pool_scope = 3

  integer :: ctx_s(nscope)
  integer :: nprow(nscope)
  integer :: npcol(nscope)
  integer :: ctx_r(nscope)
  integer :: ctx_c(nscope)

  contains
#if 0
  subroutine grab_desc(desc)
    use scalapack_mod, only : myrow, mycol, nprow, npcol, desc_sq
    implicit none
    type(mydesc), intent(inout) :: desc
    desc%desc = desc_sq
    desc%myrow = myrow
    desc%mycol = mycol
    desc%nprow = nprow
    desc%npcol = npcol
  end subroutine grab_desc
#endif
  subroutine print_desc(desc)
    implicit none
    type(mydesc), intent(in) :: desc
    write(*,'(A,4(I2,A),6(I6,A),I3)') "desc: (",desc%myrow,",",desc%mycol,") in (", &
                                           desc%nprow,",",desc%npcol,") has [", &
                                           desc%nrl,",",desc%ncl,"] of [", &
                                           desc%desc(3),",",desc%desc(4),"] in [",&
                                           desc%desc(5),",",desc%desc(6),"] on ctx ",&
                                           desc%desc(2)
  end subroutine print_desc

  subroutine print_dmat(A)
    implicit none
    type(dmat), intent(in) :: A
    write(*,'(A,4(I2,A),6(I6,A),I3)') "desc: (",A%myrow,",",A%mycol,") in (", &
                                           A%nprow,",",A%npcol,") has [", &
                                           size(A%dat,1),",",size(A%dat,2),"] of [", &
                                           A%desc(3),",",A%desc(4),"] in [",&
                                           A%desc(5),",",A%desc(6),"] on ctx ",&
                                           A%desc(2)
  end subroutine print_dmat


  subroutine srb_matrix_init()
    use mp_global, only : nproc_pot, nproc_pool
    implicit none

    integer :: i, foo, bar
    integer :: nproc(nscope)

    nproc(1) = serial_scope
    nproc(2) = nproc_pot 
    nproc(3) = nproc_pool

    ! create three contexts for each scope
    do i = 1, nscope
      ! stripes are easy
      call blacs_get( -1, 0, ctx_r(i) )
      call blacs_gridinit( ctx_r(i), 'c', nproc(i), 1 )

      call blacs_get( -1, 0, ctx_c(i) )
      call blacs_gridinit( ctx_c(i), 'c', 1, nproc(i))

      ! square takes some thinking
      nprow(i) = max(1, int(sqrt(dble(nproc(i)))))
      npcol(i) = nproc(i) / nprow(i)
      do while (nproc(i) - nprow(i) * npcol(i) > 0)
        nprow(i) = nprow(i) + 1
        npcol(i) = nproc(i)/nprow(i) 
      enddo
      call blacs_get( -1, 0, ctx_s(i) )
      call blacs_gridinit( ctx_s(i), 'c', nprow(i), npcol(i))
    enddo

  end subroutine srb_matrix_init

  subroutine setup_dmat(A, m, n, mb_in, nb_in, scope_in)
    implicit none

    type(dmat), intent(inout) :: A
    integer, intent(in) :: m, n
    integer, intent(in), optional :: mb_in, nb_in
    integer, intent(in), optional :: scope_in

    integer :: mb, nb, scope, ctx, nrl, ncl, info
    integer, external :: numroc

    ! default scope is pool
    if (present(scope_in)) then
      scope = scope_in
    else
      scope = pool_scope
    endif

    ! pick reasonable default blocks
    if (present(mb_in)) then
      mb = mb_in
    else
      if (nprow(scope) == 1) then
        mb = m
      else
        write(*,*) m, nprow(scope)
        mb = max(4,4*(m/(4*nprow(scope))))
      endif
    endif
    if (present(nb_in)) then
      nb = nb_in
    else
      if (npcol(scope) == 1) then
        nb = n
      else
        nb = max(4,4*(n/(4*npcol(scope))))
      endif
    endif

    if (nprow(scope) > 1 .and. nb == n) then
      ctx = ctx_c(scope) ! parallel over columns
    else if (npcol(scope) > 1 .and. mb == m) then
      ctx = ctx_r(scope) ! parallel over rows
    else
      ctx = ctx_s(scope) ! parallel over both
    endif 
    call blacs_gridinfo(ctx, A%nprow, A%npcol, A%myrow, A%mycol)
    nrl = numroc(m, mb, A%myrow, 0, A%nprow)
    ncl = numroc(n, nb, A%mycol, 0, A%npcol)
    if (allocated(A%dat)) deallocate(A%dat)
    allocate(A%dat(nrl, ncl)); A%dat = 0.d0
    call descinit(A%desc, m, n, mb, nb, 0, 0, ctx, max(1,nrl), info)
    if (info .ne. 0) write(*,*) "descinit failed in srb_matrix"

  end subroutine setup_dmat

  subroutine setup_desc(desc, nr, nc, blockr_in, blockc_in)
    use scalapack_mod, only : myrow_sq, mycol_sq, nprow, npcol
    use scalapack_mod, only :  ctx_sq, ctx_rex, ctx_rey
    use scalapack_mod, only : myrow_rex, mycol_rex, myrow_rey, mycol_rey
    use scalapack_mod, only : scalapack_blocksize
    use mp_global, only : nproc_pot, me_pot
    implicit none
    type(mydesc), intent(inout) :: desc
    integer, intent(in) :: nr, nc
    integer, optional, intent(in) ::  blockr_in, blockc_in

    integer, external :: numroc
    integer :: blockr, blockc
    integer :: ctx, info

    if (.not. present(blockr_in) .or. .not. present(blockc_in)) then
      call scalapack_blocksize( blockr, 32, nr, nc, nprow, npcol )
      blockc = blockr
    else
      blockc = blockc_in; blockr = blockr_in
    endif
    
    if (blockr == nr) then
      desc%nprow = 1
      desc%npcol = nproc_pot
      desc%myrow = 0
      desc%mycol = mycol_rey
      ctx = ctx_rey
    else if (blockc == nc) then
      desc%nprow = nproc_pot
      desc%npcol = 1
      desc%myrow = mycol_rex
      desc%mycol = 0
      ctx = ctx_rex
    else
      desc%nprow = nprow
      desc%npcol = npcol
      desc%myrow = myrow_sq
      desc%mycol = mycol_sq
      ctx = ctx_sq
    endif
    
    desc%nrl = numroc( nr, blockr, desc%myrow, 0, desc%nprow )
    desc%ncl = numroc( nc, blockc, desc%mycol, 0, desc%npcol )
    !write(*,*) nr, nc, blockr, blockc, desc%nrl, desc%ncl
    call descinit( desc%desc, nr, nc, blockr, blockc, 0, 0, ctx, max(1,desc%nrl), info )

  end subroutine setup_desc

#define CHUNK 32
  subroutine block_inner(n, k, alpha, A, lda, B, ldb, beta, C, C_desc)
!  subroutine block_inner(n, k, A, lda, B, ldb, C, C_desc)
    use kinds, only : DP
    use mp, only : mp_sum 
    use mp_global, only: intra_pool_comm

    implicit none

    integer,      intent(in)  :: n,k
    integer,      intent(in)  :: lda, ldb
    complex(DP),  intent(in)  :: A(*)
    complex(DP),  intent(in)  :: B(*)
    complex(DP),  intent(out) :: C(:,:)
    complex(DP),  intent(in)  :: alpha, beta
    type(mydesc), intent(in) :: C_desc

    integer :: i,j,blocki,blockj,ip,jp,i_l,j_l,prow,pcol
    complex(DP), allocatable :: Z(:,:)
    complex(DP), parameter :: one  = cmplx(1.d0,kind=DP)
    complex(DP), parameter :: zero = cmplx(0.d0,kind=DP)
    allocate(Z(CHUNK, CHUNK))

    do j = 1, n, CHUNK
      blockj = min(n-j+1, CHUNK) 
      do i = 1, j, CHUNK
        blocki = min(n-i+1, CHUNK) 
        call zgemm('C','N', blocki, blockj, k, &
                   alpha, A(1+(i-1)*lda), lda, &
                          B(1+(j-1)*ldb), ldb, &
                   zero , Z             , CHUNK)
        call mp_sum(Z, intra_pool_comm)
        do jp = 1,blockj
          do ip = 1,blocki
             call infog2l( i+ip-1, j+jp-1, &
                           C_desc%desc, C_desc%nprow, C_desc%npcol, C_desc%myrow, C_desc%mycol, &
                           i_l, j_l, prow, pcol )
              if (prow == C_desc%myrow .and. pcol == C_desc%mycol) then
                C(i_l, j_l) = Z(ip,jp) + beta * C(i_l, j_l)
              endif
          enddo
        enddo
      enddo
    enddo
    deallocate(Z)

  end subroutine block_inner

  subroutine block_inner2(n, k, alpha, A, lda, B, ldb, beta, C)
!  subroutine block_inner(n, k, A, lda, B, ldb, C, C_desc)
    use kinds, only : DP
    use mp, only : mp_sum 
    use mp_global, only: intra_pool_comm

    implicit none

    integer,      intent(in)  :: n,k
    integer,      intent(in)  :: lda, ldb
    complex(DP),  intent(in)  :: A(*)
    complex(DP),  intent(in)   :: B(*)
    type(dmat),   intent(inout) :: C
    complex(DP),  intent(in)  :: alpha, beta

    integer :: i,j,blocki,blockj,ip,jp,i_l,j_l,prow,pcol
    complex(DP), allocatable :: Z(:,:)
    complex(DP), parameter :: one  = cmplx(1.d0,kind=DP)
    complex(DP), parameter :: zero = cmplx(0.d0,kind=DP)
    allocate(Z(CHUNK, CHUNK))

    do j = 1, n, CHUNK
      blockj = min(n-j+1, CHUNK) 
      do i = 1, j, CHUNK
        blocki = min(n-i+1, CHUNK) 
        call zgemm('C','N', blocki, blockj, k, &
                   alpha, A(1+(i-1)*lda), lda, &
                          B(1+(j-1)*ldb), ldb, &
                   zero , Z             , CHUNK)
        call mp_sum(Z, intra_pool_comm)
        do jp = 1,blockj
          do ip = 1,blocki
             call infog2l( i+ip-1, j+jp-1, &
                           C%desc, C%nprow, C%npcol, C%myrow, C%mycol, &
                           i_l, j_l, prow, pcol )
              if (prow == C%myrow .and. pcol == C%mycol) then
                C%dat(i_l, j_l) = Z(ip,jp) + beta * C%dat(i_l, j_l)
              endif
          enddo
        enddo
      enddo
    enddo
    deallocate(Z)

  end subroutine block_inner2

#define CHUNK 32
  subroutine block_outer(n, k, alpha, A, lda, B, ldb, beta, C, C_desc)
!  subroutine block_inner(n, k, A, lda, B, ldb, C, C_desc)
    use kinds, only : DP
    use mp, only : mp_sum 
    use mp_global, only: intra_pool_comm

    implicit none

    integer,      intent(in)  :: n,k
    integer,      intent(in)  :: lda, ldb
    complex(DP),  intent(in)  :: A(*)
    complex(DP),  intent(in)  :: B(*)
    complex(DP),  intent(out) :: C(:,:)
    complex(DP),  intent(in)  :: alpha, beta
    type(mydesc), intent(in) :: C_desc

    integer :: i,j,blocki,blockj,ip,jp,i_l,j_l,prow,pcol
    complex(DP), allocatable :: Z(:,:)
    complex(DP), parameter :: one  = cmplx(1.d0,kind=DP)
    complex(DP), parameter :: zero = cmplx(0.d0,kind=DP)
    allocate(Z(CHUNK, CHUNK))

    do j = 1, n, CHUNK
      blockj = min(n-j+1, CHUNK) 
      do i = 1, j, CHUNK
        blocki = min(n-i+1, CHUNK) 
        call zgemm('N','C', blocki, blockj, k, &
                   alpha, A(i), lda, &
                          B(j), ldb, &
                   zero , Z             , CHUNK)
        call mp_sum(Z, intra_pool_comm)
        do jp = 1,blockj
          do ip = 1,blocki
            call infog2l( i+ip-1, j+jp-1, &
                          C_desc%desc, C_desc%nprow, C_desc%npcol, C_desc%myrow, C_desc%mycol, &
                          i_l, j_l, prow, pcol )
            if (prow == C_desc%myrow .and. pcol == C_desc%mycol) then
                C(i_l, j_l) = Z(ip,jp) + beta * C(i_l, j_l)
            endif
          enddo
        enddo
      enddo
    enddo
    deallocate(Z)

  end subroutine block_outer

  subroutine add_diag(n, alpha, A, desc)
    use kinds, only : DP
    implicit none

    integer, intent(in) :: n
    complex(DP), intent(in) :: alpha
    complex(DP), intent(inout) :: A(:,:)
    type(mydesc), intent(in) :: desc

    integer :: i, i_l, j_l, prow, pcol
    do i = 1, n
      call infog2l( i, i, &
                    desc%desc, desc%nprow, desc%npcol, desc%myrow, desc%mycol, &
                    i_l, j_l, prow, pcol )
      if (prow == desc%myrow .and. pcol == desc%mycol) then
        A(i_l, j_l) = A(i_l, j_l) + alpha
      endif
    enddo

  end subroutine add_diag

  subroutine diag(A, W, Z, B, num_in)
    use kinds, only : DP
    implicit none
    
    type(dmat) :: A
    type(dmat), target :: Z
    real(DP) :: W(:)
    type(dmat), optional :: B
    integer, optional :: num_in
    integer :: num

    logical :: use_tmp
    complex(DP), pointer :: ztmp(:,:)
    integer :: tmp_desc(9)
    integer, allocatable :: ifail(:), iclustr(:), iwork(:)
    real(DP), allocatable :: rwork(:), gap(:)
    complex(DP), allocatable :: work(:)
    integer :: info, lwork, lrwork, liwork, nv_out, ne_out
    real(DP) :: abstol
    real(DP), external :: PDLAMCH
    if (present(num_in)) then
      num = num_in
    else
      num = A%desc(3)
    endif

    ! check if we need tmp space for evecs
    use_tmp = Z%desc(4) < A%desc(4) .or. Z%desc(2) .ne. A%desc(2)
    if (use_tmp) then
      allocate(ztmp(size(A%dat, 1), size(A%dat,2)))
      tmp_desc = A%desc
    else
      ztmp => Z%dat
      tmp_desc = Z%desc
    endif 
    abstol = PDLAMCH(A%desc(2), 'U')

    if (.not. present(B)) then
      allocate(ifail(A%desc(3)))
      allocate(iclustr(2*A%nprow*A%npcol), gap(A%nprow*A%npcol))
      allocate(work(1), rwork(1), iwork(1))
      CALL pzheevx('V', 'I', 'U', A%desc(3), &
                  A%dat, 1, 1, A%desc, &
                  0, 0, 1, num, &
                  abstol, ne_out, nv_out, W, &
                  -1.d0, &
                  ztmp, 1, 1, tmp_desc, &
                  work, -1, rwork, -1, iwork, -1, &
                  ifail, iclustr, gap, info)
      lwork = work(1); deallocate(work); allocate(work(lwork))
      lrwork = rwork(1); deallocate(rwork); allocate(rwork(lrwork))
      liwork = iwork(1); deallocate(iwork); allocate(iwork(liwork))
      CALL pzheevx('V', 'I', 'U', A%desc(3), &
                   A%dat, 1, 1, A%desc, &
                   0, 0, 1, num, &
                   abstol, ne_out, nv_out, W, &
                   -1.d0, &
                   ztmp, 1, 1, tmp_desc, &
                   work, lwork, rwork, lrwork, iwork, liwork, &
                   ifail, iclustr, gap, info)
      if (use_tmp) then
        call pzgemr2d(A%desc(3), num, &
                      ztmp,  1, 1, tmp_desc, &
                      Z%dat, 1, 1, Z%desc, &
                      tmp_desc(2))
      endif
    endif

  end subroutine diag


end module srb_matrix
