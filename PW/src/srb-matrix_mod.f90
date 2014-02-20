module srb_matrix
  use kinds, only : DP
  type dmat
    integer :: desc(9)
    integer :: n,m
    complex(DP), allocatable :: dat(:,:)
    integer :: myrow
    integer :: mycol
    integer :: nprow
    integer :: npcol
    integer :: scope
  end type dmat

  integer, parameter :: nscope = 3

  integer :: serial_scope = 1
  integer :: pot_scope = 2
  integer :: pool_scope = 3

  integer :: ctx_s(nscope)
  integer :: nprow(nscope)
  integer :: npcol(nscope)
  integer :: first_proc(nscope)
  integer :: ctx_r(nscope)
  integer :: ctx_c(nscope)
  integer :: comm_scope(nscope)

  integer, external :: numroc

  contains
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
    use mp_global, only : nproc_pot, nproc_pool, my_pot_id, my_pool_id, me_pot
    use mp_global, only : intra_pot_comm, intra_pool_comm
    use mp, only: mp_barrier
    implicit none
    include 'mpif.h'
    integer :: i, prow, pcol, foo, bar
    integer :: nproc(nscope)
    integer, allocatable :: map(:,:)

    nproc(1) = serial_scope
    nproc(2) = nproc_pot 
    nproc(3) = nproc_pool
    comm_scope(1) = MPI_COMM_SELF ! not quite right
    comm_scope(2) = intra_pot_comm ! not quite right
    comm_scope(3) = intra_pool_comm ! not quite right
    first_proc(3) = 0!my_pool_id*nproc_pool
    first_proc(2) = 0!first_proc(3) + my_pot_id * nproc_pot
    first_proc(1) = 0!first_proc(2) + me_pot

    ctx_s = comm_scope
    ctx_r = comm_scope
    ctx_c = comm_scope

    allocate(map(nproc(nscope), nproc(nscope)))
    ! create three contexts for each scope
    do i = 1, nscope
      map = -1
      ! stripes are easy
      do prow = 1, nproc(i)
        map(prow, 1) = first_proc(i) + prow - 1
      enddo
      call blacs_gridmap( ctx_r(i), map, nproc(nscope), nproc(i), 1 )
      call mp_barrier(intra_pool_comm)
      do pcol = 1, nproc(i)
        map(1, pcol) = first_proc(i) + pcol - 1
      enddo
      call blacs_gridmap( ctx_c(i), map, nproc(nscope), 1, nproc(i))
      call mp_barrier(intra_pool_comm)

      ! square takes some thinking
      nprow(i) = max(1, int(0.99999 + sqrt(dble(nproc(i)))))
      npcol(i) = nproc(i) / nprow(i)
      do while (nproc(i) - nprow(i) * npcol(i) > 0)
        nprow(i) = nprow(i) + 1
        npcol(i) = nproc(i)/nprow(i) 
      enddo
      do prow = 1, nprow(i)
        do pcol = 1, npcol(i) 
          map(prow, pcol) = first_proc(i) + prow + (pcol - 1)*nprow(i) - 1
        enddo
      enddo
      !call blacs_get( -1, 0, ctx_s(i) )
      call blacs_gridmap( ctx_s(i), map, nproc(nscope), nprow(i), npcol(i))
      call mp_barrier(intra_pool_comm)
      call mp_barrier(intra_pool_comm)
      call mp_barrier(intra_pool_comm)
      call mp_barrier(intra_pool_comm)
    enddo

  end subroutine srb_matrix_init

  subroutine setup_dmat(A, m, n, mb_in, nb_in, scope_in)
    implicit none

    type(dmat), intent(inout) :: A
    integer, intent(in) :: m, n
    integer, intent(in), optional :: mb_in, nb_in
    integer, intent(in), optional :: scope_in

    integer :: mb, nb, ctx, nrl, ncl, info
   
    A%m = m; A%n = n

    ! default scope is pool
    if (present(scope_in)) then
      A%scope = scope_in
    else
      A%scope = pool_scope
    endif
    if (nprow(A%scope) * npcol(A%scope) == 1) A%scope = serial_scope

    ! pick reasonable default blocks
    if (.not. (present(mb_in) .or. present(nb_in))) then
      mb = min(m/nprow(A%scope) + 1,n/npcol(A%scope)+1)
      nb = mb
    endif
    if (present(mb_in)) then
      mb = mb_in
    else if (present(nb_in)) then
      if (nprow(A%scope) == 1) then
        mb = m
      else
        mb = max(4,4*(m/(4*nprow(A%scope))))
      endif
    endif
    if (present(nb_in)) then
      nb = nb_in
    else if (present(mb_in)) then 
      if (npcol(A%scope) == 1) then
        nb = n
      else
        nb = max(4,4*(n/(4*npcol(A%scope))))
      endif
    endif

    if (nprow(A%scope) > 1 .and. mb == m) then
      ctx = ctx_c(A%scope) ! parallel over columns
    else if (npcol(A%scope) > 1 .and. nb == n) then
      ctx = ctx_r(A%scope) ! parallel over rows
    else
      ctx = ctx_s(A%scope) ! parallel over both
    endif 
    call blacs_gridinfo(ctx, A%nprow, A%npcol, A%myrow, A%mycol)
    nrl = numroc(m, mb, A%myrow, 0, A%nprow)
    ncl = numroc(n, nb, A%mycol, 0, A%npcol)
    if (allocated(A%dat)) deallocate(A%dat)
    allocate(A%dat(nrl, ncl)); A%dat = 0.d0
    call descinit(A%desc, m, n, mb, nb, 0, 0, ctx, max(1,nrl), info)
    if (info .ne. 0) write(*,*) "descinit failed in srb_matrix"

  end subroutine setup_dmat

  subroutine copy_dmat(A, B)
    implicit none
    type(dmat), intent(inout) :: A
    type(dmat), intent(in) :: B
    integer :: nrl, ncl
    A%desc  = B%desc
    A%m     = B%m
    A%n     = B%n
    A%myrow = B%myrow
    A%mycol = B%mycol
    A%nprow = B%nprow
    A%npcol = B%npcol
    A%scope = B%scope
    nrl = numroc(A%desc(3), A%desc(5), A%myrow, 0, A%nprow)
    ncl = numroc(A%desc(4), A%desc(6), A%mycol, 0, A%npcol)
    if (allocated(A%dat)) deallocate(A%dat)
    allocate(A%dat(nrl, ncl)); A%dat = 0.d0
  end subroutine copy_dmat

#define CHUNK 32
  subroutine block_inner(n, k, alpha, A, lda, B, ldb, beta, C, comm_in)
    use kinds, only : DP
    use mp, only : mp_sum 

    implicit none

    integer,      intent(in)  :: n,k
    integer,      intent(in)  :: lda, ldb
    complex(DP),  intent(in)  :: A(*)
    complex(DP),  intent(in)   :: B(*)
    type(dmat),   intent(inout) :: C
    complex(DP),  intent(in)  :: alpha, beta
    integer, intent(in), optional :: comm_in

    integer :: i,j,blocki,blockj,ip,jp,i_l,j_l,prow,pcol
    complex(DP), allocatable :: Z(:,:)
    complex(DP), parameter :: one  = cmplx(1.d0,kind=DP)
    complex(DP), parameter :: zero = cmplx(0.d0,kind=DP)
    integer :: comm
    if (present(comm_in)) then
      comm = comm_in
    else
      call blacs_get(C%desc(2), 10, comm)
    endif
    allocate(Z(CHUNK, CHUNK))

    do j = 1, n, CHUNK
      blockj = min(n-j+1, CHUNK) 
      do i = 1, j, CHUNK
        blocki = min(n-i+1, CHUNK) 
        call zgemm('C','N', blocki, blockj, k, &
                   alpha, A(1+(i-1)*lda), lda, &
                          B(1+(j-1)*ldb), ldb, &
                   zero , Z             , CHUNK)
        call mp_sum(Z, comm)
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

  end subroutine block_inner

#define CHUNK 32
  subroutine block_outer(n, k, alpha, A, lda, B, ldb, beta, C)
    use kinds, only : DP
    use mp, only : mp_sum 

    implicit none

    integer,      intent(in)  :: n,k
    integer,      intent(in)  :: lda, ldb
    complex(DP),  intent(in)  :: A(*)
    complex(DP),  intent(in)  :: B(*)
    type(dmat),  intent(inout) :: C
    complex(DP),  intent(in)  :: alpha, beta

    integer :: i,j,blocki,blockj,ip,jp,i_l,j_l,prow,pcol
    complex(DP), allocatable :: Z(:,:)
    complex(DP), parameter :: one  = cmplx(1.d0,kind=DP)
    complex(DP), parameter :: zero = cmplx(0.d0,kind=DP)
    integer :: comm
    call blacs_get(C%desc(2), 10, comm)

    allocate(Z(CHUNK, CHUNK))

    do j = 1, n, CHUNK
      blockj = min(n-j+1, CHUNK) 
      do i = 1, j, CHUNK
        blocki = min(n-i+1, CHUNK) 
        call zgemm('N','C', blocki, blockj, k, &
                   alpha, A(i), lda, &
                          B(j), ldb, &
                   zero , Z             , CHUNK)
        call mp_sum(Z, comm)
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

  end subroutine block_outer

  subroutine add_diag(A, alpha)
    use kinds, only : DP
    implicit none

    type(dmat), intent(inout) :: A
    complex(DP), intent(in) :: alpha

    integer :: i, i_l, j_l, prow, pcol
    do i = 1, A%desc(3)
      call infog2l( i, i, &
                    A%desc, A%nprow, A%npcol, A%myrow, A%mycol, &
                    i_l, j_l, prow, pcol )
      if (prow == A%myrow .and. pcol == A%mycol) then
        A%dat(i_l, j_l) = A%dat(i_l, j_l) + alpha
      endif
    enddo

  end subroutine add_diag

  subroutine diag(A, W, Z, B, num_in, factored_in)
    use kinds, only : DP
    use mp_global, only : me_image, intra_pool_comm
    use mp, only : mp_sum
    implicit none
    
    type(dmat) :: A
    type(dmat), target :: Z
    real(DP) :: W(:)
    type(dmat), optional :: B
    integer, optional :: num_in
    logical, optional :: factored_in
    integer :: num
    logical :: factored

    logical :: use_tmp
    complex(DP), pointer :: ztmp(:,:)
    integer :: tmp_desc(9)
    integer, allocatable :: ifail(:), iclustr(:), iwork(:)
    real(DP), allocatable :: rwork(:), gap(:)
    complex(DP), allocatable :: work(:)
    integer :: ierr, lwork, lrwork, liwork, nv_out, ne_out, ctx_tmp
    complex(DP), parameter :: one = cmplx(1.d0,kind=DP)
    real(DP) :: abstol, scal
    real(DP), external :: PDLAMCH

    if (present(num_in)) then
      num = num_in
    else
      num = A%desc(3)
    endif

    if (present(factored_in)) then
      factored = factored_in
    else
      factored = .false.
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

    abstol = 2.d0*PDLAMCH(A%desc(2), 'U')

    if (present(B) .and. .not. factored) then
      call cholesky(B)
    endif

    if (present(B)) then
      allocate(work(1))
      call PZHENGST(1, 'U', A%desc(4), &
                    A%dat, 1, 1, A%desc, &
                    B%dat, 1, 1, B%desc, &
                    scal, work, -1, ierr)
      lwork = work(1); deallocate(work); allocate(work(lwork)) 
      call PZHENGST(1, 'U', A%desc(4), &
                    A%dat, 1, 1, A%desc, &
                    B%dat, 1, 1, B%desc, &
                    scal, work, lwork, ierr)
      if (ierr /= 0) write(*,*) "zhengst error: ", ierr
      deallocate(work)
    endif

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
                 ifail, iclustr, gap, ierr)
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
                 ifail, iclustr, gap, ierr)
    if (ierr /= 0) write(*,*) "pzheevx error: ", ierr
    if (present(B)) then
      call PZTRSM('L', 'U', 'N', 'N', A%desc(3), ne_out, one, &
                  B%dat, 1, 1, B%desc, &
                  ztmp, 1, 1, tmp_desc)
      if (scal .ne. one) call dscal(A%desc(3), scal, W, 1)
    endif

    if (use_tmp) then
      if ( Z%scope == 1 .and. A%scope == 3 ) then
        Z%dat = 0.d0
        ctx_tmp = Z%desc(2)
        if (me_image .ne. 0) Z%desc(2) = -1
        call pzgemr2d(A%desc(3), num, &
                      ztmp,  1, 1, tmp_desc, &
                      Z%dat, 1, 1, Z%desc, &
                      tmp_desc(2))
        call mp_sum(Z%dat, intra_pool_comm)
        Z%desc(2) = ctx_tmp  
      else
        call pzgemr2d(A%desc(3), num, &
                      ztmp, 1, 1, tmp_desc, &
                      Z%dat, 1, 1, Z%desc, &
                      tmp_desc(2))
      endif
    endif

  end subroutine diag

  subroutine g2l(M, i, j, i_l, j_l, local)
    implicit none
    type(dmat), intent(in) :: M
    integer, intent(in) :: i,j
    integer, intent(out) :: i_l, j_l
    logical, intent(out) :: local
    integer :: prow, pcol
      call infog2l(i, j, &
                   M%desc, M%nprow, M%npcol, &
                   M%myrow, M%mycol, &
                   i_l, j_l, prow, pcol)
    local = (prow == M%myrow .and. pcol == M%mycol)
  end subroutine

  subroutine l2g(M, i, j, i_g, j_g)
    implicit none
    type(dmat), intent(in) :: M
    integer, intent(in) :: i,j
    integer, intent(out) :: i_g, j_g
    integer, external :: indxl2g
    i_g = indxl2g(i, M%desc(5), M%myrow, M%desc(7), M%nprow)
    j_g = indxl2g(j, M%desc(6), M%mycol, M%desc(8), M%npcol)
  end subroutine

  subroutine collect(G, L, m, n, i, j)
    use mp_global, only : me_image, intra_pool_comm, intra_pot_comm
    use mp, only : mp_sum
    implicit none
    type(dmat), intent(in)    :: G
    type(dmat), intent(inout) :: L
    integer, intent(in) :: m, n, i, j
    integer ::  tmp
    tmp = L%desc(2)
    !write(*,*) G%mycol, G%myrow, comm_scope(G%scope), intra_pot_comm
    if (G%mycol /= 0 .or. G%myrow /= 0) L%desc(2) = -1
    call pzgemr2d(m, n, &
                  G%dat,  i, j, G%desc, &
                  L%dat, 1, 1, L%desc, &
                  G%desc(2))
    L%desc(2) = tmp
    call mp_sum(L%dat, comm_scope(G%scope))
  end subroutine collect

  subroutine distribute(L,G, who)
    implicit none
    type(dmat), intent(inout)  :: L
    type(dmat), intent(inout) :: G
    integer, intent(in) :: who
    integer tmp
    tmp =  L%desc(2) 
    if (who .ne. 0) L%desc(2) = -1
    call pzgemr2d(L%desc(3), L%desc(4), &
                  L%dat, 1, 1, L%desc,  &
                  G%dat, 1, 1, G%desc, &
                  G%desc(2))
    L%desc(2) = tmp
  end subroutine distribute

  subroutine load_from_local(A, ia, ja, B)
    implicit none
    type(dmat), intent(inout) :: A
    complex(DP), intent(in) :: B(:,:)
    integer, intent(in) :: ia, ja
    integer :: i, j, i_l, j_l, prow, pcol
    do j = 1, size(B,2)
      do i = 1, size(B,1)
        call infog2l( i+ia-1, j+ja-1, &
                     A%desc, A%nprow, A%npcol, A%myrow, A%mycol, &
                     i_l, j_l, prow, pcol )
            if (prow == A%myrow .and. pcol == A%mycol) then
                A%dat(i_l, j_l) = B(i,j)
            endif
      enddo
    enddo
  end subroutine load_from_local

  subroutine col_scal(A, B,scal)
    type(dmat), intent(in) :: A
    type(dmat), intent(inout) :: B
    real(DP), intent(in) :: scal(:)
    integer :: i, i_g

    do i = 1, size(A%dat,2)
      i_g = indxl2g(i, A%desc(6), A%mycol, A%desc(8), A%npcol)
      B%dat(:,i) = scal(i_g) * A%dat(:,i)
    enddo
  end subroutine col_scal

  subroutine parallel_inner(A, B, C, Ci)
    implicit none
    type(dmat), intent(in) :: A,B
    type(dmat), intent(inout) :: C
    integer, intent(in) :: Ci
    complex(DP), parameter :: one = cmplx(1.d0, kind=DP), zero = cmplx(0.d0, kind=DP)
    call pZGEMM('C', 'N', A%desc(4), B%desc(4), A%desc(3), &
               one,  A%dat, 1, 1, A%desc, &
                     B%dat, 1, 1, B%desc, &
               zero, C%dat, Ci, 1, C%desc)
  end subroutine parallel_inner

  subroutine cholesky(A)
    implicit none
    type(dmat), intent(inout) :: A
    integer :: ierr
    call PZPOTRF('U', A%desc(4), A%dat, 1,1, A%desc, ierr)
    if (ierr /= 0) write(*,*) "zpotrf error: ", ierr
  end subroutine cholesky

end module srb_matrix
