module srb_matrix

  type mydesc
    integer :: desc(9)
    integer :: myrow
    integer :: mycol
    integer :: nprow
    integer :: npcol
    integer :: nrl
    integer :: ncl
  end type mydesc

  contains

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

#define BLOCK 32
  subroutine block_inner(n, k, A, lda, B, ldb, C, C_desc)
    use kinds, only : DP
    use mp, only : mp_sum 
    use mp_global, only: intra_pool_comm

    implicit none

    integer,      intent(in)  :: n,k
    integer,      intent(in)  :: lda, ldb
    complex(DP),  intent(in)  :: A(*)
    complex(DP),  intent(in)  :: B(*)
    complex(DP),  intent(out) :: C(:,:)
    type(mydesc), intent(in), optional :: C_desc

    integer :: i,j,blocki,blockj,ip,jp,i_l,j_l,prow,pcol
    complex(DP), allocatable :: Z(:,:)
    complex(DP), parameter :: one  = cmplx(1.d0,kind=DP)
    complex(DP), parameter :: zero = cmplx(0.d0,kind=DP)

    allocate(Z(BLOCK, BLOCK))

    do j = 1, n, BLOCK
      blockj = min(n-j+1, BLOCK) 
      do i = 1, j, BLOCK
        blocki = min(n-i+1, BLOCK) 
        call zgemm('C','N', blocki, blockj, k, &
                   one,  A(1+(i-1)*lda), lda, &
                         B(1+(j-1)*ldb), ldb, &
                   zero, Z         , BLOCK)
        call mp_sum(Z(1:blocki,1:blockj), intra_pool_comm)
        do ip = 1,blocki
          do jp = 1,blockj
             call infog2l( i+ip-1, j+jp-1, &
                           C_desc%desc, C_desc%nprow, C_desc%npcol, C_desc%myrow, C_desc%mycol, &
                           i_l, j_l, prow, pcol )
              if (prow == C_desc%myrow .and. pcol == C_desc%mycol) then
                C(i_l, j_l) = Z(ip,jp)
              endif
          enddo
        enddo
      enddo
    enddo
    deallocate(Z)

  end subroutine block_inner

end module srb_matrix
