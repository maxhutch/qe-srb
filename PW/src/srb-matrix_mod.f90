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
    complex(DP) :: alpha_l, beta_l
    allocate(Z(CHUNK, CHUNK))

    write(*,*) "alpha,beta= ", abs(alpha-one), abs(beta-zero)

    beta_l = beta

    do j = 1, n, CHUNK
      blockj = min(n-j+1, CHUNK) 
      do i = 1, j, CHUNK
        blocki = min(n-i+1, CHUNK) 
        call zgemm('C','N', blocki, blockj, k, &
                   one    , A(1+(i-1)*lda), lda, &
                          B(1+(j-1)*ldb), ldb, &
                   zero , Z             , CHUNK)
        call mp_sum(Z, intra_pool_comm)
        do ip = 1,blocki
          do jp = 1,blockj
             call infog2l( i+ip-1, j+jp-1, &
                           C_desc%desc, C_desc%nprow, C_desc%npcol, C_desc%myrow, C_desc%mycol, &
                           i_l, j_l, prow, pcol )
              if (prow == C_desc%myrow .and. pcol == C_desc%mycol) then
                C(i_l, j_l) = Z(ip,jp) + beta_l * C(i_l, j_l)
              endif
          enddo
        enddo
      enddo
    enddo
    deallocate(Z)

  end subroutine block_inner

end module srb_matrix
