!
!> Build the coefficients of the Hamiltonian:
!! H(k) = H^(0) + H^(1) k + H^(2) k^2 + ...
!! in the basis |Bi>
! #define AVOID_FFT
!
SUBROUTINE build_h_coeff(opt_basis, V_rs, ecut_srb, nspin, ham, saved_in)
  USE kinds, ONLY : DP
  USE mp_global, ONLY: intra_pool_comm, me_pool
  USE mp, ONLY: mp_sum, mp_barrier

  USE srb_types, ONLY : basis, ham_expansion
  USE srb_matrix, only : block_inner
  USE srb, ONLY : decomp_size

  USE cell_base, ONLY : tpiba2, tpiba
  USE klist, ONLY : nkstot, xk, ngk
  USE gvect, ONLY : ngm, g
  USE gvecs, only : nls
  use wavefunctions_module, only: psic
  use fft_interfaces, only : fwfft, invfft
  use fft_base, only : dffts 
  use buffers, only : get_buffer

  IMPLICIT NONE

  ! arguments
  TYPE(basis),         INTENT(IN)  :: opt_basis        !>!< optimal basis
  REAL(DP),            INTENT(IN)  :: V_rs(:,:)         !>!< total potential in real space
  REAL(DP),            INTENT(IN)  :: ecut_srb !>!< energy cutoff to include Gamma^8
  integer, intent(IN)              :: nspin
  TYPE(ham_expansion), INTENT(inout) :: ham              !>!< Hamiltonian
  logical, intent(in), optional :: saved_in
  ! functions
  COMPLEX(DP) , EXTERNAL :: ZDOTU

  ! locals
  integer :: nband_l ! number of local bands
  complex(DP), allocatable :: ham_l(:), gtmp(:,:), buffer(:,:)
  complex(DP), allocatable :: decomp1(:,:), decomp2(:,:)
  complex(DP), allocatable :: zroot(:,:)
  real(DP), allocatable :: root(:,:) 
  integer :: npw, nbnd, nbl
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  real(DP), parameter :: k_gamma(3) = 0.d0
  logical ::  saved
  COMPLEX(DP), parameter :: zero = cmplx(0.d0, KIND=DP), one = cmplx(1.d0, kind=DP)
  COMPLEX(DP), allocatable :: ctmp(:,:)
  ! indexes 
  integer :: ip ! processor index
  integer :: i, j, ij, ixyz, k, s
  
  if (present(saved_in)) then
    saved = saved_in
  else
    saved = .false.
  endif

  ! ============================================================
  ! Allocations and setup 
  ! ============================================================
  nbnd = opt_basis%length
  ! allocate the hamiltonian (we didn't know the size until now)
  ham%length = opt_basis%length
  if (.not. saved) then
    ALLOCATE(ham%con(ham%desc%nrl, ham%desc%ncl, nspin), ham%lin(3, ham%desc%nrl, ham%desc%ncl) )
    allocate(ham%kin_con(ham%desc%nrl, ham%desc%ncl))
  endif
  allocate(igk(ngm), g2kin(ngm))
  call gk_sort(k_gamma, ngm, g, ecut_srb/tpiba2, npw, igk, g2kin) ! re-order wrt |G+k| (k = 0 in this case)
  g2kin(1:npw) = g2kin(1:npw) * tpiba2 ! change units


  ! ===========================================================================
  ! kinetic energy term
  ! ===========================================================================
  call start_clock( '  Kinetic' )

  ! constant term
  if (saved) then
    forall (s = 1:nspin) ham%con(:,:,s) = ham%kin_con(:,:)
  else
    call start_clock( '   part1' )
    allocate(buffer(npw, nbnd))

    forall(j = 1:nbnd, i = 1:npw) buffer(i, j) = opt_basis%elements(i, j) * g2kin(i)

    call block_inner(nbnd, npw, opt_basis%elements, npw, buffer, npw, ham%con(:,:,1), ham%desc)

    if (nspin == 2) ham%con(:,:,2) = ham%con(:,:,1)
    ham%kin_con = ham%con(:,:,1)

    deallocate(buffer)
    call stop_clock( '   part1' )
  endif

  if (.not. saved) then
  ! linear term; almost the same as above except gtmp instead of g2kin
    allocate(gtmp(nbnd, nbnd), buffer(npw, nbnd))
    do ixyz = 1, 3
      gtmp = cmplx(0.d0, kind=DP)
      forall(j = 1:nbnd, i = 1:npw) buffer(i, j) = opt_basis%elements(i, j) * g(ixyz, igk(i))*tpiba

      call zgemm('C', 'N', nbnd, nbnd, npw, one, &
                 opt_basis%elements(:,:), npw, &
                 buffer, npw, zero, &
                 gtmp, nbnd)

!      do j = 1, npw
!        call zher('U', nbnd, cmplx(g(ixyz, igk(j))*tpiba, kind=DP), &
!                  opt_basis%elements(j,:), 1, &
!                  gtmp, nbnd)
!      enddo
!      ham%lin(ixyz,:,:) = conjg(gtmp)
      ham%lin(ixyz,:,:) = gtmp
    enddo
    deallocate(gtmp, buffer)
    call mp_sum(ham%lin, intra_pool_comm)
  endif 
  ! quadratic term is analytic (\delta(k,k) |k|^2), so we just add it in later

  call stop_clock( '  Kinetic' )
  ! ===========================================================================
  ! local potential
  ! ===========================================================================
  call start_clock( '  local_potential' )
!  if (saved) then
  if (.false.) then
    write(*,*) "using aux basis size: ", decomp_size
    allocate(decomp1(dffts%nnr, decomp_size), decomp2(nbnd, decomp_size))
    do j = 1, decomp_size
      call get_buffer(decomp1(:,j), dffts%nnr, 128, j)
      call get_buffer(decomp2(:,j), nbnd, 129, j)
    enddo
    allocate(buffer(dffts%nnr, decomp_size))
    forall (j = 1:decomp_size, i = 1:dffts%nnr) buffer(i,j) = decomp1(i,j) * V_rs(i, 1)
    allocate(gtmp(decomp_size, decomp_size))
    call ZGEMM('C', 'N', decomp_size, decomp_size, dffts%nnr, cmplx(1./dffts%nnr, kind=DP), &
               decomp1, dffts%nnr, &
               buffer, dffts%nnr, zero, &
               gtmp, decomp_size)
    deallocate(buffer, decomp1)
    allocate(buffer(nbnd, decomp_size))
    call ZGEMM('N', 'N', nbnd, decomp_size, decomp_size, one, &
               decomp2, nbnd, &
               gtmp, decomp_size, zero, &
               buffer, nbnd)
!    deallocate(gtmp); allocate(gtmp(nbnd,nbnd))
    call ZGEMM('N', 'C', nbnd, nbnd, decomp_size, one, &
               buffer, nbnd, &
               decomp2, nbnd, one, &
               ham%con(:,:,1), nbnd)
!    call ZGEMM('N', 'N', nbnd, nbnd, decomp_size, one, &
!                decomp2, nbnd, &
!                buffer, decomp_size, one, &
!                ham%con(:,:,1), nbnd)
    deallocate(buffer, decomp2, gtmp)
  else
#define BLOCK 32
  allocate(buffer(size(psic, 1), BLOCK), gtmp(npw, BLOCK))
  spin: do s = 1, nspin
  do j = 1, nbnd, BLOCK
    !
    nbl = min(BLOCK, nbnd - j + 1)
    buffer(1:dffts%nnr, 1:nbl) = ( 0.D0, 0.D0 )
    forall (i = 1:nbl) buffer(nls(igk(1:npw)), i) = opt_basis%elements(1:npw,j+i-1)
    do i = 1, nbl
      CALL invfft ('Wave', buffer(:,i), dffts)
    enddo
    ! ... product with the potential V_rs = (vltot+v) on the smooth grid
    forall (i = 1:nbl, k = 1:dffts%nnr)
      buffer(k,i) = buffer(k,i) * V_rs(k,s)
    end forall
    ! ... back to reciprocal space
    do i = 1, nbl
      CALL fwfft ('Wave', buffer(:,i), dffts)
    enddo
    ! ... store with correct ordering
    forall (i = 1:nbl) gtmp(1:npw, i) = buffer(nls(igk(1:npw)), i)
    !
    call ZGEMM('C', 'N', j+nbl-1, nbl, npw, one, &
                opt_basis%elements(:,:), npw, &
                gtmp, npw, one, &
                ham%con(:,j,s), nbnd)
  enddo
  enddo spin
  deallocate(buffer, gtmp)
  endif
  deallocate(igk, g2kin)

  call mp_sum( ham%con, intra_pool_comm)

  call stop_clock( '  local_potential' )

END SUBROUTINE build_h_coeff
