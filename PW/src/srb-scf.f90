!----------------------------------------------------------------------------
SUBROUTINE srb_scf(evc, V_rs, rho, eband, demet, sc_error, skip)
!----------------------------------------------------------------------------
!
! Authors: Max Hutchinson, David Prendergast, PWSCF
! 
! ... calculates the symmetrized charge density using a srb interpolation scheme 
!
!
!#define SCUDA
#define DEBUG

  USE ISO_C_BINDING,        ONLY : c_ptr, C_NULL_PTR

  USE kinds,                ONLY : DP
  USE klist,                ONLY : nelec, lgauss, nks, nkstot, wk, xk
  USE control_flags,        ONLY : io_level 
  USE gvect,                ONLY : nl
  USE wavefunctions_module, ONLY : psic
  USE fft_base,             ONLY : dfftp
  USE gvecs,                ONLY : doublegrid
  USE funct,                ONLY : dft_is_meta
  USE io_global,            ONLY : stdout
  use uspp,                 only : nkb, okvan, becsum
  USE lsda_mod,             only : nspin, isk
  use wvfct,                only : wg, et

  USE srb_types,        ONLY : basis, ham_expansion, pseudop, nk_list, kproblem
  USE srb_matrix,       ONLY : grab_desc
  USE srb,              ONLY : qpoints, basis_life, freeze_basis
  USE srb,              ONLY : use_cuda, rho_reduced 
  use srb,              ONLY : states, bstates, wgq, red_basis=>scb, ets, pp=>spp
  USE srb, ONLY : build_basis, build_basis_cuda, build_h_coeff, build_h_matrix, diagonalize, build_rho
  USE srb, ONLY : build_h_coeff_cuda
  USE srb, ONLY : build_projs, build_projs_reduced, load_projs, copy_pseudo_cuda, build_projs_cuda, build_s_matrix, store_states
  use srb, only : solve_system_cuda
  use srb, only : build_rho_reduced, transform_rho, backload
  USE scf,                  ONLY : scf_type
  USE fft_types,            ONLY : fft_dlay_descriptor
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE symme,                ONLY : sym_rho
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm, me_pool, nproc_pool, me_image
  use scalapack_mod, only : scalapack_distrib
  use buffers, only : get_buffer

  !
  IMPLICIT NONE
  !
  ! ... arguments (yes, fortran allows arguments)
  !
  COMPLEX(DP),    INTENT(IN)  :: evc(:,:)   !>!< wavefunctions in PW basis
  REAL(DP),       INTENT(IN)  :: V_rs(:,:)  !>!< total potential in real space
  TYPE(scf_type), INTENT(INOUT) :: rho        !>!< electron density (and some other stuff)
  real(DP),       INTENT(OUT) :: eband      !>!< band contribution to energy
  real(DP),       INTENT(OUT) :: demet      !>!< ??? contribution to energy
  real(DP),       INTENT(IN) :: sc_error
  logical,        INTENT(OUT) :: skip
  ! Interfaces
  interface
    SUBROUTINE weights(nks, nkstot, wk, xk, et, wg, eband)
      USE kinds,                ONLY : DP
      integer, intent(in) :: nks, nkstot
      real(DP), intent(in) :: wk(:)
      real(DP), intent(in) :: xk(:,:)
      real(DP), intent(in) :: et(:,:)
      real(DP), intent(out) :: wg(:,:)
      real(DP), intent(out) :: eband
    end subroutine weights
    subroutine addusdens_g(bar, foo)
      use kinds, only : DP
      real(DP), intent(inout) :: bar(:,:,:)
      real(DP), intent(inout) :: foo(:,:)
    end subroutine addusdens_g
  end interface
  !
  ! ... local variables
  !
  TYPE(ham_expansion), save :: h_coeff !>!< Expansion of Hamiltonian wrt k
  TYPE(kproblem)            :: Hk !>!< Hamiltonain at a specific kpoint
  complex(DP), allocatable :: S_matrix2(:,:) !>!< Overlap matrix and copy
  real(DP), allocatable, target ::  energies(:,:) !>!< eigen-values (energies)
  complex(DP), allocatable :: evecs(:,:) !>!< eigenvectors of Hamiltonian
  complex(DP), allocatable :: P(:,:,:), Pinv(:,:,:) !>!< Preconditioner and inverse
  complex(DP), allocatable :: rho_srb(:,:,:) !>!< Denstiy matrix in reduced basis
  real(DP), allocatable :: wr2(:), xr2(:,:) !>!< copies of k-points and weights
  ! parameters
  integer, save :: basis_age = -1, itr_count = 0
  integer :: nbnd
  real(DP), save :: ecut_srb
  integer :: meth
  ! tmps and workspace
  complex(DP) :: ztmp
  INTEGER :: i, j, k, q, s
  complex(DP), pointer :: ptr(:,:) => NULL()
  integer, allocatable :: ipiv(:)
  complex(DP), allocatable :: work(:)

#ifdef SCUDA
  if (use_cuda) call setup_cuda()
  write(*,*) "setup cuda"
#endif

  ! Checks, checks, checks
  if ( dft_is_meta() ) then
    write(*,*) "Shirley for meta DFT not yet supported "
  end if

  call start_clock( 'srb')
  call start_clock(  ' other')

  ! inits
  nbnd = size(evc, 2)
  states%nbnd  = nbnd
  bstates%nbnd = nbnd

  !
  ! ... Do a basis calculation
  !
  if (basis_age >= basis_life - 1 .and. sc_error > freeze_basis) then
    basis_age = -1
    deallocate(red_basis%elements)
    if (allocated(red_basis%elem_rs)) deallocate(red_basis%elem_rs)
    deallocate(h_coeff%con, h_coeff%lin, h_coeff%kin_con)
    red_basis%length = -1
  endif

  call stop_clock(  ' other')
  if (basis_age == -1) then
    call start_clock( ' build_red_basis' )
    call build_basis(evc, red_basis, ecut_srb)
    call stop_clock( ' build_red_basis' )
  else
    if (me_image == 0) write(*,*) "Using a basis of length ", red_basis%length, " and age ", basis_age
  endif
  basis_age = basis_age + 1
  skip = ((basis_age < basis_life - 1) .or. sc_error <  freeze_Basis)

  !
  ! ... Construct the local Hamiltonian
  ! 
  call scalapack_distrib(red_basis%length, red_basis%length, h_coeff%desc%nrl, h_coeff%desc%ncl)
  call grab_desc(h_coeff%desc)

  call start_clock( ' build_h_coeff' )
  call build_h_coeff(red_basis, V_rs(:,:), ecut_srb, nspin, h_coeff, basis_age /= 0)
  call stop_clock( ' build_h_coeff' )
  call start_clock(  ' other')

  !
  ! ... Setup dense data structures
  !
  Hk%desc = h_coeff%desc
  allocate(Hk%H(Hk%desc%nrl, Hk%desc%ncl))
  allocate(evecs(red_basis%length, nbnd))
  if (allocated(energies)) deallocate(energies)
  allocate(energies(red_basis%length, nspin*qpoints%nred))
  energies = 0.d0

  states%nk = qpoints%nred * nspin
  allocate(pp%projs(red_basis%length, nkb))
  if (okvan) then
    allocate(Hk%S(Hk%desc%nrl,Hk%desc%ncl), S_matrix2(Hk%desc%nrl,Hk%desc%ncl))
    pp%us = .true.
    bstates%nk = qpoints%nred * nspin
  else
    bstates%nk = 0
    allocate(Hk%S(1,1), S_matrix2(1,1))
  end if

  ! Swapping is controled by the disk_io input param
  ! setting 3rd dim to size 1 signals code to swap each k-point
  if (associated(states%host_ar)) deallocate(states%host_ar)
  if (associated(bstates%host_ar)) deallocate(bstates%host_ar)
  if (io_level < 1) then
    allocate(states%host_ar(red_basis%length, states%nbnd, states%nk+nproc_pool))
    if (bstates%nk == 0) then 
      allocate(bstates%host_ar(nkb, bstates%nbnd, 0))
    else 
      allocate(bstates%host_ar(nkb, bstates%nbnd, bstates%nk+nproc_pool))
    endif
  else
    allocate(states%host_ar(red_basis%length, states%nbnd, 1))
    allocate(bstates%host_ar(nkb, bstates%nbnd, min(bstates%nk,1)))
  endif

  !
  ! ... Transform projectors
  !
  call stop_clock(  ' other')
  if (nkb > 0 .and. basis_age == 0) then
    call start_clock(' build_proj')
    if (.true.) then
      call build_projs_reduced(red_basis, qpoints%xr(:,:), qpoints%nred, pp)
    else
      call build_projs(red_basis, qpoints%xr(:,:), qpoints%nred, pp)
    endif
    call stop_clock(' build_proj')
  endif


  !
  ! ... main q-point loop
  !
  do q = 1+me_pool, qpoints%nred, nproc_pool
    !
    ! ... Build dense S matrix 
    !
    call start_clock(  ' other')
    if (nkb > 0) then
      call load_projs(q, pp)
    endif
    call stop_clock(  ' other')

    if (okvan) then
      CALL start_clock(' build_mat' )
      if (basis_age == 0) then
        call build_s_matrix(pp, (1-q)/nproc_pool - 1, Hk)
      else
        call build_s_matrix(pp, (q-1)/nproc_pool + 1, Hk)
      endif
      CALL stop_clock(' build_mat' )
    end if
#ifdef DAVID
    ! allocate space for preconditioner
    if ((q-1)/nproc_pool == 0) then
      allocate(P(red_basis%length,red_basis%length,nspin))
      allocate(Pinv(red_basis%length,red_basis%length,nspin))
      allocate(ipiv(red_basis%length), work(16*red_basis%length)) 
    endif
#endif

    ! loop over spins (which share S matrices
    do s = 1, nspin
        !
        ! ... Build dense Hamiltonian matrix
        !
        CALL start_clock(' build_mat' )
        call build_h_matrix(h_coeff, qpoints%xr(:,q), pp, s, Hk)
        CALL stop_clock(' build_mat' )

        ! 
        ! ... Diagonalize the dense Hamiltonian 
        !
        CALL start_clock( ' diagonalize' )
#ifdef DAVID
        ! invert preconditioner
        if ((q-1)/nproc_pool == 1) then
          P(:,:,s) = Hk%H 
          Pinv(:,:,s) = P(:,:,s) 
          call ZHETRF( 'U', red_basis%length, Pinv(:,:,s), red_basis%length, &
                      ipiv, work, 16*red_basis%length, i )
          call ZHETRI( 'U', red_basis%length, Pinv(:,:,s), opt_Basis%length, &
                      ipiv, work, i )
        endif

        ! try to load old states as starting guesses
        if (basis_age /= 0) then ! use last iteration's evecs
          if (size(states%host_ar, 3) == 1) then
            ptr => states%host_ar(:,:,1)
            call get_buffer(ptr, red_basis%length*nbnd, states%file_unit, &
                            (q+(s-1)*(qpoints%nred+nproc_pool)-1)/nproc_pool + 1)
          else
            ptr => states%host_ar(:,:,(q+(s-1)*(qpoints%nred+nproc_pool)-1)/nproc_pool + 1)
          endif
          evecs = ptr
          meth = 2
        else if (q /= 1+me_pool) then ! load the most recently computed q-point
          if (size(states%host_ar, 3) == 1) then
            ptr => states%host_ar(:,:,1)
            call get_buffer(ptr, red_basis%length*nbnd, states%file_unit, &
                            (q+(s-1)*(qpoints%nred+nproc_pool)-1)/nproc_pool)
          else
            ptr => states%host_ar(:,:,1)
          endif
          evecs = ptr
          meth = 2
        else ! fallback to exact diagonalization
          meth = 0
        endif
#endif
        if (okvan) then
          S_matrix2(:,:) = Hk%S(:,:)
          call diagonalize(Hk%H, energies(:,q+(s-1)*qpoints%nred), evecs, &
                           nbnd, S_matrix2, meth_opt = meth)
        else
          call diagonalize(Hk%H, energies(:,q+(s-1)*qpoints%nred), evecs, &
                           nbnd, meth_opt = meth)
        end if
        CALL stop_clock( ' diagonalize' )

        !
        ! ... Compute <\b|psi> and store it and <b|psi>
        !
        CALL start_clock( ' store' )
        call store_states(evecs, pp, (q+(s-1)*(qpoints%nred+nproc_pool)-1)/nproc_pool + 1, states, bstates)
        CALL stop_clock( ' store' )
    enddo 
  enddo 
  call start_clock(  ' other')
  call mp_sum(energies, intra_pool_comm)
  deallocate(Hk%H, Hk%S, S_matrix2, pp%projs, evecs)
#ifdef DAVID
  deallocate(P, Pinv)
  deallocate(work, ipiv)
#endif
  call stop_clock(  ' other')
  CALL start_clock( ' build_rho' )

  !
  ! ... make new weights
  !
  if (associated(wgq)) deallocate(wgq)
  nullify(wgq); allocate(wgq(red_basis%length, nspin*qpoints%nred))
  allocate(wr2(nspin * qpoints%nred), xr2(3, nspin*qpoints%nred))
  do s = 0, nspin-1
    wr2(1+s*qpoints%nred:qpoints%nred + s*qpoints%nred) = qpoints%wr(1:qpoints%nred) / nspin
    xr2(:,1+s*qpoints%nred:qpoints%nred + s*qpoints%nred) = qpoints%xr(:,1:qpoints%nred)
  enddo
  call weights(qpoints%nred*nspin, qpoints%nred*nspin, wr2, xr2, energies(1:nbnd,:), wgq(1:nbnd,:), demet)
  deallocate(wr2, xr2)

#ifdef DEBUG
  if( nspin==2 .and. me_image == 0 ) then
    write(*,*) "sum(weights): ", sum(wgq(1:nbnd,1:qpoints%nred)), sum(wgq(1:nbnd,qpoints%nred+1:2*qpoints%nred))
  else if( nspin==1 .and. me_image == 0 ) then
    write(*,*) "sum(weights): ", sum(wgq(1:nbnd,1:qpoints%nred))
  endif
#endif

  eband = sum(energies(1:nbnd, :) * wgq(1:nbnd, :)) 
  ets = energies(1:nbnd, :)

  !
  ! ...  Build the density <r|\rho|r>
  !
  rho%of_r(:,:) = 0.D0; rho%of_g(:,:) = 0.D0; becsum = 0.d0
  if (rho_reduced) then
    allocate(rho_srb(red_basis%length, red_basis%length, nspin))
    rho_srb = cmplx(0.d0, kind = DP)
    ! Build reduced density matrix <b|\rho|b'>
    call start_clock( '  reduced')
    call build_rho_reduced(states, bstates, wgq(:,:), qpoints%wr(:) / nspin, nspin, rho_srb, becsum)
    call stop_clock( '  reduced')
    ! Take SVD of \rho and transform to <r|\rho|r>
    call start_clock( '  trans')
    call transform_rho(rho_srb, red_basis, rho)
    call stop_clock( '  trans')
    deallocate(rho_srb)
  else
    ! Go straight to real-space (bad idea)
    call build_rho(states, bstates, wgq(:,:), qpoints%wr(:) / nspin, red_basis, nspin, rho, becsum)
  endif
  CALL stop_clock( ' build_rho' )
  call start_clock(  ' other')

  ! ... interpolate rho(r) if needed
  if (doublegrid) then
    do s = 1, nspin
      CALL interpolate(rho%of_r(1,s), rho%of_r(1,s), 1)
    enddo
  endif

  ! ... add ultra-soft correction
  write(*,*) shape(becsum), becsum
  if (okvan)  call addusdens_g(becsum, rho%of_r)

  ! ... symmetrize rho(G) 
  do s = 1, nspin
    psic(:) = rho%of_r(:,s)
    CALL fwfft ('Dense', psic, dfftp)
    rho%of_g(:,s) = psic(nl(:))
  enddo
  CALL sym_rho ( nspin, rho%of_g )
  do s = 1, nspin
    psic(:) = ( 0.D0, 0.D0 )
    psic(nl(:)) = rho%of_g(:,s)
    CALL invfft ('Dense', psic, dfftp)
    rho%of_r(:,s) = psic(:)
  enddo

  itr_count = itr_count + 1

!  if (.not. skip) then
!    call backload(qpoints, red_basis, states)
!  endif

  call stop_clock(  ' other')
  call stop_clock('srb')

  return

end subroutine srb_scf

