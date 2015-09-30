!
!> Build reduced basis for representing wavefunctions
!
!#define FIO
! only enable if nbnd*nks > mp_sum(npw)
! #define USE_SVD
SUBROUTINE build_basis (evc_in, opt_basis, ecut_srb )
  USE kinds, ONLY : DP
  use constants, only : tpi
  USE mp_global, ONLY : me_pool, root_pool, nproc_pool, intra_pool_comm
  USE mp, only : mp_sum, mp_max, mp_barrier
  USE srb_types, ONLY : basis, kmap
  USE input_parameters, ONLY : ntrans, trace_tol, max_basis_size
  USE srb, only : srb_debug
  use srb_matrix, only : dmat, setup_dmat, block_inner, diag, pool_scope, serial_scope, print_dmat
  use srb_matrix, only : collect
  use buffers, only : open_buffer, save_buffer

  use lsda_mod, only : nspin
  USE cell_base, only : at, tpiba2, bg, omega, tpiba
  USE gvect, only : g, ngm, ig_l2g
  USE gvecs, only : nls
  USE klist, ONLY : xk_int => xk, ngk_int => ngk, nks_int => nks
  USE wvfct, only: ecutwfc_int => ecutwfc, npwx_int =>npwx
  use fft_base, only : dfftp_int => dfftp, dffts_int => dffts
  use mp_wave, only : mergewf, splitwf
  use wavefunctions_module, only : psic
  USE fft_interfaces, only: fwfft, invfft
  use fft_types , only : fft_dlay_descriptor
  use scalapack_mod, only : scalapack_distrib, scalapack_distrib_re, scalapack_localindex, scalapack_diag, scalapack_svd
  use buffers, only : get_buffer
  use io_files, only : iunwfc, nwordwfc
  use wavefunctions_module, only : evc_int=>evc
  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN) :: evc_in(:,:)
  TYPE(basis), INTENT(INOUT) :: opt_basis !>!< optimal basis
  real(DP), intent(out) :: ecut_srb
  ! Interfaces
  interface
  subroutine kpoint_expand(ndim, xk, at, xk_map )
  use kinds, only: DP
  use srb_types, only: kmap
  integer, intent(in) :: ndim(3)
  real(DP), intent(in) :: xk(:,:)
  real(DP), intent(in) :: at(3,3)
  type(kmap), pointer, intent(inout) :: xk_map(:)
  end subroutine
  SUBROUTINE s_psi( lda, n, m, psi, spsi )
  use kinds, only : DP
  USE noncollin_module, ONLY: npol
  INTEGER, INTENT(IN) :: lda, n, m
  COMPLEX(DP), INTENT(IN) :: psi(lda*npol,m)
  COMPLEX(DP), INTENT(OUT)::spsi(lda*npol,m)
  end subroutine s_psi
  subroutine gk_shuffle(gk_array, igk, npw, phase)
  use kinds, ONLY : DP
  COMPLEX(DP), intent(inout) :: gk_array(:)
  integer, intent(in) :: igk(:)
  integer, intent(in) :: npw
  REAL(DP), intent(in) :: phase(3)
  end subroutine gk_shuffle
  end interface
  COMPLEX(DP), external :: ZDOTC
  REAL(DP), external :: DZNRM2
  !
  ! locals
  !
  complex(dp),parameter :: iota=cmplx(0.d0,1.d0)
  COMPLEX(DP), allocatable :: evc_all(:,:,:), evc_tmp(:), evc(:,:), spsi(:)
  type(kmap), pointer :: xk_map(:)
  real(DP), allocatable :: xk(:,:), wk(:)
  real(DP), allocatable ::  xk_orig(:,:)
  real(DP) :: ecutwfc

  !> number of local plane waves [ngk(k)]
  integer               :: npw_gamma,        npw_orig,          npw_tmp
  !> max (across processors) number of local plane waves [mp_max(max(ngk))]
  integer               :: npwx_gamma,       npwx_orig,         npwx_tmp, ngk_gamma
  !> number of plane waves per kpoint
  integer,  allocatable ::                   ngk_orig(:),       ngk_tmp(:)
  !> maps k+G to G
  integer,  allocatable :: igk_gamma(:),     igk_orig(:,:),     igk_tmp(:)
  !> |k+G|
  real(DP), allocatable :: g2kin_gamma(:),   g2kin_orig(:,:),   g2kin_tmp(:)
  !> maps k+G to G (global)
  integer,  allocatable :: igk_l2g_gamma(:), igk_l2g_orig(:,:), igk_l2g_tmp(:,:)

  real(DP), allocatable :: xr(:,:)
  complex(DP), allocatable ::  expikr(:), ztmp3(:,:,:)
  complex(DP), allocatable :: S_g(:), S_g2(:,:), S_l(:,:), B_l(:,:), ztmp(:,:)
  real(DP), allocatable :: eigU(:)
  complex(DP) :: lprod, trace, ptrace, anti_trace
  real(DP) :: inv_norm
  integer :: nbasis, nbasis_trunc, nbasis_lr, nbasis_lc
  integer :: i, j, ij, k, k1, k2, l, iband, jband, i_l, j_l
!  integer :: a, b, c
!  integer :: x, y, z
  logical :: trivial, islocal
  ! shortcuts
  integer :: nks_orig, nks, nbnd, npw, igwx
  COMPLEX(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0,kind=DP)
  integer, allocatable, save :: shuffle(:,:)
  type(dmat) :: S, B, Z
  nbnd = size(evc_int,2)

  call start_clock('  mirrors')

  ! copy original k-point set information
  nks_orig = nks_int
  call n_plane_waves(ecutwfc_int, tpiba2, nks_int, xk_int, g, ngm, npwx_int, ngk_int) 
  allocate(xk_orig(3,nks_orig))
  xk_orig = xk_int(:,1:nks_orig)
  allocate( ngk_orig(nks_orig) )
  call n_plane_waves(ecutwfc_int, tpiba2, nks_orig, xk_orig, g, ngm, npwx_orig, ngk_orig) 
  allocate( igk_orig(npwx_orig,nks_orig) )
  allocate( g2kin_orig(npwx_orig,nks_orig) )
  allocate( igk_l2g_orig(npwx_orig,nks_orig) )
  igk_l2g_orig = 0
  do k=1,nks_orig
    call gk_sort (xk_orig(:,k), ngm, g, ecutwfc_int/tpiba2, npw_orig, igk_orig(:,k), g2kin_orig(:,k))
    igk_l2g_orig(1:npw_orig,k) = ig_l2g(igk_orig(1:npw_orig,k))
  enddo
  allocate(evc(npwx_orig, nbnd))


  ! apply translational symmetry to the list of kpoints
  nullify(xk_map)
  call kpoint_expand(ntrans, xk_orig, at, xk_map) 
  nks = sum(xk_map(:)%nmap)  
  ! make a new list of coodinates
  allocate( xk(3, nks), wk(nks) )
  j=0
  do k=1,nks_orig
    do i=1,xk_map(k)%nmap
      j = j + 1
      if (xk_map(k)%invert(i)) then 
        xk(:,j) = matmul( bg, xk_map(k)%xmap(:,i) ) - xk_orig(:,k)
      else
        xk(:,j) = matmul( bg, xk_map(k)%xmap(:,i) ) + xk_orig(:,k)
      endif
      wk(j) = 0.d0 ! I'm zeroing this just in case
    enddo
  enddo

  ! expand the cutoff to include all kpoints; resize ffts acordingly
  ! note that the process of resizing the fft's has tons of side effects
  ! Note: we're not actually doing this anymore
  ecutwfc = ecutwfc_int
  ngk_gamma = ngk_orig(1); npwx_gamma = npwx_orig
  allocate(igk_gamma(npwx_gamma), g2kin_gamma(npwx_gamma), igk_l2g_gamma(npwx_gamma))
  igk_gamma = igk_orig(:,1); g2kin_gamma = g2kin_orig(:,1); igk_l2g_gamma = igk_l2g_orig(:,1)
  ecut_srb = ecutwfc

  !==================
  ! construct the mirrors of the gamma point
  ! make space for the wavefunctions
  allocate(ngk_tmp(nks))
  call n_plane_waves(ecutwfc, tpiba2, nks, xk, g, ngm, npwx_tmp, ngk_tmp) ! get the sizes ngk_max_tmp and ngk_tmp(:)
  call mp_max(npwx_tmp, intra_pool_comm)
  if (nproc_pool == 1 .and. .not. allocated(shuffle)) then
    allocate(shuffle(npwx_tmp, nks))
    shuffle = -1
  endif
  allocate(igk_tmp(npwx_tmp), g2kin_tmp(npwx_tmp), igk_l2g_tmp(npwx_tmp, nks))
  allocate(evc_all(npwx_tmp, nbnd, nks))
  ! put real-space coordinates into a single dim array
  allocate(xr(dffts_int%nnr, 3), expikr(dffts_int%nnr))
  call rgrid(xr)

  j = 0
  k_point: do k = 1, nks_orig
    if (k == 1 .and. nspin == 1) then
      evc(1:ngk_orig(k),:) = evc_in(1:ngk_orig(k),:)
    else
      CALL get_buffer ( evc, nwordwfc, iunwfc, k )
    end if

    do i = 1, xk_map(k)%nmap
      j = j + 1 ! global kpoint index

      ! is this a trivial phase?
      trivial = .not. (any(abs(xk_map(k)%xmap(:,i)) > 0.00000001) .or. xk_map(k)%invert(i))

      ! make the phase
      ! e^(-i 2 pi xk_map.r ) - r and xk in crystal units
      if (.not. trivial .and. nproc_pool /= 0 ) then 
        expikr(:) = matmul( xr, xk_map(k)%xmap(:,i) )
        expikr(:) = exp( ( - iota * tpi ) * expikr(:) )
      end if

      ! make igk and igk_l2g
      npwx_int = npwx_tmp
      igk_tmp = 0
      call gk_sort(xk(:,j), ngm, g, ecutwfc / tpiba2, npw_tmp, igk_tmp, g2kin_tmp)
      npwx_int = npwx_orig
      igk_l2g_tmp(:, j) = 0
      igk_l2g_tmp(1:npw_tmp,j) = ig_l2g(igk_tmp(1:npw_tmp))
      igwx = max(maxval(igk_l2g_orig), maxval(igk_l2g_tmp(:,j)))
      call mp_max(igwx, intra_pool_comm)
      allocate(evc_tmp(igwx))

      band: do l = 1, nbnd
        iband = l

        ! reorder wavefunction by new global many-gamma ordering
        evc_all(:, iband, j) = cmplx(0.d0, kind=DP)
        if (.not. trivial ) then
          evc_tmp = cmplx(0.d0, kind=DP);
          call mergewf(evc(:, iband), evc_tmp, ngk_orig(k), igk_l2g_orig(:,k), &
                       me_pool, nproc_pool, root_pool, intra_pool_comm)
          call splitwf(evc_all(1:npw_tmp, iband, j), evc_tmp, npw_tmp,  igk_l2g_tmp(:,j), &
                       me_pool, nproc_pool, root_pool, intra_pool_comm)
        else
          evc_all(1:ngk_orig(k), iband, j) = evc(1:ngk_orig(k), iband)
        endif

        ! apply the phase shift in real space
        if (.not. trivial ) then
          ! if evc is on a single proc, just shuffle
          if ((.not. xk_map(k)%invert(i)) .and.  nproc_pool == 1) then
            if (shuffle(1, j) == -1) then
              psic = cmplx(0.d0, kind=DP)
              forall (i=1:npw_tmp) psic(nls(igk_tmp(i))) = i
              call invfft('Wave', psic, dffts_int)
              psic(1:dffts_int%nnr) = psic(1:dffts_int%nnr) * expikr(1:dffts_int%nnr) ! actual application
              call fwfft('Wave', psic, dffts_int)
              shuffle(1:npw_tmp, j) = abs(psic(nls(igk_tmp(1:npw_tmp)))) + .5
            endif
            psic(1) = cmplx(0.d0)
            psic(2:npw_tmp+1) = evc_all(1:npw_tmp, iband, j)
            evc_all(1:npw_tmp, iband, j) = psic(shuffle(1:npw_tmp, j)+1) 
          ! otherwise, do the full FFT
          else
            psic = cmplx(0.d0, kind=DP)
            psic(nls(igk_tmp(1:npw_tmp))) = evc_all(1:npw_tmp, iband, j)
            call invfft('Wave', psic, dffts_int)
            ! flip
            if (xk_map(k)%invert(i)) psic(1:dffts_int%nnr) = conjg(psic(1:dffts_int%nnr))  
            ! translate
            psic(1:dffts_int%nnr) = psic(1:dffts_int%nnr) * expikr(1:dffts_int%nnr) ! actual application
            call fwfft('Wave', psic, dffts_int)
            evc_all(1:npw_tmp, iband, j) = psic(nls(igk_tmp(1:npw_tmp)))
          endif
        endif

        ! reorder based on uniform gamma ordering
        if (any(abs(xk(:,j)) > 0.00000001)) then
          evc_tmp = cmplx(0.d0, kind=DP) 
          call mergewf(evc_all(1:npw_tmp, iband, j),      evc_tmp, npw_tmp,    igk_l2g_tmp(:,j), &
                       me_pool, nproc_pool, root_pool, intra_pool_comm)
          evc_all(:, iband, j) = cmplx(0.d0, kind=DP)
          call splitwf(evc_all(1:ngk_gamma, iband, j), evc_tmp, ngk_gamma,  igk_l2g_gamma,  &
                       me_pool, nproc_pool, root_pool, intra_pool_comm)
        endif

      enddo band
      deallocate(evc_tmp)

    if (j /= 1 .and. all(xk_orig(:,k) == xk_orig(:,1))) then
        lprod = sum(abs(evc_all(1:ngk_gamma,:,j) - evc_all(1:ngk_gamma,:,1)))
        if (abs(lprod) < 1.d-10) then
          nks = j - 1
          exit
        end if
    end if

    enddo
  enddo k_point
  deallocate(evc)
  deallocate(xr, expikr)

  call stop_clock('  mirrors')
  call start_clock('  make_S')
  !=====================
  ! constuct the overlap matrix
  nbasis = nbnd * nks ! global array size
  call setup_dmat(S, nbasis, nbasis, scope_in = pool_scope)
  call setup_dmat(B, nbasis, nbasis, scope_in = pool_scope)
  allocate(eigU(nbasis))

  call block_inner(nbnd*nks, ngk_gamma, &
                   one,  evc_all, npwx_tmp, &
                         evc_all, npwx_tmp, &
                   zero, S)

  call stop_clock('  make_S')

  !============================
  ! diagonalize it
  call start_clock('  diag_S')
  call diag(S, eigU, B)
  call stop_clock('  diag_S')
  !============================
  ! Truncate it
  if (max_basis_size < 0) max_basis_size = HUGE(max_basis_size)
#if 0
  open(unit=1234, file='vars', access='APPEND', status='replace')
  do i=1,nks*nbnd
  write(1234, *) i, eigU(i)
  enddo
  close(1234)
#endif
  if (trace_tol <= 1.d-16) then
    nbasis_trunc = min(nbasis, max_basis_size)
    ptrace = sum(abs(eigU(1:nbasis-nbasis_trunc))) / sum(abs(eigU))
  else
    trace = sum(abs(eigU))
    eigU = abs(eigU)/trace
    ptrace = eigU(1); nbasis_trunc = 1
    do while (abs(ptrace) < abs(trace_tol))
      nbasis_trunc = nbasis_trunc + 1
      ptrace = ptrace + eigU(nbasis_trunc)
    enddo
    nbasis_trunc = nbasis - nbasis_trunc
    nbasis_trunc = min(nbasis_trunc, max_basis_size)
    ptrace = sum(abs(eigU(1:nbasis-nbasis_trunc))) / sum(abs(eigU))
  endif
  if (me_pool == 0) then
    write(*,'(5X,A,I7,A,I7,A,F10.6,A)') "Made basis: ", nbasis_trunc, " of ", nks*nbnd, " elements (", abs(ptrace)*100, "% error)"
    if (nbasis_trunc / nks*nbnd > 0.9) then
      write(*,'(A,F5.1,A)') "WARNING! truncation removed only ", &
                            100.*(nks*nbnd - nbasis_trunc)/(nks*nbnd), "% of degrees of freedom."
      write(*,'(A)') "Consider increasing the number of q-points to improve statistics"
    endif
  endif
  !=============================
  ! Express the basis in the gamma-point pw basis <G|B>
  call start_clock('  expand')
  allocate(opt_basis%elements(ngk_gamma, nbasis_trunc)); opt_basis%length = nbasis_trunc

  call setup_dmat(Z, nbasis, nbasis_trunc, scope_in = serial_scope)
  call collect(B, Z, nbasis, nbasis_trunc, 1, nbasis-nbasis_trunc+1)

  ! transform <G|psi><psi|B>
  call ZGEMM('N', 'N', ngk_gamma, nbasis_trunc, nbasis, one, & 
             evc_all, npwx_tmp, &
             Z%dat, nbasis, zero, &
             opt_basis%elements, ngk_gamma)


  deallocate(eigU)
  deallocate(evc_all)

  ! Normalize 
  do i = 1, nbasis_trunc
    inv_norm = DZNRM2(ngk_gamma, opt_basis%elements(:,i), 1)
    inv_norm = inv_norm * inv_norm
    call mp_sum(inv_norm, intra_pool_comm)
    inv_norm = 1.d0 / sqrt(inv_norm)
    call zdscal(ngk_gamma, inv_norm, opt_basis%elements(:,i), 1)
  enddo
  call stop_clock('  expand')

  ! prepare to write to file
#ifdef FIO
  if (opt_basis%elem_rs_unit < -1) then
    opt_basis%elem_rs_unit = - opt_basis%elem_rs_unit
    call open_buffer(opt_basis%elem_rs_unit, 'elem_rs', dffts_int%nnr, 1, islocal)
  endif
#endif


  ! Bring to real-space if we'll need it later
  if (.false.) then
#ifdef FIO
    allocate(opt_basis%elem_rs(dffts_int%nnr, 1))
    do i = 1, nbasis_trunc
      psic(:) = ( 0.D0, 0.D0 )
      psic(nls(igk_gamma(1:ngk_gamma))) = opt_basis%elements(1:ngk_gamma,i)
      CALL invfft ('Wave', psic, dffts_int)
      opt_basis%elem_rs(:, 1) = psic(1:dffts_int%nnr)
      call save_buffer(opt_basis%elem_rs(:, 1), dffts_int%nnr, opt_basis%elem_rs_unit, i)
    enddo
    deallocate(opt_basis%elem_rs)
#else
    allocate(opt_basis%elem_rs(dffts_int%nnr, nbasis_trunc))
    do i = 1, nbasis_trunc
      psic(:) = ( 0.D0, 0.D0 )
      psic(nls(igk_gamma(1:ngk_gamma))) = opt_basis%elements(1:ngk_gamma,i)
      CALL invfft ('Wave', psic, dffts_int)
      opt_basis%elem_rs(:, i) = psic(1:dffts_int%nnr)
    enddo
#endif
  endif


  if (srb_debug) then
  do i = 1, nbasis_trunc
    do j = 1, nbasis_trunc
      lprod = dot_product(opt_basis%elements(:,i), opt_basis%elements(:,j))
      call mp_sum(lprod, intra_pool_comm)
      if (me_pool == 0 .and. i == j .and. abs(lprod - 1.d0) > .001) write(*,*) "WARN: <G|B> not orthonormal", i, j, lprod
      if (me_pool == 0 .and. i /= j .and. abs(lprod) > 0.001) write(*,*) "WARN: <G|B> not orthonormal", i, j, lprod
    enddo
  enddo
  endif

  deallocate(xk_orig, xk_map, xk, wk)! xk_orig)
  deallocate(ngk_orig, ngk_tmp, igk_gamma, igk_orig, igk_tmp)
  deallocate(g2kin_gamma, g2kin_orig, g2kin_tmp, igk_l2g_gamma, igk_l2g_orig, igk_l2g_tmp)

END SUBROUTINE build_basis

!==============================================================================
! Expand the kpoints using the translation group
!==============================================================================
subroutine kpoint_expand(ndim_in, xk, at, xk_map )
  use kinds, only: DP
  use srb_types, only: kmap

  IMPLICIT NONE

  integer, intent(in) :: ndim_in(3)
  real(DP), intent(in) :: xk(:,:)
  real(DP), intent(in) :: at(3,3)
  type(kmap), pointer, intent(inout) :: xk_map(:)

  real(dp),parameter :: eps=1.d-12
  real(dp), allocatable :: xk_latt(:,:)
  real(dp) :: xk_tmp(3)
  real(dp) :: disp
  integer :: ndim(3), nmin(3) = (/0,0,0/)
  integer :: izero(3), kg, jg, ig, ik, flip, counter
  integer :: nkstot
  real(dp), allocatable :: xk_all(:,:)

  ! now pad the BZ with edges
  ! and collect the final number of mapped kpoints

  ndim = abs(ndim_in)
  do ik=1,3
    if (ndim_in(ik) < 0) nmin(ik) = ndim_in(ik) 
  enddo

  nkstot = size(xk, 2) 
  allocate(xk_latt(3, nkstot))
  allocate(xk_all(3,10000))
  xk_all(:,:) = -10000000

  ! transform kpoint indicies to reciprical lattice coordinates
  do ik=1,nkstot
    xk_latt(:,ik) = matmul( transpose(at), xk(:,ik) )
  enddo

  ! init
  if( associated(xk_map) ) then
    do ik=1,size(xk_map)
      if (allocated(xk_map(ik)%xmap)) deallocate( xk_map(ik)%xmap )
    enddo
    deallocate(xk_map)
  endif
  allocate( xk_map(nkstot) )
  xk_map(:)%nmap = 0

    ! count translations
    do flip = 1,-1,-2
    do kg=nmin(3),ndim(3)
      do jg=nmin(2),ndim(2)
        do ig=nmin(1),ndim(1)
  do ik=1,nkstot
          xk_tmp = xk_latt(:,ik) * flip + (/ig, jg, kg/)
          if ((any(xk_tmp > ndim + 0.01)) &
             .or. (any(xk_tmp  < nmin - 0.01))) cycle
          counter = 1
          do while (all(xk_all(:,counter) > xk_all(:,10000)))
            if (all(abs(xk_tmp - xk_all(:,counter)) < eps)) exit
            counter = counter + 1
          enddo
          if (all(abs(xk_all(:,counter)  - xk_all(:,10000)) < eps)) then
            xk_all(:,counter) = xk_tmp
            xk_map(ik)%nmap = xk_map(ik)%nmap + 1
          endif 
  enddo
        enddo
      enddo
    enddo
    enddo

  ! allocate a map for each symmetry
  do ik=1,nkstot
    allocate( xk_map(ik)%xmap( 3, xk_map(ik)%nmap ) )
    allocate( xk_map(ik)%invert( xk_map(ik)%nmap ) )
  enddo

  ! store displacements
  xk_map(:)%nmap = 0
  xk_all(:,:) = -10000000
    do flip = 1,-1,-2
    do kg=nmin(3),ndim(3)
      do jg=nmin(2),ndim(2)
        do ig=nmin(1),ndim(1)
  do ik=1,nkstot
          xk_tmp = xk_latt(:,ik) * flip + (/ig, jg, kg/)
          if ((any(xk_tmp > ndim + 0.01)) &
             .or. (any(xk_tmp  < nmin - 0.01))) cycle
          counter = 1
          do while (all(xk_all(:,counter) > xk_all(:,10000)))
            if (all(abs(xk_tmp - xk_all(:,counter)) < eps)) exit
            counter = counter + 1
          enddo
          if (all(abs(xk_all(:,counter)  - xk_all(:,10000)) < eps)) then
            xk_all(:,counter) = xk_tmp
            xk_map(ik)%nmap = xk_map(ik)%nmap + 1
            xk_map(ik)%xmap(:,xk_map(ik)%nmap) = (/ ig, jg, kg /)
            xk_map(ik)%invert(xk_map(ik)%nmap) = (flip == -1)
!            write(*,'(A,3F,A,3F)') "Mirror ", xk_latt(:,ik), " to ", xk_tmp
          endif 
  enddo
        enddo
      enddo
    enddo
    enddo

end subroutine kpoint_expand

