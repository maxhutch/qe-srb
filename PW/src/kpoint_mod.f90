!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE kpoint
  !
  ! ... Basic variables for k-point generations, as read from input
  !
  USE kinds,      ONLY : DP
  USE parameters,         ONLY : npk

  ! types!
  type k_list
    integer :: n
    integer :: ntot
    integer :: nred
    REAL(DP), pointer :: x(:,:)
    REAL(DP), pointer :: xr(:,:)
    REAL(DP), pointer :: w(:)
    REAL(DP), pointer :: wr(:)
    integer, pointer :: equiv(:)
  end type k_list
  !
  ! ... uniform k-point grid parameters
  !
  CONTAINS
 
  SUBROUTINE generate_k (style, nk1, nk2, nk3, k1, k2, k3, nrot, s, t_rev, at, bg, nkstot, xk, wk, kpoints ) 
    !
    ! initialize the grid of k points 
    !
    character(len=80), intent(IN) :: style
    INTEGER, INTENT (IN) :: nk1, nk2, nk3, k1, k2, k3 ! size of generated grid
    integer, intent(in) :: nrot, t_rev(48), s(3,3,48)
    real(DP), intent(in) :: at(3,3), bg(3,3)
    integer, intent(in) :: nkstot
    real(DP), intent(in) :: xk(:,:), wk(:)
    type(k_list), INTENT (inout) :: kpoints

    ! locals
    real(DP), parameter :: eps=1.0d-5
    logical, parameter :: time_reversal = .true.
    real(DP) :: xkr(3), xx, yy, zz, fact
    integer, allocatable :: map(:), equiv_tmp(:)
    integer :: i, j, k, n, nk, ns, stat
    logical :: in_list


    select case ( style ) 
    case ('automatic')    
      kpoints%ntot = nk1 * nk2 * nk3
      kpoints%n = kpoints%ntot
      deallocate(kpoints%x, kpoints%xr, kpoints%w, kpoints%wr, equiv_tmp, stat = stat) ! stat catches unallocated
      allocate(kpoints%x(3,kpoints%n), kpoints%w(kpoints%n), equiv_tmp(kpoints%ntot))
 
      do i = 0, nk1-1
        do j = 0, nk2-1
          do k = 0, nk3-1
            kpoints%x(1,i*nk3*nk2 + j * nk3 + k + 1) = dble(k1)/(2.*nk1) + dble(i)/nk1
            kpoints%x(2,i*nk3*nk2 + j * nk3 + k + 1) = dble(k2)/(2.*nk2) + dble(j)/nk2
            kpoints%x(3,i*nk3*nk2 + j * nk3 + k + 1) = dble(k3)/(2.*nk3) + dble(k)/nk3
          enddo
        enddo
      enddo
      kpoints%w = 2.d0 / kpoints%ntot
    !
    forall( nk = 1: kpoints%ntot ) equiv_tmp(nk) = nk
    kpoints%nred = kpoints%ntot
#if 1
    DO nk=1,kpoints%ntot
    !  check if this k-point has already been found equivalent to another
      IF (equiv_tmp(nk) == nk) THEN
        !  check if there are equivalent k-point to this in the list
        !  (excepted those previously found to be equivalent to another)
        !  check both k and -k
!    write(*,*) "boring", nrot  
        DO ns=1,nrot
           xkr = matmul(s(:,:,ns), kpoints%x(:,nk))
           xkr(:) = xkr(:) - nint( xkr(:) )
           IF(t_rev(ns)==1) xkr = -xkr
           xx = xkr(1)*nk1 - 0.5d0*k1
           yy = xkr(2)*nk2 - 0.5d0*k2
           zz = xkr(3)*nk3 - 0.5d0*k3
           in_list = abs(xx-nint(xx))<=eps .and. &
                         abs(yy-nint(yy))<=eps .and. &
                         abs(zz-nint(zz))<=eps
           IF (in_list) THEN
              i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
              j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
              k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              IF (n>nk .and. equiv_tmp(n)==n) THEN
                 equiv_tmp(n) = nk
                 kpoints%nred = kpoints%nred - 1 
              ELSE
                 IF (equiv_tmp(n)/=nk .or. n<nk ) CALL errore('generate_k', &
                    'something wrong in the checking algorithm',1)
              ENDIF
           ENDIF    
           IF ( time_reversal ) THEN
              xx =-xkr(1)*nk1 - 0.5d0*k1
              yy =-xkr(2)*nk2 - 0.5d0*k2
              zz =-xkr(3)*nk3 - 0.5d0*k3
              in_list=abs(xx-nint(xx))<=eps.and.abs(yy-nint(yy))<=eps &
                                                 .and. abs(zz-nint(zz))<=eps
              IF (in_list) THEN
                 i = mod ( nint (-xkr(1)*nk1 - 0.5d0 * k1 + 2*nk1), nk1 ) + 1
                 j = mod ( nint (-xkr(2)*nk2 - 0.5d0 * k2 + 2*nk2), nk2 ) + 1
                 k = mod ( nint (-xkr(3)*nk3 - 0.5d0 * k3 + 2*nk3), nk3 ) + 1
                 n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                 IF (n>nk .and. equiv_tmp(n)==n) THEN
                    equiv_tmp(n) = nk
                    kpoints%nred = kpoints%nred - 1
                 ELSE
                    IF (equiv_tmp(n)/=nk.or.n<nk) CALL errore('generate_k', &
                    'something wrong in the checking algorithm',2)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
      ENDIF
    ENDDO
#endif

    allocate(kpoints%xr(3, npk), kpoints%wr(npk))
    allocate(map(kpoints%ntot), kpoints%equiv(kpoints%ntot))
    kpoints%wr = 0.; k = 0
    !  count irreducible points and order them
    DO nk=1,kpoints%ntot
       IF (equiv_tmp(nk)==nk) THEN
         k = k + 1
         map(nk) = k
         kpoints%xr(:,k) = kpoints%x(:,nk) - floor(kpoints%x(:,nk))
       ENDIF
       kpoints%equiv(nk) = map(equiv_tmp(nk))
       kpoints%wr(kpoints%equiv(nk)) = kpoints%wr(kpoints%equiv(nk)) + 1
    ENDDO
    kpoints%wr = 2.d0 * kpoints%wr / sum(kpoints%wr)
    !
    case('crystal')
      kpoints%ntot = nkstot
      kpoints%n = kpoints%ntot
      kpoints%nred = kpoints%ntot
      deallocate(kpoints%x, kpoints%w, stat = stat) ! stat catches unallocated
      allocate(kpoints%x(3,kpoints%n), kpoints%w(kpoints%n))
      allocate(kpoints%xr(3,npk), kpoints%wr(npk))
      
      kpoints%x = xk
      kpoints%w = wk
      kpoints%xr = xk
      kpoints%wr = wk

    case('crystal_b')
      kpoints%ntot = sum(wk(1:nkstot-1)) + 1
      kpoints%n = kpoints%ntot
      kpoints%nred = kpoints%ntot
      deallocate(kpoints%x, kpoints%w, stat = stat) 
      allocate(kpoints%x(3,kpoints%n), kpoints%w(kpoints%n))
      allocate(kpoints%xr(3,npk), kpoints%wr(npk))
      
      j = 1
      do k = 1, nkstot - 1
        do i = 0, wk(k) - 1
          kpoints%x(:,j) = xk(:,k) + i * ((xk(:,k+1) - xk(:,k))/wk(k))
          j = j + 1
        enddo
      enddo
      kpoints%x(:,j) = xk(:,nkstot)
      kpoints%w(:) = 0.d0
      
      kpoints%xr = kpoints%x
      kpoints%wr = kpoints%w

    case('tpiba')
      kpoints%ntot = nkstot
      kpoints%n = kpoints%ntot
      kpoints%nred = kpoints%ntot
      deallocate(kpoints%x, kpoints%w, stat = stat) ! stat catches unallocated
      allocate(kpoints%x(3,kpoints%n), kpoints%w(kpoints%n))
      allocate(kpoints%xr(3,npk), kpoints%wr(npk))

      do i = 1, kpoints%n        
        kpoints%x(:,i) = matmul( transpose(at), xk(:,i) ) 
      enddo
      
      kpoints%w = 2* wk / sum(wk)
      kpoints%xr(:,1:kpoints%n) = kpoints%x
      kpoints%wr(1:kpoints%n) = kpoints%w

    case('tpiba_b')
      kpoints%ntot = sum(wk(1:nkstot-1)) + 1
      kpoints%n = kpoints%ntot
      kpoints%nred = kpoints%ntot
      deallocate(kpoints%x, kpoints%w, stat = stat) 
      allocate(kpoints%x(3,kpoints%n), kpoints%w(kpoints%n))
      allocate(kpoints%xr(3,npk), kpoints%wr(npk))
      
      j = 1
      do k = 1, nkstot - 1
        do i = 0, wk(k) - 1
          kpoints%x(:,j) = xk(:,k) + i * ((xk(:,k+1) - xk(:,k))/wk(k))
          j = j + 1
        enddo
      enddo
      kpoints%x(:,j) = xk(:,nkstot)

      do i = 1, kpoints%n        
        kpoints%x(:,i) = matmul( transpose(at), kpoints%x(:,i) ) 
      enddo

      kpoints%w(:) = 0.d0
      
      kpoints%xr = kpoints%x
      kpoints%wr = kpoints%w

 
    case default
      write(*,*) "kpoint style unhandled: ", style
    end select

  !  go to cartesian axis (in units 2pi/a0)
  ! CALL cryst_to_cart(kpoints%n,kpoints%x,bg,1)
  ! CALL cryst_to_cart(kpoints%nred,kpoints%xr,bg,1)
  ! normalize weights to one


    return
  END SUBROUTINE generate_k 

END MODULE kpoint

