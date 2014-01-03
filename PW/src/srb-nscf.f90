!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
#define __SSDIAG
!-----------------------------------------------------------------------
  SUBROUTINE srb_nscf ( V_rs)
  !-----------------------------------------------------------------------
  !
  ! ... diagonalization of the KS hamiltonian in the non-scf case
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE bp,                   ONLY : lelfield, lberry, lorbm
  USE check_stop,           ONLY : stopped_by_user
  USE control_flags,        ONLY : io_level, conv_elec
  USE ener,                 ONLY : ef
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : iunwfc, nwordwfc, iunefield
  USE buffers,              ONLY : save_buffer
  USE klist,                ONLY : xk, wk, nks, nkstot
  USE lsda_mod,             ONLY : lsda, nspin
  USE wvfct,                ONLY : nbnd, et, npwx
  USE wavefunctions_module, ONLY : evc

  use srb, only : qpoints
  use srb,              ONLY : opt_basis=>scb
  use uspp, only : nkb, okvan
  use uspp_param, only : nhm
  use ions_base, only : nat
  use lsda_mod, only : nspin
  use cell_base, only : bg

  use srb_types, only : basis, ham_expansion, pseudop, kproblem
  use srb, only : build_basis, build_h_coeff, build_h_matrix, diagonalize
  use srb, only : build_projs, build_projs_reduced, load_projs, build_s_matrix
  use mp, only : mp_sum
  use mp_global, only : intra_pool_comm, me_pool, nproc_pool
  use scalapack_mod, only : scalapack_distrib

  !
  IMPLICIT NONE
  !
  REAL(DP), intent(in) :: V_rs(:,:)
  !
  ! ... local variables
  !
  REAL(DP), EXTERNAL :: get_clock
  INTEGER :: iter = 1, i, ik, q, j, s, meth = 0
  REAL(DP) :: dr2 = 0.d0
!  TYPE(basis) :: opt_basis  !>!< optimal basis; looks like wavefunctions
  TYPE(ham_expansion) :: ham
  TYPE(pseudop) :: projs
  TYPE(kproblem) :: Hk
  complex(DP), allocatable :: S_matrix2(:,:), evecs(:,:)
  real(DP), allocatable ::  energies(:,:)
  real(DP) :: ecut_srb
  character(255) :: fmtstr
  !
  !
  CALL start_clock( 'electrons' )
  iter = 1
  !
  WRITE( stdout, 9002 )
  CALL flush_unit( stdout )
  !
  IF ( lelfield) THEN
     !
     CALL c_bands_efield ( iter )
     !
  ELSE
     !
     CALL c_bands_nscf ( )
     !
  END IF

  call start_clock( 'srb')

  call start_clock( ' build_opt_basis' )
  ! Do a basis calculation
  call build_basis(evc, opt_basis, ecut_srb)
  call stop_clock( ' build_opt_basis' )

  ! Construct the Hamiltonian
  call start_clock( ' build_h_coeff' )
  write(*,*) "making coeffs"
  call build_h_coeff(opt_basis, V_rs, ecut_srb, nspin, ham)
  write(*,*) "made coeffs"
  call stop_clock(' build_h_coeff')

  ! Diagonalize the Hamiltonian, producing states
#ifdef __SSDIAG
  allocate(Hk%H(opt_basis%length, opt_basis%length))
  allocate(projs%projs(opt_basis%length, nkb))
  if (okvan) then
    allocate(Hk%S(opt_basis%length, opt_basis%length), S_matrix2(opt_basis%length, opt_basis%length))
    projs%us = .true.
  else
    allocate(Hk%S(1,1), S_matrix2(1,1))
  end if
#else
  call scalapack_distrib(opt_basis%length, opt_basis%length, i, j)
  allocate(Hk%H(i, j))
  if (okvan) then
    allocate(Hk%S(i, j), S_matrix2(i,j))
    allocate(projs%projs(opt_basis%length, nkb))
    projs%us = .true.
  else
    allocate(Hk%S(1,1), S_matrix2(1,1), projs%projs(1,1))
  end if
#endif
  allocate(evecs(opt_basis%length, nbnd))
  allocate(energies(opt_basis%length, nspin*qpoints%nred))


  if (nkb > 0) then
    call build_projs_reduced(opt_basis, qpoints%xr(:,:), qpoints%nred, projs)
  endif

  energies = 0. 
  do q = 1+me_pool, qpoints%nred, nproc_pool
    write(*,*) q
    if (nkb > 0) call load_projs(q, projs)
    write(*,*) q

    if (okvan) then
      call start_clock(' build_proj')
      call build_s_matrix(projs, (1-q)/nproc_pool - 1, Hk)
      call stop_clock(' build_proj')
    end if
    write(*,*) q
    spin: do s = 1, nspin
      CALL start_clock(' build_mat' )
      call build_h_matrix(ham, qpoints%xr(:,q), projs, s, Hk%H)
      CALL stop_clock(' build_mat' )

    write(*,*) q
      CALL start_clock( ' diagonalize' )
        if (okvan) then
          S_matrix2(:,:) = Hk%S(:,:)
          call diagonalize(Hk%H, energies(:,q+(s-1)*qpoints%nred), evecs, & 
                           nbnd, S_matrix2, meth_opt = meth)
        else
          call diagonalize(Hk%H, energies(:,q+(s-1)*qpoints%nred), evecs, &
                           nbnd, meth_opt = meth)
        end if
      CALL stop_clock( ' diagonalize' )
    enddo spin
    write(*,*) q
  enddo
  call mp_sum(energies, intra_pool_comm)
  deallocate(Hk%H, ham%con, ham%lin, Hk%S, S_matrix2, projs%projs)
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  WRITE( stdout, 9102 )

  open(unit=66, file="bands.dat")
  write(fmtstr,'(a,i8,a)') '(i8,3e20.12,',opt_basis%length,'e20.12)'
  do q = 1, nspin * qpoints%nred
    write(66,fmtstr) q, qpoints%xr(:,MOD(q-1,qpoints%nred)+1), energies(:,q) * rytoev
  enddo
  close(66)
 !
  CALL stop_clock( 'electrons' )
  !
9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9002 FORMAT(/'     Band Structure Calculation' )
9102 FORMAT(/'     End of band structure calculation' )
  !
END SUBROUTINE srb_nscf

