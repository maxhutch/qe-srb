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
  use srb,              ONLY : red_basis=>scb
  use uspp, only : nkb, okvan
  use uspp_param, only : nhm
  use ions_base, only : nat
  use lsda_mod, only : nspin
  use cell_base, only : bg

  use srb_types, only : basis, ham_expansion, pseudop, kproblem
  use srb, only : build_basis, build_h_coeff, build_h_matrix, diagonalize
  use srb, only : build_projs_reduced, load_projs, build_s_matrix
  use srb_matrix, only : setup_dmat, dmat, pot_scope, copy_dmat
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
!  TYPE(basis) :: red_basis  !>!< optimal basis; looks like wavefunctions
  TYPE(ham_expansion) :: h_coeff
  TYPE(pseudop) :: projs
  TYPE(kproblem) :: Hk
  complex(DP), allocatable :: S_matrix2(:,:)
  type(dmat) :: evecs
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

  call start_clock( ' build_red_basis' )
  ! Do a basis calculation
  call build_basis(evc, red_basis, ecut_srb)
  call stop_clock( ' build_red_basis' )

  ! Construct the Hamiltonian
  ! set up the blacs descriptors
  if (.not. allocated(h_coeff%con)) allocate(h_coeff%con(nspin))
  do s = 1,nspin
    call setup_dmat(h_coeff%con(s), red_basis%length, red_basis%length)
  enddo
  call setup_dmat(h_coeff%kin_con, red_basis%length, red_basis%length)

  call start_clock( ' build_h_coeff' )
  write(*,*) "making coeffs"
  call build_h_coeff(red_basis, V_rs, ecut_srb, nspin, h_coeff)
  write(*,*) "made coeffs"
  call stop_clock(' build_h_coeff')

  ! Diagonalize the Hamiltonian, producing states
  call setup_dmat(evecs, red_basis%length, nbnd, red_basis%length,min(16,nbnd),pot_scope)
  allocate(energies(red_basis%length, nspin*qpoints%nred))
  call setup_dmat(Hk%H, red_basis%length, red_basis%length)
  if (okvan) then
    call copy_dmat(Hk%S, Hk%H)
    allocate(S_matrix2(size(Hk%S%dat, 1), size(Hk%S%dat,2)))
    projs%us = .true.
    Hk%generalized = .true.
  else
    allocate(S_matrix2(1,1))
  end if

  if (nkb > 0) then
    call build_projs_reduced(red_basis, qpoints%xr(:,:), qpoints%nred, projs)
  endif

  energies = 0. 
  do q = 1+me_pool, qpoints%nred, nproc_pool
    if (nkb > 0) call load_projs(q, projs)

    if (okvan) then
      call start_clock(' build_proj')
      call build_s_matrix(projs, (1-q)/nproc_pool - 1, Hk)
      call stop_clock(' build_proj')
    end if
    spin: do s = 1, nspin
      CALL start_clock(' build_mat' )
      call build_h_matrix(h_coeff, qpoints%xr(:,q), projs, s, Hk)
      CALL stop_clock(' build_mat' )

      CALL start_clock( ' diagonalize' )
        if (okvan) then
          Hk%generalized = .true.
          S_matrix2(:,:) = Hk%S%dat
        else
          Hk%generalized = .false.
        end if
        call diagonalize(Hk, energies(:,q+(s-1)*qpoints%nred), evecs, &
                         nbnd, meth_opt = meth)
        if (Hk%generalized) Hk%S%dat = S_matrix2
      CALL stop_clock( ' diagonalize' )
    enddo spin
  enddo
  call mp_sum(energies, intra_pool_comm)
  deallocate(h_coeff%lin, S_matrix2)
  !
  WRITE( stdout, 9000 ) get_clock( 'PWSCF' )
  !
  WRITE( stdout, 9102 )

  open(unit=66, file="bands.dat")
  write(fmtstr,'(a,i8,a)') '(i8,3e20.12,',red_basis%length,'e20.12)'
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

