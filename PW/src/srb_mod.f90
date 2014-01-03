!
! Copyright (C) 2012 LBL
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
module srb_types
  USE ISO_C_BINDING, only : c_ptr, C_NULL_PTR
  USE kinds, only : DP
  USE srb_matrix, only : mydesc

  TYPE ham_expansion
     COMPLEX(DP), ALLOCATABLE :: con(:,:,:) ! H(i,j,k) is <i|H^(0)|j> (spin k)
     COMPLEX(DP), ALLOCATABLE :: lin(:,:,:) ! H(i,j) is <i|H^(0)|j>
     COMPLEX(DP), ALLOCATABLE :: kin_con(:,:) ! constant part of kinetic energy
     TYPE(mydesc)             :: desc
     TYPE(c_ptr)              :: con_d = C_NULL_PTR
     TYPE(c_ptr)              :: lin_d = C_NULL_PTR
     INTEGER                  :: length
  END TYPE ham_expansion
  !
  TYPE kproblem
     COMPLEX(DP), allocatable :: H(:,:)
     COMPLEX(DP), allocatable :: S(:,:)
     TYPE(mydesc)             :: desc
     integer                  :: k 
  END TYPE kproblem
  !
  TYPE basis
     COMPLEX(DP), ALLOCATABLE :: elements(:,:) ! basis_elements(i,j) is <B_j|psi_i>
     TYPE(c_ptr)              :: elements_d = C_NULL_PTR
     COMPLEX(DP), ALLOCATABLE :: elem_rs(:,:) ! <R|psi>
     INTEGER                  :: elem_gs_unit = -4471
     INTEGER                  :: elem_rs_unit = -7225
     INTEGER                  :: length ! size of basis
  END TYPE basis
  !
  TYPE pseudop
    COMPLEX(DP), ALLOCATABLE :: projs(:,:)
    TYPE(c_ptr)              :: projs_d = C_NULL_PTR
    TYPE(c_ptr)              :: Q_d = C_NULL_PTR
    TYPE(c_ptr)              :: D_d = C_NULL_PTR
    TYPE(c_ptr)              :: S_d = C_NULL_PTR
    INTEGER                  :: projs_unit = -3997
    INTEGER                  :: s_unit = -3998
    integer                  :: b_unit = -3999
    integer                  :: b_size(16)
    INTEGER                  :: nkb
    logical                  :: us = .false.
  END TYPE pseudop
  !
  TYPE nk_list
    COMPLEX(DP), pointer :: host_ar(:,:,:)
    TYPE(c_ptr)          :: device_ptr = C_NULL_PTR
    INTEGER              :: file_unit = -1
    integer              :: nbnd, nk
  END TYPE nk_list
  !
  type kmap
    integer               :: nmap
    real(dp), allocatable :: xmap(:,:)
  end type kmap
end module srb_types

MODULE srb
  !
  ! ... srb contains data structures and data for srb interpolation in SCF
  !
  USE kinds,            ONLY : DP
  USE srb_types, ONLY : ham_expansion, basis, pseudop, nk_list
  USE kpoint, only : k_list
  !
  INTEGER, allocatable :: band_index(:)
  integer, save :: basis_life
  real(DP), save :: freeze_basis
  logical :: rho_reduced
  logical :: srb_debug
  logical :: use_cuda
  REAL(DP) :: trkin2nlp(3,3), trnlp2kin(3,3)

  type(k_list), save :: qpoints
  REAL(DP), pointer, save :: wgq(:,:) => NULL()
  real(DP), allocatable, save :: ets(:,:)
  type(basis),   save :: scb  
  type(pseudop), save :: spp 
  type(nk_list), save :: states
  type(nk_list), save :: bstates
  integer, save :: decomp_size = -1

  ! interfaces
  interface
     SUBROUTINE srb_scf(evc, V_rs, rho, eband, demet, sc_error, skip)
        use kinds, only : DP
        use scf, only : scf_type
        COMPLEX(DP),    INTENT(IN)  :: evc(:,:)   !>!< wavefunctions in PW basis
        REAL(DP),    INTENT(IN)  :: V_rs(:,:) !>!< total potential in real space
        TYPE(scf_type), INTENT(INOUT) :: rho
        real(DP),       INTENT(OUT) :: eband
        real(DP),       INTENT(OUT) :: demet
        real(DP),       INTENT(in) :: sc_error
        logical,        intent(out) :: skip
     end subroutine srb_scf

     SUBROUTINE srb_nscf(V_rs)
        use kinds, only : DP
        REAL(DP),    INTENT(IN)  :: V_rs(:,:) !>!< total potential in real space
     end subroutine srb_nscf

    SUBROUTINE build_basis (evc, opt_basis, ecut_srb )
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis
      USE fft_types, ONLY : fft_dlay_descriptor
      COMPLEX(DP), INTENT(IN) :: evc(:,:)
      TYPE(basis), INTENT(OUT) :: opt_basis !>!< optimal basis
      real(DP), intent(out) :: ecut_srb
    end subroutine build_basis
  
    SUBROUTINE build_basis_cuda(evc, opt_basis, ecut_srb )
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis
      USE fft_types, ONLY : fft_dlay_descriptor
      COMPLEX(DP), INTENT(IN) :: evc(:,:)
      TYPE(basis), INTENT(OUT) :: opt_basis !>!< optimal basis
      real(DP), intent(out) :: ecut_srb
    end subroutine build_basis_cuda

    SUBROUTINE build_h_coeff(opt_basis, V_rs, ecut_srb, nspin, ham, saved)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis, ham_expansion
      TYPE(basis),         INTENT(IN)  :: opt_basis        !>!< optimal basis
      REAL(DP),         INTENT(IN)  :: V_rs(:,:)         !>!< total potential in real space
      REAL(DP),            INTENT(IN)  :: ecut_srb !>!< energy cutoff to include Gamma^8
      INTEGER,             INTENT(IN)  :: nspin
      TYPE(ham_expansion), INTENT(inout) :: ham              !>!< Hamiltonian
      logical, intent(in), optional :: saved
    end subroutine build_h_coeff

    SUBROUTINE build_h_coeff_cuda(opt_basis, V_rs, ecut_srb, nspin, ham, saved)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis, ham_expansion
      TYPE(basis),         INTENT(IN)  :: opt_basis        !>!< optimal basis
      REAL(DP),         INTENT(IN)  :: V_rs(:,:)         !>!< total potential in real space
      REAL(DP),            INTENT(IN)  :: ecut_srb !>!< energy cutoff to include Gamma^8
      INTEGER,             INTENT(IN)  :: nspin
      TYPE(ham_expansion), INTENT(inout) :: ham              !>!< Hamiltonian
      logical, intent(in), optional :: saved
    end subroutine build_h_coeff_cuda

    subroutine copy_pseudo_cuda(pp)
      USE srb_types, only : pseudop
      type(pseudop), intent(inout) :: pp
    end subroutine copy_pseudo_cuda

    SUBROUTINE build_h_matrix(ham, qpoint, pp, spin, ham_matrix)
      USE kinds,   ONLY : DP
      USE srb_types, ONLY : basis, ham_expansion, pseudop
      TYPE(ham_expansion),      INTENT(in)  :: ham
      REAL(DP),                 INTENT(in)  :: qpoint(3)
      type(pseudop), intent(inout) :: pp
      integer, intent(in)     :: spin
      COMPLEX(DP), INTENT(OUT) :: ham_matrix(:,:)
    end subroutine build_h_matrix

    subroutine build_projs(opt_basis, xq, nq, pp)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis, pseudop
      TYPE(basis), intent(in)  :: opt_basis
      REAL(DP),    intent(in)  :: xq(:,:)
      integer,     intent(in)  :: nq
      type(pseudop), intent(inout) :: pp
    end subroutine build_projs

    subroutine build_projs_reduced(opt_basis, xq, nq, pp)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis, pseudop
      TYPE(basis), intent(in)  :: opt_basis
      REAL(DP),    intent(in)  :: xq(:,:)
      integer,     intent(in)  :: nq
      type(pseudop), intent(inout) :: pp
    end subroutine build_projs_reduced

    subroutine build_projs_cuda(opt_basis, qpoint, q, projs, S_matrix, saved_in)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis, pseudop
      TYPE(basis), intent(in)  :: opt_basis
      REAL(DP),    intent(in)  :: qpoint(3)
      integer,     intent(in)  :: q
      TYPE(pseudop), intent(inout) :: projs
      complex(DP), intent(out) :: S_matrix(:,:)
      logical, optional, intent(in) :: saved_in
    end subroutine build_projs_cuda

    subroutine load_projs(q, pp)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : pseudop
      integer,     intent(in)  :: q
      type(pseudop), intent(inout) :: pp
    end subroutine load_projs

    subroutine build_s_matrix(pp, q, Hk)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : pseudop, kproblem
      type(pseudop), intent(inout) :: pp
      integer, intent(in) :: q
      type(kproblem), intent(inout) :: Hk
    end subroutine build_s_matrix

    SUBROUTINE diagonalize (matrix, evals, evecs, num_opt, smatrix, meth_opt, P, Pinv, btype_opt)
      USE kinds, ONLY: DP
      IMPLICIT NONE
      COMPLEX(DP), intent(INOUT) :: matrix(:,:)
      COMPLEX(DP), intent(INOUT) :: evecs(:,:)
      REAL(DP),    intent(OUT)   :: evals(:)
      COMPLEX(DP), intent(INOUT), optional :: smatrix(:,:)
      integer,     intent(in), optional :: num_opt
      integer,     intent(in), optional :: meth_opt
      COMPLEX(DP), intent(IN), optional :: P(:,:)
      COMPLEX(DP), intent(IN), optional :: Pinv(:,:)
      integer, intent(IN), optional :: btype_opt(:)
    end subroutine diagonalize

    SUBROUTINE solve_system_cuda(ham, qpoint, projs, k, spin, states, betawfc, energies)
      USE kinds,   ONLY : DP
      USE srb_types, ONLY : basis, ham_expansion, pseudop, nk_list
      TYPE(ham_expansion),      INTENT(in)  :: ham
      REAL(DP),                 INTENT(in)  :: qpoint(3)
      TYPE(pseudop), intent(in) :: projs
      integer, intent(in)     :: k
      integer, intent(in)     :: spin
      type(nk_list), intent(inout)   :: states
      type(nk_list), intent(inout)   :: betawfc
      REAL(DP), INTENT(OUT) :: energies(:)
    end subroutine solve_system_cuda

    subroutine store_states(wfc, projs, k, states, betawfc)
      use kinds, only : DP
      use srb_types, only : pseudop, nk_list
      complex(DP), intent(in)    :: wfc(:,:)
      type(pseudop), intent(in)    :: projs
      integer, intent(in) :: k
      type(nk_list), intent(inout)   :: states
      type(nk_list), intent(inout)   :: betawfc
    end subroutine store_states

    subroutine build_rho(states, betawfc, wg, wq, opt_basis, nspin, rho, becsum)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis, nk_list
      USE scf, ONLY : scf_type
      type(nk_list), intent(in) :: states
      type(nk_list), intent(in) :: betawfc
      REAL(DP), intent(in) :: wg(:,:)
      REAL(DP), intent(in) :: wq(:)
      type(basis), intent(in) :: opt_basis
      integer, intent(in)   :: nspin
      type(scf_type), intent(inout) :: rho
      real(DP), intent(inout) :: becsum(:,:,:)
    end subroutine build_rho

    subroutine build_rho_reduced(states, betawfc, wg, wq, nspin, rho, becsum)
      USE kinds, ONLY : DP
      use srb_types, only : nk_list
      type(nk_list), intent(in) :: states
      type(nk_list), intent(in) :: betawfc
      REAL(DP), intent(in) :: wg(:,:)
      REAL(DP), intent(in) :: wq(:)
      integer, intent(in)   :: nspin
      COMPLEX(DP), intent(inout) :: rho(:,:,:)
      real(DP), intent(inout) :: becsum(:,:,:)
    end subroutine build_rho_reduced

    subroutine transform_rho(rho_compact, opt_basis, rho)
      USE kinds, ONLY : DP
      USE srb_types, ONLY : basis
      USE scf, ONLY : scf_type
      COMPLEX(DP), intent(inout) :: rho_compact(:,:,:)
      type(basis), intent(in) :: opt_basis
      type(scf_type), intent(inout) :: rho
    end subroutine transform_rho

    subroutine backload(qpoints, srb, states)
      USE srb_types, ONLY : basis, nk_list
      USE kpoint, only : k_list
      type(k_list), intent(in) :: qpoints
      type(basis), intent(in) :: srb
      type(nk_list), intent(inout) :: states
    end subroutine backload

  end interface


CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE init_srb()
  !-----------------------------------------------------------------------
    !_
    USE mp_global, ONLY: intra_bgrp_comm
    USE mp, ONLY: mp_size, mp_rank, mp_get_comm_null
    use control_flags, only : tr2

    USE input_parameters, ONLY : wq
    USE input_parameters, ONLY : rho_reduced_in=>rho_reduced
    use input_parameters, ONLY :  basis_life_in=>basis_life
    use input_parameters, ONLY :  freeze_basis_in=>freeze_basis
    use input_parameters, ONLY :  srb_debug_in=>srb_debug
    use input_parameters, ONLY :  use_cuda_in=>use_cuda

    use scalapack_mod, only : scalapack_init
    USE wvfct, ONLY: nbnd
    USE uspp, ONLY: nkb
    USE cell_base, ONLY : at, bg, tpiba
    USE klist, ONLY: nkstot
    USE lsda_mod, only : nspin
    !
    IMPLICIT NONE
    integer i

    ! module transfer
    basis_life = basis_life_in
    if (freeze_basis_in < 0) then
      freeze_basis = sqrt(tr2)
    else
      freeze_basis = freeze_basis_in
    endif
    rho_reduced = rho_reduced_in
    srb_debug = srb_debug_in
    use_cuda = use_cuda_in

    ! default to all bands
    allocate(band_index(nbnd))
    do i = 1, nbnd
     band_index(i) = i
    enddo

    ! transformation matrices
    trkin2nlp = transpose(at)/tpiba
    trnlp2kin = bg * tpiba

    wq = wq / sum(wq)
    if (nspin == 1) wq = wq * 2.d0

    allocate(ets(nbnd, qpoints%nred*nspin))

    !==========================================================================
    ! Init splines
    !==========================================================================
    !call init_splines()

    call scalapack_init()


  END SUBROUTINE init_srb
END MODULE srb


