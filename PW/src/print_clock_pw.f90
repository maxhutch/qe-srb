!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock_pw()
   !---------------------------------------------------------------------------
   !
   ! ... this routine prints out the clocks at the end of the run
   ! ... it tries to construct the calling tree of the program.
   !
   USE io_global,          ONLY : stdout
   USE control_flags,      ONLY : isolve, iverbosity, gamma_only
   USE paw_variables,      ONLY : okpaw
   USE uspp,               ONLY : okvan
   USE realus,             ONLY : real_space
   USE ldaU,               ONLY : lda_plus_U
   USE funct,              ONLY : dft_is_hybrid
   use input_parameters,   ONLY : use_srb
#ifdef __ENVIRON
   USE environ_base,       ONLY : do_environ
   USE environ_info,       ONLY : environ_clock
#endif
   !
   IMPLICIT NONE
   !
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'init_run' )
   CALL print_clock( 'electrons' )
   CALL print_clock( 'update_pot' )
   CALL print_clock( 'forces' )
   CALL print_clock( 'stress' )
   !
   WRITE( stdout, '(/5x,"Called by init_run:")' )
   CALL print_clock( 'wfcinit' )
   CALL print_clock( 'potinit' )
   CALL print_clock( 'realus' )
   IF ( iverbosity > 0 ) THEN
      CALL print_clock( 'realus:boxes' )
      CALL print_clock( 'realus:spher' )
      CALL print_clock( 'realus:qsave' )
   END IF
   !
   WRITE( stdout, '(/5x,"Called by electrons:")' )
   CALL print_clock( 'c_bands' )
   CALL print_clock( 'sum_band' )
   CALL print_clock( 'v_of_rho' )
   IF ( iverbosity > 0 ) THEN
      CALL print_clock( 'v_h' )
      CALL print_clock( 'v_xc' )
      CALL print_clock( 'v_xc_meta' )
   END IF
   CALL print_clock( 'newd' )
   IF ( iverbosity > 0 ) THEN
      CALL print_clock( 'newd:fftvg' )
      CALL print_clock( 'newd:qvan2' )
      CALL print_clock( 'newd:int1' )
      CALL print_clock( 'newd:int2' )
   END IF
   CALL print_clock( 'mix_rho' )

   CALL print_clock( 'vdW_energy' )
   CALL print_clock( 'vdW_ffts' )
   CALL print_clock( 'vdW_v' )
   
   !
   WRITE( stdout, '(/5x,"Called by c_bands:")' )
   CALL print_clock( 'init_us_2' )
   IF ( isolve == 0 ) THEN
      IF ( gamma_only ) THEN
         CALL print_clock( 'regterg' )
      ELSE
         CALL print_clock( 'cegterg' )
      ENDIF
   ELSE 
      IF ( gamma_only ) THEN
         CALL print_clock( 'rcgdiagg' )
      ELSE
         CALL print_clock( 'ccgdiagg' )
      ENDIF
      CALL print_clock( 'wfcrot' )
   ENDIF
   !
   IF ( iverbosity > 0)  THEN
      WRITE( stdout, '(/5x,"Called by sum_band:")' )
      CALL print_clock( 'sum_band:becsum' )
      CALL print_clock( 'addusdens' )
      CALL print_clock( 'addus:qvan2' )
      CALL print_clock( 'addus:strf' )
      CALL print_clock( 'addus:aux2' )
      CALL print_clock( 'addus:aux' )
      CALL print_clock( 'rho_finish' )
   ENDIF
   !
   IF ( isolve == 0 ) THEN
      WRITE( stdout, '(/5x,"Called by *egterg:")' )
   ELSE 
      WRITE( stdout, '(/5x,"Called by *cgdiagg:")' )
   END IF
   !
   IF (real_space ) THEN
    WRITE( stdout, '(/5x,"Called by real space routines:")' )
    CALL print_clock ( 'realus' )
    CALL print_clock ( 'betapointlist' )
    CALL print_clock ( 'addusdens' )
    CALL print_clock ( 'calbec_rs' )
    CALL print_clock ( 's_psir' )
    CALL print_clock ( 'add_vuspsir' )
    CALL print_clock ( 'fft_orbital' )
    CALL print_clock ( 'bfft_orbital' )
    CALL print_clock ( 'v_loc_psir' )
   ELSE
    CALL print_clock( 'h_psi' )
    CALL print_clock( 's_psi' )
    CALL print_clock( 'g_psi' )
   ENDIF
   IF ( gamma_only ) THEN
      CALL print_clock( 'rdiaghg' )
      IF ( iverbosity > 0 )  THEN
         CALL print_clock( 'regterg:overlap' )
         CALL print_clock( 'regterg:update' )
         CALL print_clock( 'regterg:last' )
         CALL print_clock( 'rdiaghg:choldc' )
         CALL print_clock( 'rdiaghg:inversion' )
         CALL print_clock( 'rdiaghg:paragemm' )
      ENDIF
   ELSE
      CALL print_clock( 'cdiaghg' )
      IF ( iverbosity > 0 )  THEN
         CALL print_clock( 'cegterg:overlap' )
         CALL print_clock( 'cegterg:update' )
         CALL print_clock( 'cegterg:last' )
         CALL print_clock( 'cdiaghg:choldc' )
         CALL print_clock( 'cdiaghg:inversion' )
         CALL print_clock( 'cdiaghg:paragemm' )
      END IF
   END IF
   !
   WRITE( stdout, '(/5x,"Called by h_psi:")' )
   IF ( iverbosity > 0 )  THEN
      CALL print_clock( 'h_psi:init' )
      CALL print_clock( 'h_psi:vloc' )
      CALL print_clock( 'h_psi:vnl' )
   END IF
   CALL print_clock( 'add_vuspsi' )
   CALL print_clock( 'vhpsi' )
   CALL print_clock( 'h_psi_meta' )
   !
   WRITE( stdout, '(/5x,"Called by force/stress:")' )
   CALL print_clock( 'k_stress' )
   CALL print_clock( 'k_force' )
   !
   WRITE( stdout, '(/5X,"General routines")' )
   !
   CALL print_clock( 'calbec' )
   CALL print_clock( 'fft' )
   CALL print_clock( 'ffts' )
   CALL print_clock( 'fftw' )
   CALL print_clock( 'interpolate' )
   CALL print_clock( 'davcio' )
   !    
   WRITE( stdout, * )
   !
#if defined (__MPI)
   WRITE( stdout, '(5X,"Parallel routines")' )
   !
   CALL print_clock( 'reduce' )
   CALL print_clock( 'fft_scatter' )
   CALL print_clock( 'ALLTOALL' )
#endif
   !
   IF ( lda_plus_U ) THEN
      WRITE( stdout, '(5X,"Hubbard U routines")' )
      CALL print_clock( 'new_ns' )
      CALL print_clock( 'vhpsi' )
      CALL print_clock( 'force_hub' )
      CALL print_clock( 'stres_hub' )
   ENDIF
   !
   IF ( dft_is_hybrid() ) THEN
      WRITE( stdout, '(/,5X,"EXX routines")' )
      CALL print_clock( 'exx_grid' )
      CALL print_clock( 'exxinit' )
      CALL print_clock( 'vexx' )
      !CALL print_clock( 'vexx_ngmloop' )
      CALL print_clock( 'exxenergy' )
      CALL print_clock( 'exxen2' )
      !CALL print_clock( 'exxen2_ngmloop' )
      CALL print_clock ('cycleig')
      IF( okvan) THEN
        WRITE( stdout, '(/,5X,"EXX+US routines")' )
        CALL print_clock( 'becxx' )
        CALL print_clock( 'addusxx' )
        CALL print_clock( 'newdxx' )
        CALL print_clock( 'nlxx_pot' )
      ENDIF
   ENDIF
   !
   IF ( okpaw ) THEN
      WRITE( stdout, '(/,5X,"PAW routines")' )
      ! radial routines:
      CALL print_clock ('PAW_pot')
      CALL print_clock ('PAW_newd')
      CALL print_clock ('PAW_int')
      CALL print_clock ('PAW_ddot')
      CALL print_clock ('PAW_rad_init')
      CALL print_clock ('PAW_energy')
      CALL print_clock ('PAW_symme')
      ! second level routines:
      CALL print_clock ('PAW_rho_lm')
      CALL print_clock ('PAW_h_pot')
      CALL print_clock ('PAW_xc_pot')
      CALL print_clock ('PAW_lm2rad')
      CALL print_clock ('PAW_rad2lm')
      ! third level, or deeper:
      CALL print_clock ('PAW_rad2lm3')
      CALL print_clock ('PAW_gcxc_v')
      CALL print_clock ('PAW_div')
      CALL print_clock ('PAW_grad')
      IF ( dft_is_hybrid() ) THEN
        WRITE( stdout, '(/,5X,"PAW+EXX routines")' )
        CALL print_clock("PAW_newdxx")
        CALL print_clock("PAW_xx_nrg")
        CALL print_clock('PAW_keeq')
      ENDIF
   END IF
   !
   CALL plugin_clock()
   !
   IF ( use_srb ) THEN
      write(stdout, *)
      write(stdout, '(5X,"SHIRLEY routines")' )
      call print_clock( 'srb' )
      call print_clock( ' build_red_basis' )
      call print_clock( '  mirrors' )
      call print_clock( '  SVD' )
      call print_clock( '  make_S' )
      call print_clock( '  diag_S' )
      call print_clock( '  expand' )
      call print_clock( ' build_h_coeff' )
      call print_clock( '  local_potential' )
      call print_clock( '  Kinetic' )
      call print_clock( '   part1' )
      call print_clock( ' build_proj' )
      call print_clock( '  make_vbas')
      call print_clock( '  make_St')
      call print_clock( '  make_vkb')
      call print_clock( '  apply_S')
      call print_clock( '  proj_init')
      call print_clock( '  proj_gemm')
      call print_clock( '  proj_load')
      call print_clock( '  proj_save')
      call print_clock( ' build_mat' )
      call print_clock( ' diagonalize' )
      call print_clock( ' store' )
      call print_clock( ' build_rho' )
      call print_clock( '  rho_matrix')
      call print_clock( '  addproj')
      call print_clock( '  svd')
      call print_clock( '  gemm')
      call print_clock( '  fft')
      call print_clock( '  acc')
      call print_clock( '  transpose')
      call print_clock( ' finish_rho' )
      call print_clock( ' other' )
   ENDIF
#ifdef __ENVIRON
   IF ( do_environ ) call environ_clock( stdout )
#endif
   !
   RETURN
   !
END SUBROUTINE print_clock_pw
