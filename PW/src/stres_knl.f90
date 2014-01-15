!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine stres_knl (sigmanlc, sigmakin)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi, e2
  USE cell_base,            ONLY: omega, alat, at, bg, tpiba, tpiba2
  USE gvect,                ONLY: g, ngm
  USE klist,                ONLY: nks, xk, ngk
  USE uspp,                 ONLY: nkb, vkb, okvan
  USE becmod,               ONLY: becp, calbec, allocate_bec_type, deallocate_bec_type
  USE io_files,             ONLY: iunwfc, nwordwfc, iunigk
  USE buffers,              ONLY: get_buffer
  USE symme,                ONLY: symmatrix
  USE wvfct,                ONLY: npw, npwx, igk, wg, qcutz, ecfixed, q2sigma
  USE control_flags,        ONLY: gamma_only
  USE noncollin_module,     ONLY: noncolin, npol
  USE wavefunctions_module, ONLY: evc
  USE wvfct,                ONLY: ecutwfc
  use lsda_mod, only : nspin
  USE mp_pools,             ONLY: inter_pool_comm
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp_global,            ONLY: intra_pool_comm
  use mp_global,            only: me_pot, npot, my_pot_id, inter_pot_comm
  USE mp,                   ONLY: mp_sum
  use input_parameters, only : use_srb
  use srb, only : qpoints, states, bstates, wgq, scb
  use srb_matrix, only : dmat, copy_dmat, setup_dmat, serial_scope

  implicit none
  real(DP) :: sigmanlc (3, 3), sigmakin (3, 3)
  real(DP), allocatable :: gk (:,:), kfac (:)
  real(DP) :: twobysqrtpi, gk2, arg
  real(DP), parameter :: k_gamma(3) = 0.d0
  integer, allocatable :: igk2(:)
  real(DP), allocatable :: g2kin(:)
  real(DP) :: xk_tmp(3)
  integer :: nbasis
  integer :: ik, l, m, i, ibnd, ibnd_g, is, s, q, nbnd, nbnd_l
  complex(DP), parameter :: zero = (0.d0, 0.d0), one = (1.d0, 0.d0)
  type(dmat) :: tmp_mat, serial_mat


  allocate (gk(  3, npwx))    
  allocate (kfac(   npwx))    

  sigmanlc(:,:) =0.d0
  sigmakin(:,:) =0.d0
  twobysqrtpi = 2.d0 / sqrt (pi)

  kfac(:) = 1.d0

  if (use_srb)  then
    nbnd = states%host_ar(1)%desc(4)
    nbnd_l = size(states%host_ar(1)%dat,2)
    allocate(igk2(ngm), g2kin(ngm))
    call gk_sort(k_gamma, ngm, g, ecutwfc/tpiba2, npw, igk2, g2kin)
    
    igk = igk2

    npw = size(scb%elements, 1)
    nbasis = size(scb%elements, 2) 
    deallocate(gk, kfac)
    allocate(gk(3,npw), kfac(npw))
    kfac(:) = 1.d0
  do ik = 1, qpoints%nred
     ! setup some q-point related stuff
     xk_tmp = matmul( bg, qpoints%xr(:,ik) - floor(qpoints%xr(:,ik)) )
     do i = 1, npw
        gk (:, i) = (xk_tmp (:) + g (:, igk2 (i) ) ) * tpiba
        if (qcutz.gt.0.d0) then 
           gk2 = gk (1, i) **2 + gk (2, i) **2 + gk (3, i) **2
           arg = ( (gk2 - ecfixed) / q2sigma) **2
           kfac (i) = 1.d0 + qcutz / q2sigma * twobysqrtpi * exp ( - arg)
        endif
     enddo
    do s = 1, nspin
     call copy_dmat(tmp_mat, states%host_ar(1))
     ! pull the wavefunctions
     q = (ik+(s-1)*(qpoints%nred+npot)-1)/npot + 1
     if (size(states%host_ar) == 1) then
       if (MOD(ik-1, npot) == my_pot_id) then
         call get_buffer(tmp_mat%dat, size(tmp_mat%dat), states%file_unit,q)
       endif
     else
       if (MOD(ik-1, npot) == my_pot_id) then
         tmp_mat%dat = states%host_ar(q)%dat
       endif
     endif
     call setup_dmat(serial_mat, nbasis, tmp_mat%desc(4), scope_in = serial_scope)
     if (MOD(ik-1,npot) == my_pot_id) then
       if (me_pot /= 0) serial_mat%desc(2) = -1
       call pzgemr2d(nbasis, tmp_mat%desc(4), &
                     tmp_mat%dat, 1, 1, tmp_mat%desc, &
                     serial_mat%dat, 1, 1, serial_mat%desc, &
                     tmp_mat%desc(2))
     endif
     call mp_sum(serial_mat%dat, intra_pool_comm)

     ! transform them into the PW basis
     call ZGEMM('N', 'N', npw, tmp_mat%desc(4), nbasis, one, &
                scb%elements, npw, &
                serial_mat%dat, nbasis, zero, &
                evc, npwx)
     !
     !   kinetic contribution
     !
     do l = 1, 3
        do m = 1, l
           do ibnd = 1, nbnd
              do i = 1, npw
                 if (noncolin) then
                    sigmakin (l, m) = sigmakin (l, m) + wgq (ibnd, ik+(s-1)*qpoints%nred) * &
                     gk (l, i) * gk (m, i) * kfac (i) * &
                     ( DBLE (CONJG(evc(i     ,ibnd))*evc(i     ,ibnd)) + &
                       DBLE (CONJG(evc(i+npwx,ibnd))*evc(i+npwx,ibnd)))
                 else  
                    sigmakin (l, m) = sigmakin (l, m) + wgq(ibnd, ik+(s-1)*qpoints%nred) * &
                        gk (l, i) * gk (m, i) * kfac (i) * &
                          DBLE (CONJG(evc (i, ibnd) ) * evc (i, ibnd) )
                 end if   
              enddo
           enddo
        enddo
     enddo
     !
     !  contribution from the  nonlocal part
     !
     CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm )
     if (okvan) then
     call copy_dmat(tmp_mat, bstates%host_ar(1))
     if (size(bstates%host_ar) == 1) then
       if (MOD(ik-1, npot) == my_pot_id) then
         call get_buffer(tmp_mat%dat, nkb*nbnd_l, bstates%file_unit, q)
       else
         tmp_mat%dat = cmplx(0.d0, kind=DP)
       endif
     else
       if (MOD(ik-1, npot) == my_pot_id) then
         tmp_mat%dat = bstates%host_ar(q)%dat ! memory leak
       else
         tmp_mat%dat = cmplx(0.d0, kind=DP)
       endif
     endif
     !call mp_sum(tmp_mat%dat, inter_pot_comm)
     call setup_dmat(serial_mat, nkb, tmp_mat%desc(4), scope_in = serial_scope)
     if (MOD(ik-1,npot) == my_pot_id) then
       if (me_pot /= 0) serial_mat%desc(2) = -1
       call pzgemr2d(nkb, tmp_mat%desc(4), &
                     tmp_mat%dat, 1, 1, tmp_mat%desc, &
                     serial_mat%dat, 1, 1, serial_mat%desc, &
                     tmp_mat%desc(2))
     endif
     call mp_sum(serial_mat%dat, intra_pool_comm)
     becp%k = serial_mat%dat
     endif 
     call stress_us_srb (ik+(s-1)*qpoints%nred, gk, sigmanlc)
     CALL deallocate_bec_type ( becp )
    enddo !spin
  enddo
  !========
  else
  !========

  if (nks.gt.1) rewind (iunigk)
  do ik = 1, nks
     npw = ngk(ik)
     if (nks > 1) then
        read (iunigk) igk
        call get_buffer (evc, nwordwfc, iunwfc, ik)
     endif
     do i = 1, npw
        gk (1, i) = (xk (1, ik) + g (1, igk (i) ) ) * tpiba
        gk (2, i) = (xk (2, ik) + g (2, igk (i) ) ) * tpiba
        gk (3, i) = (xk (3, ik) + g (3, igk (i) ) ) * tpiba
        if (qcutz.gt.0.d0) then
           gk2 = gk (1, i) **2 + gk (2, i) **2 + gk (3, i) **2
           arg = ( (gk2 - ecfixed) / q2sigma) **2
           kfac (i) = 1.d0 + qcutz / q2sigma * twobysqrtpi * exp ( - arg)
        endif
     enddo
     !
     !   kinetic contribution
     !
     do l = 1, 3
        do m = 1, l
           do ibnd = 1, nbnd
              do i = 1, npw
                 if (noncolin) then
                    sigmakin (l, m) = sigmakin (l, m) + wg (ibnd, ik) * &
                     gk (l, i) * gk (m, i) * kfac (i) * &
                     ( DBLE (CONJG(evc(i     ,ibnd))*evc(i     ,ibnd)) + &
                       DBLE (CONJG(evc(i+npwx,ibnd))*evc(i+npwx,ibnd)))
                 else
                    sigmakin (l, m) = sigmakin (l, m) + wg (ibnd, ik) * &
                        gk (l, i) * gk (m, i) * kfac (i) * &                       
                          DBLE (CONJG(evc (i, ibnd) ) * evc (i, ibnd) )
                 end if
              enddo
           enddo
        enddo

     enddo
     !
     !  contribution from the  nonlocal part
     !
!     CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm ) 
!     if (.not. gamma_only) then
!       IF ( nks > 1 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
!       CALL calbec( npw, vkb, evc, becp )
!     endif
     call stres_us (ik, gk, sigmanlc)
!     CALL deallocate_bec_type ( becp ) 
  enddo
  endif 
  !
  ! add the US term from augmentation charge derivatives
  !
  call addusstres (sigmanlc)
  !
  call mp_sum( sigmakin, intra_bgrp_comm )
  call mp_sum( sigmanlc, intra_bgrp_comm )
  call mp_sum( sigmakin, inter_pool_comm )
  call mp_sum( sigmanlc, inter_pool_comm )
  !
  do l = 1, 3
     do m = 1, l - 1
        sigmanlc (m, l) = sigmanlc (l, m)
        sigmakin (m, l) = sigmakin (l, m)
     enddo
  enddo
  !
  write(*,*) sigmakin
  if (gamma_only) then
     sigmakin(:,:) = 2.d0 * e2 / omega * sigmakin(:,:)
  else
     sigmakin(:,:) = e2 / omega * sigmakin(:,:)
  end if
  sigmanlc(:,:) = -1.d0 / omega * sigmanlc(:,:)
  !
  ! symmetrize stress
  !
  call symmatrix ( sigmakin )
  call symmatrix ( sigmanlc )

  deallocate(kfac)
  deallocate(gk)
  return
end subroutine stres_knl

