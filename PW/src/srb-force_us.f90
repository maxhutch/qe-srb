!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE force_us_srb( forcenl )
  !----------------------------------------------------------------------------
  !
  ! ... nonlocal potential contribution to forces
  ! ... wrapper
  !
  USE kinds,                ONLY : DP
  use constants, only : rytoev
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE klist,                ONLY : nks, xk, ngk
  USE gvect,                ONLY : g
  USE uspp,                 ONLY : nkb, qq, deeq, qq_so, deeq_nc, okvan
  USE uspp_param,           ONLY : upf, nh, newpseudo, nhm
  USE wvfct,                ONLY : wg, et, ecutwfc
  USE lsda_mod,             ONLY : lsda, current_spin, isk, nspin
  USE symme,                ONLY : symvector
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : npol, noncolin
  USE spin_orb,             ONLY : lspinorb
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE buffers,              ONLY : get_buffer, open_buffer, save_buffer
  USE becmod,               ONLY : bec_type, becp, allocate_bec_type, deallocate_bec_type
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm, intra_pool_comm
  USE mp_global,            ONLY : npot, inter_pot_comm
  USE mp,                   ONLY : mp_sum, mp_get_comm_null
  !
  IMPLICIT NONE
  !
  ! ... the dummy variable
  !
  REAL(DP) :: forcenl(3,nat)
  ! output: the nonlocal contribution
  !
  IF ( gamma_only ) THEN
     !
     CALL force_us_gamma_srb( forcenl )
     !
  ELSE
     !
     CALL force_us_k_srb( forcenl )
     !
  END IF  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_gamma_srb( forcenl )
       !-----------------------------------------------------------------------
       !
       ! ... calculation at gamma
       !
       USE becmod, ONLY : calbec
       use wvfct, only : igk, npw, npwx, nbnd
       use uspp, only : vkb
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       TYPE(bec_type) :: rdbecp (3)
       ! auxiliary variable, contains <dbeta|psi>
       COMPLEX(DP), ALLOCATABLE :: vkb1(:,:)
       ! auxiliary variable contains g*|beta>
       REAL(DP) :: ps
       INTEGER       :: ik, ipol, ibnd, ibnd_loc, ig, ih, jh, na, nt, ikb, jkb, ijkb0
       ! counters
       !
       !
       forcenl(:,:) = 0.D0
       !
       DO ipol = 1, 3
          CALL allocate_bec_type ( nkb, nbnd, rdbecp(ipol), intra_bgrp_comm )   
       END DO
       ALLOCATE( vkb1(  npwx, nkb ) ) 
       !   
       IF ( nks > 1 ) REWIND iunigk
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk (ik)
          IF ( nks > 1 ) THEN
             READ( iunigk ) igk
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
             IF ( nkb > 0 ) &
                CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          CALL calbec ( npw, vkb, evc, becp )
          !
          DO ipol = 1, 3
             DO jkb = 1, nkb
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb) * (0.D0,-1.D0) * g(ipol,igk(ig))
                END DO
!$OMP END PARALLEL DO
             END DO
             !
             CALL calbec ( npw, vkb1, evc, rdbecp(ipol) )
             !
          END DO
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd_loc = 1, becp%nbnd_loc
                         ibnd = ibnd_loc + becp%ibnd_begin - 1
                         ps = deeq(ih,ih,na,current_spin) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         DO ipol = 1, 3
                            forcenl(ipol,na) = forcenl(ipol,na) - &
                                       ps * wg(ibnd,ik) * 2.D0 * tpiba * &
                                       rdbecp(ipol)%r(ikb,ibnd_loc) *becp%r(ikb,ibnd_loc)
                         END DO
                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd_loc = 1, becp%nbnd_loc
                               ibnd = ibnd_loc + becp%ibnd_begin - 1
                               ps = deeq(ih,jh,na,current_spin) - &
                                    et(ibnd,ik) * qq(ih,jh,nt)
                               DO ipol = 1, 3
                                  forcenl(ipol,na) = forcenl(ipol,na) - &
                                     ps * wg(ibnd,ik) * 2.d0 * tpiba * &
                                     (rdbecp(ipol)%r(ikb,ibnd_loc) *becp%r(jkb,ibnd_loc) + &
                                      rdbecp(ipol)%r(jkb,ibnd_loc) *becp%r(ikb,ibnd_loc) )
                               END DO
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          END DO
       END DO
       !
       IF( becp%comm /= mp_get_comm_null() ) CALL mp_sum( forcenl, becp%comm )
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
       CALL addusforce( forcenl )
       !
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces
       !
       CALL symvector ( nat, forcenl )
       !
       DEALLOCATE( vkb1 )
       DO ipol = 1, 3
          CALL deallocate_bec_type ( rdbecp(ipol) )   
       END DO
       !
       RETURN
       !
     END SUBROUTINE force_us_gamma_srb
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE force_us_k_srb( forcenl )
       !-----------------------------------------------------------------------
       !  
       USE becmod, ONLY : calbec
       use srb, only : qpoints, states, bstates, wgq, scb, ets, spp
       use srb_matrix, only : dmat, copy_dmat
       use gvect, only : ngm
       use cell_base, only : tpiba2, tpiba
       use wvfct, only : npwx_int => npwx
       use mp_global, only :  npot, my_pot_id, intra_pot_comm, inter_pot_comm
       IMPLICIT NONE
       !
       REAL(DP) :: forcenl(3,nat)
       REAL(DP) :: qcart(3)
       real(DP), parameter :: k_gamma(3) = 0.d0
       integer, allocatable :: igk(:)
       real(DP), allocatable :: g2kin(:)
       COMPLEX(DP), ALLOCATABLE :: dbecp(:,:,:), dbecp_nc(:,:,:,:)
       ! auxiliary variable contains <beta|psi> and <dbeta|psi>
       COMPLEX(DP), ALLOCATABLE :: vkb(:,:), vkb1(:,:)
       COMPLEX(DP), allocatable :: dprojs(:,:,:)
       ! auxiliary variable contains g*|beta>
       COMPLEX(DP) :: psc(2,2), fac
       COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
       REAL(DP), ALLOCATABLE :: deff(:,:,:)
       REAL(DP) :: ps
       INTEGER       :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0, &
                        is, js, ijs, s
       integer :: ibnd_g
       COMPLEX(DP), parameter :: zero = cmplx(0.d0, kind=DP), one = cmplx(1.d0, kind=DP)
       integer :: npw, npwx_tmp, nbnd
       logical :: stat
       integer, external :: indxl2g
       type(dmat) :: tmp_mat
       ! counters
       !
       !
       forcenl(:,:) = 0.D0
       if (nkb < 1) return
       !
       nbnd = size(states%host_ar(1)%dat,2)
       !
       IF (noncolin) then
          ALLOCATE( dbecp_nc(nkb,npol,nbnd,3) )
          ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
       ELSE
          ALLOCATE( dbecp( nkb, nbnd, 3 ) )    
          ALLOCATE( deff(nhm,nhm,nat) )
       ENDIF
       npw = size(scb%elements, 1)
       ALLOCATE(vkb( npw, nkb ) )
       ALLOCATE(vkb1( npw, nkb ) )
       allocate(dprojs(scb%length, nkb, 0:3))
       allocate(becp%k(nkb, nbnd))
       !
       ! ... the forces are a sum over the K points and the bands
       !
       allocate(igk(ngm), g2kin(ngm))
       call gk_sort(k_gamma(:), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
       DO ik = 1, qpoints%nred
          ! transform the d beta/ dR into the SCB
          ! get <G|beta(q)>
          qcart = matmul( bg, qpoints%xr(:,ik) - floor(qpoints%xr(:,ik)))
          npwx_tmp = npwx_int; npwx_int = npw
          CALL init_us_2_srb(npw, igk, qcart, -1, vkb )
          npwx_int = npwx_tmp

          ! make <G| d beta(q)/dR> = <G | g | beta(q)>
          DO ipol = 0, 3
             if (ipol == 0) then 
               if (size(bstates%host_ar) > 0) cycle 
               vkb1 = vkb
             else
             DO jkb = 1, nkb
                DO ig = 1, npw
                   vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*(g(ipol,igk(ig)) + qcart(ipol))
                END DO
             END DO
             endif
             ! transform into the SCB
             call ZGEMM('C', 'N', scb%length, nkb, npw, one, &
                        scb%elements, npw, & 
                        vkb1, npw, zero, &
                        dprojs(:,:,ipol), scb%length)
          END DO
          call mp_sum(dprojs, intra_pool_comm)
        do s = 1, nspin
          ! grab the states <B|psi>
          call copy_dmat(tmp_mat, states%host_ar(1))
          if (size(states%host_ar) == 1) then
            if (MOD(ik-1, npot) == my_pot_id) then
              call get_buffer(tmp_mat%dat, scb%length*nbnd, states%file_unit, (ik+(s-1)*(qpoints%nred+npot)-1)/npot + 1)
            endif
          else
            if (MOD(ik-1, npot) == my_pot_id) then
              tmp_mat%dat = states%host_ar((ik+(s-1)*(qpoints%nred+npot)-1)/npot + 1)%dat
            endif
          endif
          call mp_sum(tmp_mat%dat, inter_pot_comm)

          ! grab the projector states <beta|psi>
          becp%k = cmplx(0.d0, kind=DP)
          if (size(bstates%host_ar) == 0) then
            call ZGEMM('C', 'N', nkb, nbnd, scb%length, one, &
                       dprojs(:,:,0), scb%length, &
                       tmp_mat%dat, scb%length, zero, &
                       becp%k, nkb)
            if (bstates%file_unit < 0) then
              bstates%file_unit = 4827
              call open_buffer(bstates%file_unit, 'bstates', size(becp%k), 1, stat)
            endif
            if (MOD(ik-1, npot) == my_pot_id) then
              call save_buffer(becp%k, size(becp%k), bstates%file_unit,(ik+(s-1)*(qpoints%nred+npot)-1)/npot + 1)
            else
              becp%k = cmplx(0.d0, kind=DP)
            endif
          else if (size(bstates%host_ar) == 1) then
            if (MOD(ik-1, npot) == my_pot_id) then
              call get_buffer(becp%k, nkb*nbnd, bstates%file_unit, (ik+(s-1)*(qpoints%nred+npot)-1)/npot + 1)
            endif
          else
            if (MOD(ik-1, npot) == my_pot_id) then
              becp%k = bstates%host_ar((ik+(s-1)*(qpoints%nred+npot)-1)/npot + 1)%dat ! memory leak
            endif
          endif
          call mp_sum(becp%k, inter_pot_comm)
          
          ! make the product d <beta|psi>/dR = <(d beta)/dR|psi>
          do ipol = 1, 3
          call ZGEMM('C', 'N', nkb, nbnd, scb%length, one, &
                     dprojs(:,:,ipol), scb%length, &
                     tmp_mat%dat, scb%length, zero, &
                     dbecp(:,:,ipol), nkb)
          enddo
          DO ibnd = 1, nbnd 
             ibnd_g = indxl2g(ibnd, &
                              states%host_ar(1)%desc(6), states%host_ar(1)%mycol, &
                              states%host_ar(1)%desc(8), states%host_ar(1)%npcol)
             ! scale the D and S matrices by the energy
             IF (noncolin) THEN
                CALL compute_deff_nc(deff_nc,et(ibnd_g,ik))
             ELSE
                current_spin = s
                CALL compute_deff(deff,ets(ibnd_g,ik+(s-1)*qpoints%nred) )
             ENDIF
             ! grab the right weight
             fac=wgq(ibnd_g,ik+(s-1)*qpoints%nred)*tpiba
             ! from here, it is a literal copy of force_us.f90
             ijkb0 = 0
             DO nt = 1, ntyp
                DO na = 1, nat
                   IF ( ityp(na) == nt ) THEN
                      DO ih = 1, nh(nt)
                         ikb = ijkb0 + ih
                         IF (noncolin) THEN
                            DO ipol=1,3
                               ijs=0
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     forcenl(ipol,na) = forcenl(ipol,na)- &
                                         deff_nc(ih,ih,na,ijs)*fac*( &
                                         CONJG(dbecp_nc(ikb,is,ibnd,ipol))* &
                                         becp%nc(ikb,js,ibnd)+ &
                                         CONJG(becp%nc(ikb,is,ibnd))* &
                                         dbecp_nc(ikb,js,ibnd,ipol) )
                                  END DO
                               END DO
                            END DO
                         ELSE
                            DO ipol=1,3
                               forcenl(ipol,na) = forcenl(ipol,na) - &
                                  2.D0 * fac * deff(ih,ih,na)*&
                                      DBLE( CONJG( dbecp(ikb,ibnd,ipol) ) * &
                                            becp%k(ikb,ibnd) )
                            END DO
                         END IF
                         !
                         IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                            DO jh = ( ih + 1 ), nh(nt)
                               jkb = ijkb0 + jh
                               IF (noncolin) THEN
                                  DO ipol=1,3
                                     ijs=0
                                     DO is=1,npol
                                        DO js=1,npol
                                           ijs=ijs+1
                                           forcenl(ipol,na)=forcenl(ipol,na)- &
                                           deff_nc(ih,jh,na,ijs)*fac*( &
                                          CONJG(dbecp_nc(ikb,is,ibnd,ipol))* &
                                                 becp%nc(jkb,js,ibnd)+ &
                                          CONJG(becp%nc(ikb,is,ibnd))* &
                                                dbecp_nc(jkb,js,ibnd,ipol))- &
                                           deff_nc(jh,ih,na,ijs)*fac*( &
                                          CONJG(dbecp_nc(jkb,is,ibnd,ipol))* &
                                                becp%nc(ikb,js,ibnd)+ &
                                          CONJG(becp%nc(jkb,is,ibnd))* &
                                                dbecp_nc(ikb,js,ibnd,ipol) )
                                        END DO
                                     END DO
                                  END DO
                               ELSE
                                  DO ipol = 1, 3
                                     forcenl(ipol,na) = forcenl (ipol,na) - &
                                          2.D0 * fac * deff(ih,jh,na)* &
                                       DBLE( CONJG( dbecp(ikb,ibnd,ipol) ) * &
                                             becp%k(jkb,ibnd) +       &
                                             dbecp(jkb,ibnd,ipol) * &
                                             CONJG( becp%k(ikb,ibnd) ) )
                                  END DO
                               END IF
                            END DO !jh
                         END IF ! tvanp
                      END DO ! ih = 1, nh(nt)
                      ijkb0 = ijkb0 + nh(nt)
                   END IF ! ityp(na) == nt
                END DO ! nat
             END DO ! ntyp
          END DO ! nbnd
        end do !spin
       END DO ! nks
       !
       CALL mp_sum(  forcenl , intra_pot_comm )
       !
       DEALLOCATE( vkb )
       DEALLOCATE( vkb1 )
       deallocate( dprojs )
       deallocate( becp%k )
       deallocate(igk, g2kin)
       IF (noncolin) THEN
          DEALLOCATE( dbecp_nc )
          DEALLOCATE( deff_nc )
       ELSE
          DEALLOCATE( dbecp )
          DEALLOCATE( deff )
       ENDIF
       !
       ! ... The total D matrix depends on the ionic position via the
       ! ... augmentation part \int V_eff Q dr, the term deriving from the 
       ! ... derivative of Q is added in the routine addusforce
       !
       CALL addusforce( forcenl ) 
       !
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( forcenl, inter_pool_comm )
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces.
       !
       CALL symvector ( nat, forcenl )
       !
       RETURN
       !
     END SUBROUTINE force_us_k_srb
     !     
END SUBROUTINE force_us_srb
