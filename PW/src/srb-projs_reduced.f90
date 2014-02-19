subroutine build_projs_reduced(opt_basis, xq, nq, pp)
  USE kinds, ONLY : DP
  use constants, only : tpi
  USE srb_types, ONLY : basis, pseudop
  use srb_matrix, only : setup_dmat, load_from_local, pot_scope, print_dmat, g2l

  use input_parameters, only : aux_tol, min_aux_size
  USE us, ONLY : dq, tab
  USE uspp, only : nkb, nhtol, nhtolm, indv
  use uspp_param, only : nh, nhm, lmaxkb, upf
  USE gvect, ONLY : ngm, g, eigts1, eigts2, eigts3, mill
  USE cell_base, ONLY : tpiba, tpiba2, bg
  use ions_base, only : nat, ntyp => nsp, ityp, tau
  use wvfct, only : npwx_int => npwx, ecutwfc
  use mp, only : mp_sum, mp_max
  use mp_global, only : intra_pool_comm, me_pool, nproc_pool
  use mp_global, only : intra_pot_comm, me_pot, nproc_pot, npot, my_pot_id
  use buffers, only : open_buffer, save_buffer, get_buffer, close_buffer
  use srb, only : me_tub, nproc_tub, my_tub_id

  IMPLICIT NONE

  TYPE(basis), intent(in) :: opt_basis
  REAL(DP),    intent(in) :: xq(:,:)
  integer,     intent(in) :: nq
  type(pseudop), intent(inout) :: pp

  ! locals
  integer q, q2, left, nbasis, npw, npwx, npwx_tmp, nqx, stat, vbs, fq
  integer :: a, b, t, ig, ih, i, j, i_l, j_l, jkb, lm
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  real(DP) :: k_gamma(3)
  complex(DP), allocatable :: S(:,:,:), Stmp(:,:), sk(:)
  complex(DP), allocatable :: zbuffer(:,:), vkb(:,:), vkb1z(:,:), vkb3(:,:)
  real(DP), allocatable :: buffer(:,:), vbas(:,:), C(:,:), vkb2(:,:)
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:), sv(:)
  real(DP) :: px, ux, vx, wx, arg
  integer :: i0, i1, i2, i3
  complex(DP) :: phase, prefac
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP)
  complex(DP), parameter :: one = cmplx(1.d0, kind=DP)
  real(DP) :: qcart(3), qtmp(3), inv_norm
  logical ::  info, local
  real(DP), allocatable :: work(:)
  integer, allocatable :: iwork(:)
  integer :: lwork, liwork
  integer :: counter, prow, pcol
  character(len=10) :: fname

  integer,save :: old_basis_length

  COMPLEX(DP), external :: ZDOTC
  REAL(DP), external :: DNRM2

  interface
subroutine build_centered_projs(qcart, t, gk, qg, vq, ylm, igk, vkb0)
  use kinds, only : DP
  real(DP), intent(in) :: qcart(3)
  integer, intent(in) :: t
  real(DP), intent(out) :: gk(:,:), qg(:), vq(:), ylm(:,:)
  integer, intent(in) :: igk(:)
  real(DP), intent(out) :: vkb0(:,:)
end subroutine
  end interface

  ! initial setup
  nbasis = opt_basis%length
  npw = size(opt_basis%elements, 1)
  npwx = npw
  call mp_max(npwx, intra_pool_comm)
  nqx = size(tab,1)
  k_gamma = 0.d0
  allocate(igk(ngm), g2kin(ngm))
  call gk_sort(k_gamma(:), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)

  if (.not. allocated(pp%projs)) then
    ! check ordering
    do a = 2, nat
      if (ityp(a) < ityp(a-1)) then
        CALL errore('srb-projs_reduced','atomic positions out of order', a)
      endif
    enddo

    allocate(pp%projs(ntyp), pp%p_unit(ntyp))
    allocate(pp%na(ntyp), pp%na_off(nat), pp%nt_off(ntyp))
    pp%na = 0; pp%na_off = 0
    do a = 1, nat
      pp%na(ityp(a)) = pp%na(ityp(a)) + 1
    enddo
    pp%nkb_l = 0
    do t = 1, ntyp    
      call setup_dmat(pp%projs(t), nbasis, nh(t)*pp%na(t), nbasis, nh(t), pot_scope)
!      call print_dmat(pp%projs(t))
      pp%nkb_l = pp%nkb_l + size(pp%projs(t)%dat,2) 
      pp%p_unit(t) = 3000+t
      write(fname, '(A6,I4)') "projs_", t
      if (size(pp%projs(t)%dat) > 0) then
        call open_buffer(pp%p_unit(t), trim(fname), size(pp%projs(t)%dat), 1, info)
      endif
    enddo
    jkb = 1
    ! this makes assumptions about atom ordering by type
    pp%nt_off(1) = 1
    do t = 2, ntyp
      pp%nt_off(t) = pp%nt_off(t-1) + pp%na(t-1)
    enddo
    do t = 1, ntyp
      counter = 0
      do a = 1, nat
        if (ityp(a) /= t) cycle
        pp%na_off(a) = jkb
        jkb = jkb + nh(t)
        counter = counter + 1
      enddo
    enddo
    pp%ntyp = ntyp
    pp%nat  = nat
    pp%nkb  = nkb
    old_basis_length = opt_basis%length
  else
    ! resize nbasis
    do t = 1, ntyp    
      call setup_dmat(pp%projs(t), nbasis, nh(t)*pp%na(t), nbasis, nh(t), pot_scope)
    enddo
  endif

  allocate (ylm(npw, (lmaxkb + 1)**2), gk(3, npw), qg(npw), vq(npw))
  ! if we have a longer buffer - adjust size of file dynamically
  do t = 1, ntyp    
    if (opt_basis%length > old_basis_length .and. size(pp%projs(t)%dat) > 0) then
      call close_buffer(pp%p_unit(t),'delete')
      write(fname, '(A6,I4)') "projs_", t
      call open_buffer(pp%p_unit(t), trim(fname), size(pp%projs(t)%dat), 1, info)
      old_basis_length = opt_basis%length
    endif
  enddo

  ! if the auxillary basis isn't in file, make it
  if (pp%b_unit < 0) then 
    pp%b_unit = - pp%b_unit
    pp%v_unit = - pp%v_unit
    call open_buffer(pp%b_unit, 'vbas', npwx*nhm*nq, 1, info)
    call open_buffer(pp%v_unit, 'vkb_x', nhm*nq*nhm, 1, info)

    ! auxillary q-grid
    allocate(vbas(npw, nhm*nq))
    allocate(sv(nhm*nq))
    do t = 1, ntyp ! loop over types
      if (pp%na(t) < min_aux_size) cycle
      call start_clock('  make_vbas')
      allocate(vkb1(npw, nh(t)), vkb1z(npw, nh(t)))
      do q = 1, nq ! loop over aux q-grid
        qtmp = xq(:,q) - floor(xq(:,q))
        qcart = matmul( bg, qtmp )
        call build_centered_projs(qcart, t, gk, qg, vq, ylm, igk, vkb1)
        vbas(:,(q-1)*nh(t)+1:q*nh(t)) = vkb1
      enddo 
      vbs = nh(t) * nq
      ! SVD!
      allocate(C(vbs, vbs))
      call dsyrk('U', 'C', vbs, npw, one, &
                 vbas, npw, zero, &
                 C, vbs)
      call mp_sum(C, intra_pool_comm)
      lwork = 1+ 6*vbs + 2 * vbs * vbs; liwork = 3+5*vbs
      allocate(work(lwork), iwork(liwork))
      sv = 0.
      call dsyevd('V', 'U', vbs, &
              C, vbs, &
              sv, &
              work, lwork, &
              iwork, liwork, stat)
      deallocate(work, iwork)
      if (stat /= 0) write(*,*) "zheevd failed: ", stat
      sv = sqrt(abs(sv))
      fq = vbs
      do ih = 1, fq
        vbs = ih
        if (abs(1 - sum(sv(fq-ih+1:fq))/sum(sv(1:fq))) < aux_tol) exit
      enddo
      allocate(buffer(npw, vbs), zbuffer(npw, vbs))
      call dgemm('N', 'N', npw, vbs, fq, one, &
                 vbas, npw, &
                 C(:,fq-vbs+1:fq), fq, zero, &
                 buffer, npw)
      do ih = 1, vbs
        inv_norm = dnrm2 (npw, buffer(:,ih), 1)
        inv_norm = inv_norm * inv_norm
        call mp_sum(inv_norm, intra_pool_comm)
        inv_norm = 1.d0/sqrt(inv_norm)
        call dscal(npw, inv_norm, buffer(:,ih), 1)
      enddo

      ! end SVD
      call davcio ( buffer, npw*vbs, pp%b_unit, t, +1 )
      deallocate(C)

      if (me_pool == 0) write(*,'(5X,A,F7.2,A,I5,A,I2)') "Aux basis uses ", (1.*vbs)/nh(t), " of ", nq, "; type ", t
      pp%b_size(t) = vbs
      call stop_clock('  make_vbas')


      allocate(vkb2(vbs, nh(t)))
      call start_clock('  make_vkb')
      do q = 1, nq ! loop over aux q-grid
        qtmp = xq(:,q) - floor(xq(:,q))
        qcart = matmul( bg, qtmp )
        call build_centered_projs(qcart, t, gk, qg, vq, ylm, igk, vkb1)

        ! transform vkb1 into the auxillary basis
        vkb1z = cmplx(vkb1, kind=DP)
        call DGEMM('C', 'N', vbs, nh(t), npw, one, &
                 buffer, npw, &
                 vkb1, npw, zero, &
                 vkb2, vbs)
      call mp_sum(vkb2, intra_pool_comm)
      call davcio ( vkb2, vbs*nh(t), pp%v_unit, (t-1)*nq+q, +1 )
      enddo
      call stop_clock('  make_vkb')
      deallocate(buffer, zbuffer)
      deallocate(vkb1, vkb1z, vkb2)
    enddo
    deallocate(vbas, sv)
  endif

  ! transform the projectors
  types: do t = 1, ntyp ! loop over types
    if (pp%na(t) < min_aux_size) then
      allocate(zbuffer(npw, nh(t)*pp%na(t)), vkb3(nbasis, nh(t)*pp%na(t)))
      do q = 1, nq
        call start_clock('  proj_init')
        qtmp = xq(:,q) - floor(xq(:,q))
        qcart = matmul( bg, qtmp )

        npwx_tmp = npwx_int; npwx_int = npw
        call init_us_2_srb(npw, igk, qcart, t, zbuffer)
        npwx_int = npwx_tmp
        call stop_clock('  proj_init')

        call start_clock('  proj_gemm')
        call ZGEMM('C', 'N', nbasis, nh(t)*pp%na(t), npw, one, &
                   opt_basis%elements, npw, &
                   zbuffer, npw, zero, &
                   vkb3, nbasis)
        call mp_sum(vkb3, intra_pool_comm )
        jkb = 1
        call load_from_local(pp%projs(t), 1, 1, vkb3)
        call stop_clock('  proj_gemm')
        ! save!
        call start_clock('  proj_save')
        if (MOD(q-1, npot) == my_pot_id .and. size(pp%projs(t)%dat) > 0) then
          call save_buffer(pp%projs(t)%dat, size(pp%projs(t)%dat), pp%p_unit(t), (q-1)/npot+1)
        endif
        call stop_clock('  proj_save')
      enddo 
      deallocate(zbuffer, vkb3)
      cycle
    endif

    vbs = pp%b_size(t)
    allocate(vbas(npw, vbs))
    call davcio ( vbas, npw*vbs, pp%b_unit, t, -1 )

    ! make structure factors in mixed basis (<opt_basis|S|vbas>)
    call start_clock('  make_St')
    allocate(S(nbasis, vbs, size(pp%projs(t)%dat,2)/nh(t)), Stmp(nbasis, vbs), sk(npw), zbuffer(npw, vbs))
    do a = 1, pp%na(t)
      do ig = 1, npw
        sk (ig) = eigts1 (mill(1,igk(ig)), a+pp%nt_off(t)-1) * &
                  eigts2 (mill(2,igk(ig)), a+pp%nt_off(t)-1) * &
                  eigts3 (mill(3,igk(ig)), a+pp%nt_off(t)-1)
      enddo
      forall (ig = 1:npw, j = 1:vbs) zbuffer(ig, j) = vbas(ig,j) * sk(ig)
      call ZGEMM('C', 'N', nbasis, vbs, npw, one, &
               opt_basis%elements, npw, &
               zbuffer, npw, zero, &
               Stmp, nbasis)
      call mp_sum(Stmp, intra_pool_comm)
      call g2l(pp%projs(t), 1, 1+(a-1)*nh(t), i_l, j_l, local)
      if (local) S(:,:,(j_l-1)/nh(t)+1) = Stmp
    enddo
    call stop_clock('  make_St')
    deallocate(sk, zbuffer, Stmp)

    ! construct and transform origin-centered projectors
    allocate(vkb2(vbs, nh(t)))
    allocate(zbuffer(vbs,nh(t)))
    qs2: do q = 1+my_pot_id, nq, npot
     qtmp = xq(:,q) - floor(xq(:,q))
     qcart = matmul( bg, qtmp )
     call davcio ( vkb2, vbs*nh(t), pp%v_unit, (t-1)*nq+q, -1 )
     zbuffer = vkb2

      ! apply structure factor to take to atomic position
      call start_clock('  apply_S')
#if 1
      pp%projs(t)%dat = zero
      atom: do a = 1, pp%na(t)
        call g2l(pp%projs(t), 1, 1+(a-1)*nh(t), i_l, j_l, local)
        if (.not. local) cycle
        arg = tpi * sum(qcart(:)*tau(:,a+pp%nt_off(t)-1))
        phase = CMPLX (cos(arg), - sin(arg))
        do ih = 1, nh(t)
          prefac = (0.d0, -1.d0)**nhtol(ih, t) * phase
          call zgemv('N', nbasis, vbs, prefac, &
                     S(:,:,(j_l-1)/nh(t)+1), nbasis, &
                     zbuffer(:,ih), 1, zero, &
                     pp%projs(t)%dat(1,j_l + ih-1), 1)
        enddo
      enddo atom
#endif 
!      call mp_sum(pp%projs(t)%dat, intra_pool_comm)
      call stop_clock('  apply_S')
      ! save!
      call start_clock('  proj_save')
      if (MOD(q-1, npot) == my_pot_id .and. size(pp%projs(t)%dat) > 0) then
        call save_buffer(pp%projs(t)%dat, size(pp%projs(t)%dat), pp%p_unit(t), (q-1)/npot+1)
      endif
      call stop_clock('  proj_save')
    enddo qs2
    deallocate(vkb2, zbuffer)

    deallocate(vbas)
    deallocate(S)
  enddo types

  deallocate(qg, vq, ylm, gk)
  deallocate(igk, g2kin)

  return

end subroutine build_projs_reduced


subroutine build_centered_projs(qcart, t, gk, qg, vq, ylm, igk, vkb0)
  use kinds, only : DP
  USE uspp, only : nkb, nhtol, nhtolm, indv
  use uspp_param, only : nh, nhm, lmaxkb, upf
  use gvect, only : g
  USE us, ONLY : dq, tab
  use cell_base, only : tpiba

  real(DP), intent(in) :: qcart(3)
  integer, intent(in) :: t
  real(DP), intent(out) :: gk(:,:), qg(:), vq(:), ylm(:,:)
  integer, intent(in) :: igk(:)
  real(DP), intent(out) :: vkb0(:,:)

  integer :: ig, ih, b, lm
  integer :: npw, nqx
  real(DP) :: px, ux, vx, wx, arg
  integer :: i0, i1, i2, i3

  nqx = size(tab,1)
  npw = size(gk, 2)

  ! make ylm
  do ig = 1, npw
    gk (:,ig) = qcart(:) + g(:, igk(ig) )
    qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  call ylmr2 ((lmaxkb+1)**2, npw, gk, qg, ylm)

  ! make vkb0
  forall (ig = 1:npw) qg(ig) = sqrt(qg(ig))*tpiba

  ! calculate beta in G-space using an interpolation table
  chan2: do b = 1, upf(t)%nbeta
    do ig = 1, npw
      px = qg(ig)/dq - int(qg(ig)/dq)
      ux = 1.d0 - px; vx = 2.d0 - px; wx = 3.d0 - px
      i0 = INT(qg(ig)/dq) + 1
      i1 = i0 + 1; i2 = i0 + 2; i3 = i0 + 3
      if (i0 > nqx .or. i1 > nqx .or. i2 > nqx .or. i3 > nqx) then
        vq (ig ) = 0.d0
        cycle
      endif
      vq (ig) = tab (i0, b, t) * ux * vx * wx / 6.d0 + &
                tab (i1, b, t) * px * vx * wx / 2.d0 - &
                tab (i2, b, t) * px * ux * wx / 2.d0 + &
                tab (i3, b, t) * px * ux * vx / 6.d0
    enddo

    ! add spherical harmonic part
    do ih = 1, nh(t)
      if (b .eq. indv(ih, t)) then
        lm = nhtolm(ih, t)
        forall(ig = 1:npw) vkb0(ig, ih) = vq(ig) * ylm(ig, lm)
      endif
    enddo
    enddo chan2

  return

end subroutine build_centered_projs

