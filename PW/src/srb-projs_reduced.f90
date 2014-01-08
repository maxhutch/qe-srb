#define THRESH 5 

subroutine build_projs_reduced(opt_basis, xq, nq, pp)
  USE kinds, ONLY : DP
  use constants, only : tpi
  USE srb_types, ONLY : basis, pseudop
  use srb_matrix, only : setup_desc

  use input_parameters, only : proj_tol
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
  integer :: a, b, t, ig, ih, i, j, jkb, jkb_old, lm
  integer, allocatable :: igk(:)
  real(DP), allocatable :: g2kin(:)
  real(DP) :: k_gamma(3)
  complex(DP), allocatable :: S(:,:,:), Stmp(:,:), sk(:)
  complex(DP), allocatable :: buffer(:,:), vbas(:,:), vkb(:,:), vkb1z(:,:), vkb2(:,:), vkb3(:,:), C(:,:)
  complex(DP), allocatable :: projs_l(:,:)
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:), sv(:)
  real(DP) :: px, ux, vx, wx, arg
  integer :: i0, i1, i2, i3
  complex(DP) :: phase, prefac
  complex(DP), parameter :: zero = cmplx(0.d0, kind=DP)
  complex(DP), parameter :: one = cmplx(1.d0, kind=DP)
  real(DP) :: qcart(3), qtmp(3), inv_norm
  logical ::  info
  complex(DP), allocatable :: work(:)
  real(DP), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  integer :: lwork, lrwork, liwork
  integer :: counter

  integer,save :: old_basis_length

  COMPLEX(DP), external :: ZDOTC
  REAL(DP), external :: DZNRM2

  ! initial setup
  nbasis = opt_basis%length
  npw = size(opt_basis%elements, 1)
  npwx = npw
  call mp_max(npwx, intra_pool_comm)
  nqx = size(tab,1)
  k_gamma = 0.d0
  allocate(igk(ngm), g2kin(ngm))
  call gk_sort(k_gamma(:), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
  write(*,*) allocated(pp%desc), size(pp%desc)
  if (.not. allocated(pp%desc)) then
    allocate(pp%desc(ntyp))
    allocate(pp%na(ntyp), pp%na_off(nat), pp%nt_off(ntyp))
    pp%na = 0; pp%na_off = 0
    do a = 1, nat
      pp%na(ityp(a)) = pp%na(ityp(a)) + 1
    enddo
    pp%nkb_l = 0
    do t = 1, ntyp    
      call setup_desc(pp%desc(t), nbasis, nh(t)*pp%na(t), nbasis, nh(t))
      pp%nkb_l = pp%nkb_l + pp%desc(t)%ncl 
    enddo
    jkb = 1
    ! this makes assumptions about atom ordering by type
    do t = 1, ntyp
      counter = 0
      pp%nt_off(t) = jkb
      do a = 1, nat
        if (ityp(a) /= t) cycle
        if (MOD(counter, nproc_pot) == me_pot) then
          pp%na_off(a) = jkb
          jkb = jkb + nh(t)
        else
          pp%na_off(a) = -1
        endif
        counter = counter + 1
      enddo
    enddo
    pp%ntyp = ntyp
    pp%nat  = nat
    pp%nkb  = nkb
    write(*,*) "na_off: ", pp%na_off
    write(*,*) "nt_off: ", pp%nt_off
    pp%nkb_max = pp%nkb_l
    call mp_max(pp%nkb_max, intra_pool_comm) 
  else
    ! resize nbasis
    do t = 1, ntyp
      call setup_desc(pp%desc(t), nbasis, nh(t)*pp%na(t), nbasis, nh(t))
    enddo
  endif

  if (allocated(pp%projs)) deallocate(pp%projs)
  allocate(pp%projs(nbasis,pp%nkb_l))

  allocate (ylm(npw, (lmaxkb + 1)**2), gk(3, npw), qg(npw), vq(npw))

  ! if the projector save file hasn't been opened, open one
  if (pp%projs_unit < 0) then
    pp%projs_unit = - pp%projs_unit
    call open_buffer(pp%projs_unit, 'projs', nbasis * pp%nkb_max, 1, info)
    old_basis_length = nbasis
  endif

  ! if we have a longer buffer - adjust size of file dynamically
  if (opt_basis%length > old_basis_length) then
    call close_buffer(pp%projs_unit,'delete')
    call open_buffer(pp%projs_unit, 'projs', opt_basis%length * pp%nkb_max, 1, info)
    old_basis_length = opt_basis%length
  endif

  ! if the auxillary basis isn't in file, make it
  if (pp%b_unit < 0) then 
    call start_clock('  make_vbas')
    pp%b_unit = - pp%b_unit
    ! auxillary q-grid
    call open_buffer(pp%b_unit, 'vbas', npwx*nhm*nq, 1, info)

    allocate(vbas(npw, nhm*nq))
    allocate(sv(nhm*nq))
    do t = 1, ntyp ! loop over types
      if (pp%na(t) < THRESH) cycle
      vbs = 0
      do q = 1, nq ! loop over aux q-grid
        qtmp = xq(:,q) - floor(xq(:,q))
        qcart = matmul( bg, qtmp )

        ! make ylm
        forall (ig = 1:npw)
          gk (:,ig) = qcart(:) + g(:, igk(ig) )
          qg (ig) = dot_product(gk(:, ig), gk(:,ig))
        end forall
        call ylmr2 ((lmaxkb+1)**2, npw, gk, qg, ylm)
        ! make vkb1
        forall (ig = 1:npw) qg(ig) = sqrt(qg(ig))*tpiba
        ! calculate beta in G-space using an interpolation table
        do b = 1, upf(t)%nbeta
          vq = 0.d0
          do ig = 1, npw
            px = qg(ig)/dq - int(qg(ig)/dq)
            ux = 1.d0 - px; vx = 2.d0 - px; wx = 3.d0 - px
            i0 = INT(qg(ig)/dq) + 1
            i1 = i0 + 1; i2 = i0 + 2; i3 = i0 + 3
            if (i0 > nqx .or. i1 > nqx .or. i2 > nqx .or. i3 > nqx) cycle
            vq (ig) = tab (i0, b, t) * ux * vx * wx / 6.d0 + &
                      tab (i1, b, t) * px * vx * wx / 2.d0 - &
                      tab (i2, b, t) * px * ux * wx / 2.d0 + &
                      tab (i3, b, t) * px * ux * vx / 6.d0
          enddo
          ! add spherical harmonic part
          do ih = 1, nh(t)
            if (b .eq. indv(ih, t)) then
              vbs = vbs + 1 
              lm = nhtolm(ih, t)
              forall(ig = 1:npw) vbas(ig, vbs) = vq(ig) * ylm(ig, lm)
            endif
          enddo
        enddo 
      enddo 
      ! SVD!
      allocate(C(vbs, vbs))
      call zherk('U', 'C', vbs, npw, one, &
                 vbas, npw, zero, &
                 C, vbs)
      call mp_sum(C, intra_pool_comm)
      lwork = 2 * vbs * vbs; lrwork = 2*lwork; liwork = 3+5*vbs
      allocate(work(lwork), rwork(lrwork), iwork(liwork))
      sv = 0.
      call zheevd('V', 'U', vbs, &
              C, vbs, &
              sv, &
              work, lwork, &
              rwork, lrwork, &
              iwork, liwork, stat)
      deallocate(work, rwork, iwork)
      if (stat /= 0) write(*,*) "zheevd failed: ", stat
      sv = sqrt(abs(sv))
      fq = vbs
      do ih = 1, fq
        vbs = ih
        if (abs(1 - sum(sv(fq-ih+1:fq))/sum(sv(1:fq))) < proj_tol) exit
      enddo
      allocate(buffer(npw, vbs))
      call zgemm('N', 'N', npw, vbs, fq, one, &
                 vbas, npw, &
                 C(:,fq-vbs+1:fq), fq, zero, &
                 buffer, npw)
      do ih = 1, vbs
        inv_norm = dznrm2 (npw, buffer(:,ih), 1)
        inv_norm = inv_norm * inv_norm
        call mp_sum(inv_norm, intra_pool_comm)
        inv_norm = 1.d0/sqrt(inv_norm)
        call zdscal(npw, inv_norm, buffer(:,ih), 1)
      enddo
      ! end SVD
      call save_buffer(buffer, npw*vbs, pp%b_unit, t)
      deallocate(buffer, C)
      write(*,*) "VBS: ", (1.*vbs)/nh(t), nq, t
      pp%b_size(t) = vbs
    enddo
    deallocate(vbas, sv)
    call stop_clock('  make_vbas')
  endif

  ! transform the projectors
  allocate(vkb1(npw, nhm), vkb1z(npw, nhm))
  jkb_old = 0
  types: do t = 1, ntyp ! loop over types
    if (pp%na(t) < THRESH) then
      allocate(vkb2(npw, nh(t)*pp%na(t)), vkb3(nbasis, nh(t)*pp%na(t)))
      do q = 1, nq
        call start_clock('  proj_init')
        qtmp = xq(:,q) - floor(xq(:,q))
        qcart = matmul( bg, qtmp )

        npwx_tmp = npwx_int; npwx_int = npw
        call init_us_2_srb(npw, igk, qcart, t, vkb2)
        npwx_int = npwx_tmp
        call stop_clock('  proj_init')

        call start_clock('  proj_save')
        if (t > 1 .and. MOD(q-1, npot) == my_pot_id) then
          call get_buffer(pp%projs, nbasis*pp%nkb_l, pp%projs_unit, (q-1)/npot+1)
        endif
        call stop_clock('  proj_save')

        call start_clock('  proj_gemm')
        call ZGEMM('C', 'N', nbasis, nh(t)*pp%na(t), npw, one, &
                   opt_basis%elements, npw, &
                   vkb2, npw, zero, &
                   vkb3, nbasis)
        call mp_sum(vkb3, intra_pool_comm )
        jkb = 1
        do a = 1, nat
          if (ityp(a) /= t) cycle
          if (pp%na_off(a) > 0) then
            pp%projs(:,pp%na_off(a):pp%na_off(a)+nh(t)-1) = vkb3(:,jkb:jkb+nh(t)-1)
          endif
          jkb = jkb + nh(t)
        enddo
        call stop_clock('  proj_gemm')
        ! save!
        call start_clock('  proj_save')
        if (MOD(q-1, npot) == my_pot_id) then
          call save_buffer(pp%projs, nbasis * pp%nkb_l, pp%projs_unit, (q-1)/npot+1)
        endif
        call stop_clock('  proj_save')

      enddo 
      jkb_old = jkb_old + nh(t) * pp%na(t)
      deallocate(vkb2, vkb3)
      cycle
    endif

    vbs = pp%b_size(t)
    allocate(vbas(npw, vbs))
    call get_buffer(vbas, npw*vbs, pp%b_unit, t) ! read the aux basis

    ! make structure factors in mixed basis (<opt_basis|S|vbas>)
    call start_clock('  make_St')
    allocate(S(nbasis, vbs, (nat-1)/nproc_pool+1), Stmp(nbasis, vbs), sk(npw), buffer(npw, vbs))
    do a = 1, nat
      if (ityp(a) /= t) cycle
      do ig = 1, npw
        sk (ig) = eigts1 (mill(1,igk(ig)), a) * &
                  eigts2 (mill(2,igk(ig)), a) * &
                  eigts3 (mill(3,igk(ig)), a)
      enddo
      forall (ig = 1:npw, j = 1:vbs) buffer(ig, j) = vbas(ig,j) * sk(ig)
      call ZGEMM('C', 'N', nbasis, vbs, npw, one, &
               opt_basis%elements, npw, &
               buffer, npw, zero, &
               Stmp, nbasis)
      call mp_sum(Stmp, intra_pool_comm)
      if (MOD(a-1,nproc_pool) == me_pool) S(:,:,(a-1)/nproc_pool+1) = Stmp
    enddo
    call stop_clock('  make_St')
    deallocate(sk, buffer, Stmp)

    ! construct and transform origin-centered projectors
    allocate(vkb2(vbs, nh(t)))
    qs2: do q = 1, nq
      call start_clock('  make_vkb')
      qtmp = xq(:,q) - floor(xq(:,q))
      qcart = matmul( bg, qtmp )

      ! make ylm
      do ig = 1, npw
         gk (:,ig) = qcart(:) + g(:, igk(ig) )
         qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
      enddo
      !
      call ylmr2 ((lmaxkb+1)**2, npw, gk, qg, ylm)
      do ig = 1, npw
        qg(ig) = sqrt(qg(ig))*tpiba
      enddo
      ! make vkb1
      jkb = jkb_old
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
            forall(ig = 1:npw) vkb1(ig, ih) = vq(ig) * ylm(ig, lm)
          endif
        enddo
      enddo chan2
      ! transform vkb1 into the auxillary basis
      vkb1z = cmplx(vkb1, kind=DP)
      call ZGEMM('C', 'N', vbs, nh(t), npw, one, &
                 vbas, npw, &
                 vkb1z, npw, zero, &
                 vkb2, vbs)
      call mp_sum(vkb2(:,1:nh(t)), intra_pool_comm)
      call stop_clock('  make_vkb')
      
      call start_clock('  proj_save')
      if (t > 1 .and. MOD(q-1, npot) == my_tub_id) then
        call get_buffer(pp%projs, nbasis*pp%nkb_l, pp%projs_unit, (q-1)/nproc_pool+1)
      else
        pp%projs = cmplx(0.d0, kind=DP)
      endif
      call stop_clock('  proj_save')

      ! apply structure factor to take to atomic position
      call start_clock('  apply_S')
      atom: do a = 1+me_pool, nat, nproc_pool
        if (ityp(a) /= t) cycle
        arg = tpi * sum(qcart(:)*tau(:,a))
        phase = CMPLX (cos(arg), - sin(arg))
        do ih = 1, nh(t)
          jkb = jkb + 1
          prefac = (0.d0, -1.d0)**nhtol(ih, t) * phase
          call zgemv('N', nbasis, vbs, prefac, &
                     S(:,:,(a-1)/nproc_pool+1), nbasis, &
                     vkb2(:,ih), 1, zero, &
                     pp%projs(:,pp%na_off(a) + ih), 1)
        enddo
      enddo atom
!      write(*,*) jkb_old, nh(t), na(t), nkb
      call mp_sum(pp%projs(:,jkb_old+1:jkb_old+nh(t)*pp%na(t)), intra_pool_comm)
      call stop_clock('  apply_S')

      ! save!
      call start_clock('  proj_save')
      if (MOD(q-1, npot) == my_pot_id) then
        call save_buffer(pp%projs, nbasis * pp%nkb_l, pp%projs_unit, (q-1)/nproc_pool+1)
      endif
      call stop_clock('  proj_save')

    enddo qs2
    jkb_old = jkb_old + nh(t) * pp%na(t)
    deallocate(vbas)
    deallocate(S)
    deallocate(vkb2)

  enddo types

  deallocate(vkb1, vkb1z, qg, vq, ylm, gk)
  deallocate(igk, g2kin)

  return

end subroutine build_projs_reduced

