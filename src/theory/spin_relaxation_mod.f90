module spin_relaxation_mod
  use precision_mod, only: wp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use para_mod, only: wp, SDIV, MDIV, s_gp, &
                      rho, gama, ww, sphi, &
                      has_scalar, mphi_r, timing
  use spin_workspace_mod, only: metric_method, scalar_method
  use anderson_optimized_mod, only: anderson_accel_optimized
  use aitken_mod, only: aitken_delta2, aitken_reset
  implicit none

contains

  subroutine relaxation(target_rho, target_gama, target_ww, target_sphi, root_mphi_re, n_of_it, dif)
    real(wp), intent(in) :: target_rho(SDIV,MDIV), target_gama(SDIV,MDIV)
    real(wp), intent(in) :: target_ww(SDIV,MDIV),  target_sphi(SDIV,MDIV)
    real(wp), intent(in) :: root_mphi_re, dif
    integer, intent(in) :: n_of_it

    real(wp), parameter :: W_PICARD = 0.7_wp, W_MIX_MIN = 5.e-2_wp, W_MIX_DECAY = 7.5e-1_wp
    real(wp), parameter :: PICARD_THRESH = 5.e-1_wp, CHEB_THRESH = 1.e-1_wp
    real(wp), parameter :: RHO_LOCK_TOL = 2.e-2_wp, CYCLE_TOL = 1.e-3_wp
    real(wp), parameter :: DIF_TAIL_THRESH = 1.e-5_wp, ANDERSON_TAIL_THRESH = 1.e-4_wp
    integer,  parameter :: M_HIST = 3, PICARD_STALL_LIMIT = 8
    integer,  parameter :: N_CHEB = 3, N_ANDERSON = 5, N_RHO_LOCK = 3, N_CHEB_FAIL = 3
    integer,  parameter :: N_ANDERSON_FAIL = 2, ANDERSON_SEQ_MAX = 2
    integer,  parameter :: AITKEN_COOLDOWN = 3, ANDERSON_COOLDOWN = 30
    integer,  parameter :: N_CYCLE_LAGS = 8, CYCLE_LAG_MAX = 200
    integer,  parameter :: CYCLE_LAGS(N_CYCLE_LAGS) = [6, 12, 18, 24, 50, 100, 150, 200]

    real(wp), save :: cheb_w = 1.0_wp, rho_spec_est = 0.7_wp, rho_spec_prev = 0.7_wp
    real(wp), save :: prev_dif_local = -1.0_wp, dif_prev = -1.0_wp
    real(wp), save :: picard_best_dif = huge(1.0_wp), w_mix = W_PICARD
    real(wp), save :: dif_min_seen = huge(1.0_wp), dif_ring(CYCLE_LAG_MAX) = -1.0_wp
    integer,  save :: n_consec_decrease = 0, n_picard_stall = 0, n_rho_locked = 0
    integer,  save :: n_cheb_hurt = 0, n_anderson_hurt = 0, n_anderson_metric = 0, n_anderson_scalar = 0
    integer,  save :: n_aitken_cooldown = 0, n_anderson_cooldown = 0, ring_pos = 0
    logical,  save :: last_was_cheb = .false., last_was_aitken = .false., last_was_anderson = .false.
    real(wp), allocatable, save :: prev_rho(:,:), prev_gama(:,:), prev_ww(:,:), prev_sphi(:,:)
    real(wp), allocatable, save :: hist_rho(:,:,:), hist_gama(:,:,:), hist_ww(:,:,:), hist_sphi(:,:,:)

    real(wp) :: x_k_rho(SDIV,MDIV), x_k_gama(SDIV,MDIV), x_k_ww(SDIV,MDIV), x_k_sphi(SDIV,MDIV)
    real(wp) :: t0, t1, dt_prep, dt_aitken, dt_metric, dt_scalar, dt_floor
    integer, parameter :: timing_calls = 5
    integer, save :: relaxation_call_count = 0
    real(wp), save :: sum_dt_prep = 0.0_wp, sum_dt_aitken = 0.0_wp, sum_dt_metric = 0.0_wp
    real(wp), save :: sum_dt_scalar = 0.0_wp, sum_dt_floor = 0.0_wp
    logical  :: aitken_fired, use_picard, use_chebys

    if (timing) then
      relaxation_call_count = relaxation_call_count + 1
      dt_prep = 0.0_wp
      dt_aitken = 0.0_wp
      dt_metric = 0.0_wp
      dt_scalar = 0.0_wp
      dt_floor = 0.0_wp
      call cpu_time(t0)
    end if

    call ensure_allocated
    if (n_of_it == 0) call reset_accel_state
    call update_diagnostics

    if (timing) then
      call cpu_time(t1); dt_prep = t1 - t0; call cpu_time(t0)
    end if

    call try_aitken

    if (timing) then
      call cpu_time(t1); dt_aitken = t1 - t0; call cpu_time(t0)
    end if

    if (.not. aitken_fired) then
      call select_tier
      call apply_metric_update
    end if

    if (timing) then
      call cpu_time(t1); dt_metric = t1 - t0; call cpu_time(t0)
    end if

    last_was_cheb     = (metric_method == 'Chebys')
    last_was_aitken   = aitken_fired
    last_was_anderson = (metric_method == 'Anders')
    if (n_aitken_cooldown > 0)  n_aitken_cooldown  = n_aitken_cooldown  - 1
    if (n_anderson_cooldown > 0) n_anderson_cooldown = n_anderson_cooldown - 1

    call check_divergence

    if (.not. has_scalar) then
      if (timing) then
        sum_dt_prep = sum_dt_prep + dt_prep
        sum_dt_aitken = sum_dt_aitken + dt_aitken
        sum_dt_metric = sum_dt_metric + dt_metric
        sum_dt_scalar = sum_dt_scalar + dt_scalar
        sum_dt_floor = sum_dt_floor + dt_floor
        if (mod(relaxation_call_count, timing_calls) == 0) then
          write(*,'(A,I0,A)') 'relaxation running avg over ', relaxation_call_count, ':'
          write(*,'(A,1X,ES12.5)') '  prep',   sum_dt_prep   / relaxation_call_count
          write(*,'(A,1X,ES12.5)') '  aitken', sum_dt_aitken / relaxation_call_count
          write(*,'(A,1X,ES12.5)') '  metric', sum_dt_metric / relaxation_call_count
          write(*,'(A,1X,ES12.5)') '  scalar', sum_dt_scalar / relaxation_call_count
          write(*,'(A,1X,ES12.5)') '  floor',  sum_dt_floor  / relaxation_call_count
          write(*,'(A,1X,ES12.5)') '  total', &
            (sum_dt_prep + sum_dt_aitken + sum_dt_metric + sum_dt_scalar + sum_dt_floor) / relaxation_call_count
        end if
      end if
      return
    end if

    if (.not. aitken_fired) call apply_scalar_update
    where(ieee_is_nan(sphi)) sphi = 0.0_wp

    if (timing) then
      call cpu_time(t1); dt_scalar = t1 - t0; call cpu_time(t0)
    end if

    call enforce_sphi_floor

    if (timing) then
      call cpu_time(t1); dt_floor = t1 - t0
      sum_dt_prep = sum_dt_prep + dt_prep
      sum_dt_aitken = sum_dt_aitken + dt_aitken
      sum_dt_metric = sum_dt_metric + dt_metric
      sum_dt_scalar = sum_dt_scalar + dt_scalar
      sum_dt_floor = sum_dt_floor + dt_floor
      if (mod(relaxation_call_count, timing_calls) == 0) then
        write(*,'(A,I0,A)') 'relaxation running avg over ', relaxation_call_count, ':'
        write(*,'(A,1X,ES12.5)') '  prep',   sum_dt_prep   / relaxation_call_count
        write(*,'(A,1X,ES12.5)') '  aitken', sum_dt_aitken / relaxation_call_count
        write(*,'(A,1X,ES12.5)') '  metric', sum_dt_metric / relaxation_call_count
        write(*,'(A,1X,ES12.5)') '  scalar', sum_dt_scalar / relaxation_call_count
        write(*,'(A,1X,ES12.5)') '  floor',  sum_dt_floor  / relaxation_call_count
        write(*,'(A,1X,ES12.5)') '  total', &
          (sum_dt_prep + sum_dt_aitken + sum_dt_metric + sum_dt_scalar + sum_dt_floor) / relaxation_call_count
      end if
    end if

  contains

    subroutine ensure_allocated
      if (allocated(prev_rho) .and. size(prev_rho,1) == SDIV .and. size(prev_rho,2) == MDIV) return
      if (allocated(prev_rho)) deallocate(prev_rho, prev_gama, prev_ww, prev_sphi)
      if (allocated(hist_rho)) deallocate(hist_rho, hist_gama, hist_ww, hist_sphi)
      allocate(prev_rho(SDIV,MDIV),  source=rho)
      allocate(prev_gama(SDIV,MDIV), source=gama)
      allocate(prev_ww(SDIV,MDIV),   source=ww)
      allocate(prev_sphi(SDIV,MDIV), source=sphi)
      allocate(hist_rho(SDIV,MDIV,M_HIST),  source=0.0_wp)
      allocate(hist_gama(SDIV,MDIV,M_HIST), source=0.0_wp)
      allocate(hist_ww(SDIV,MDIV,M_HIST),   source=0.0_wp)
      allocate(hist_sphi(SDIV,MDIV,M_HIST), source=0.0_wp)
      cheb_w = 1.0_wp; rho_spec_est = 0.7_wp; rho_spec_prev = 0.7_wp
      prev_dif_local = -1.0_wp; dif_prev = -1.0_wp
      n_consec_decrease = 0; n_picard_stall = 0; n_rho_locked = 0; n_cheb_hurt = 0
      n_anderson_hurt = 0; n_anderson_metric = 0; n_anderson_scalar = 0
      n_aitken_cooldown = 0; n_anderson_cooldown = 0
      last_was_cheb = .false.; last_was_aitken = .false.; last_was_anderson = .false.
      picard_best_dif = huge(1.0_wp); w_mix = W_PICARD
      dif_min_seen = huge(1.0_wp); dif_ring = -1.0_wp; ring_pos = 0
    end subroutine ensure_allocated

    subroutine reset_accel_state
      prev_rho = rho;  prev_gama = gama;  prev_ww = ww;  prev_sphi = sphi
      hist_rho = 0.0_wp; hist_gama = 0.0_wp; hist_ww = 0.0_wp; hist_sphi = 0.0_wp
      cheb_w = 1.0_wp; rho_spec_est = 0.7_wp; rho_spec_prev = 0.7_wp
      prev_dif_local = -1.0_wp; dif_prev = -1.0_wp
      n_consec_decrease = 0; n_picard_stall = 0; n_rho_locked = 0; n_cheb_hurt = 0
      n_anderson_hurt = 0; n_anderson_metric = 0; n_anderson_scalar = 0
      n_aitken_cooldown = 0; n_anderson_cooldown = 0
      last_was_cheb = .false.; last_was_aitken = .false.; last_was_anderson = .false.
      picard_best_dif = huge(1.0_wp); w_mix = W_PICARD
      dif_min_seen = huge(1.0_wp); dif_ring = -1.0_wp; ring_pos = 0
      call aitken_reset
    end subroutine reset_accel_state

    subroutine update_diagnostics
      real(wp) :: rho_obs, dif_ref
      integer :: lag_idx, hist_pos
      logical :: cycle_hit

      dif_min_seen = min(dif_min_seen, dif)
      if (dif_min_seen < 1.e-3_wp .and. dif > dif_min_seen * 1.e2_wp .and. n_of_it > 100) &
        w_mix = max(W_MIX_MIN, W_MIX_DECAY * w_mix)

      ring_pos = mod(ring_pos, CYCLE_LAG_MAX) + 1
      cycle_hit = .false.
      if (dif < DIF_TAIL_THRESH) then
        do lag_idx = 1, N_CYCLE_LAGS
          hist_pos = modulo(ring_pos - CYCLE_LAGS(lag_idx) - 1, CYCLE_LAG_MAX) + 1
          dif_ref = dif_ring(hist_pos)
          if (dif_ref > 0.0_wp .and. abs(dif - dif_ref) <= CYCLE_TOL * max(dif, dif_ref)) then
            cycle_hit = .true.; exit
          end if
        end do
      end if
      if (cycle_hit) then
        w_mix = max(W_MIX_MIN, 0.5_wp * w_mix)
        n_anderson_cooldown = ANDERSON_COOLDOWN
        n_anderson_hurt = max(n_anderson_hurt, N_ANDERSON_FAIL)
        n_anderson_metric = 0; n_anderson_scalar = 0
        hist_rho = 0.0_wp; hist_gama = 0.0_wp; hist_ww = 0.0_wp; hist_sphi = 0.0_wp
      end if
      dif_ring(ring_pos) = dif

      if (dif_prev > 0.0_wp) then
        if (last_was_aitken .and. dif > dif_prev) then
          n_aitken_cooldown = AITKEN_COOLDOWN; n_rho_locked = 0
        end if
        if (dif < dif_prev) then
          n_consec_decrease = n_consec_decrease + 1
          if (last_was_cheb) n_cheb_hurt = max(0, n_cheb_hurt - 1)
          if (last_was_anderson) n_anderson_hurt = max(0, n_anderson_hurt - 1)
        else if (last_was_anderson) then
          n_consec_decrease = 0
          n_anderson_hurt = n_anderson_hurt + 1
          n_anderson_cooldown = ANDERSON_COOLDOWN
          n_anderson_metric = 0; n_anderson_scalar = 0
        else if (last_was_cheb) then
          n_cheb_hurt = n_cheb_hurt + 1
        else
          n_consec_decrease = 0
        end if
      end if
      dif_prev = dif

      if (abs(rho_spec_est - rho_spec_prev) < RHO_LOCK_TOL .and. n_consec_decrease >= N_ANDERSON) then
        n_rho_locked = n_rho_locked + 1
      else
        n_rho_locked = 0
      end if
      rho_spec_prev = rho_spec_est
    end subroutine update_diagnostics

    subroutine try_aitken
      call aitken_delta2(rho, gama, ww, sphi, prev_rho, prev_gama, prev_ww, prev_sphi, &
                        merge(n_rho_locked, 0, n_aitken_cooldown <= 0), N_RHO_LOCK, has_scalar, aitken_fired)
      if (.not. aitken_fired) return
      metric_method = 'Aitken'
      if (has_scalar) scalar_method = 'Aitken'
      call reset_picard_state
      n_anderson_metric = 0; n_anderson_scalar = 0
    end subroutine try_aitken

    subroutine select_tier
      logical :: allow_anderson

      if (dif <= ANDERSON_TAIL_THRESH) then
        use_picard = .true.
        use_chebys = .false.
        return
      end if

      allow_anderson = (n_anderson_cooldown <= 0) .and. &
                       (n_anderson_hurt < N_ANDERSON_FAIL) .and. &
                       (n_anderson_metric < ANDERSON_SEQ_MAX)

      use_picard = (dif > PICARD_THRESH .or. n_consec_decrease < N_CHEB .or. n_anderson_cooldown > 0)
      use_chebys = (.not. use_picard) .and. &
        ((dif > CHEB_THRESH .or. n_consec_decrease < N_ANDERSON .or. .not. allow_anderson) .and. &
         n_cheb_hurt < N_CHEB_FAIL)
      if (.not. allow_anderson .and. .not. use_picard .and. .not. use_chebys) use_picard = .true.
    end subroutine select_tier

    subroutine apply_metric_update
      real(wp) :: rho_obs
      if (use_picard) then
        metric_method = 'Picard'
        call update_spectral_radius(0.6_wp, 0.4_wp, 1.e-1_wp, 9.9e-1_wp)
        call adapt_picard_weight
        cheb_w = 1.0_wp
        rho  = (1.0_wp - w_mix) * rho  + w_mix * target_rho
        gama = (1.0_wp - w_mix) * gama + w_mix * target_gama
        ww   = (1.0_wp - w_mix) * ww   + w_mix * target_ww
        prev_rho = rho;  prev_gama = gama;  prev_ww = ww
        n_anderson_metric = 0
      else if (use_chebys) then
        metric_method = 'Chebys'
        call reset_picard_state
        if (n_consec_decrease >= N_CHEB) &
          call update_spectral_radius(0.7_wp, 0.3_wp, 3.e-1_wp, 9.8e-1_wp)
        prev_dif_local = dif
        cheb_w = 1.0_wp / (1.0_wp - (rho_spec_est**2 / 4.0_wp) * cheb_w)
        cheb_w = min(cheb_w, 2.0_wp - 2.0_wp*epsilon(cheb_w))
        x_k_rho = rho;  x_k_gama = gama;  x_k_ww = ww
        rho  = prev_rho  + cheb_w * (target_rho  - prev_rho)
        gama = prev_gama + cheb_w * (target_gama - prev_gama)
        ww   = prev_ww   + cheb_w * (target_ww   - prev_ww)
        prev_rho = x_k_rho;  prev_gama = x_k_gama;  prev_ww = x_k_ww
        n_anderson_metric = 0
      else
        metric_method = 'Anders'
        call reset_picard_state
        cheb_w = 1.0_wp
        n_anderson_metric = n_anderson_metric + 1
        call anderson_accel_optimized(rho,  target_rho,  hist_rho,  n_anderson_metric, M_HIST)
        call anderson_accel_optimized(gama, target_gama, hist_gama, n_anderson_metric, M_HIST)
        call anderson_accel_optimized(ww,   target_ww,   hist_ww,   n_anderson_metric, M_HIST)
        prev_rho = rho;  prev_gama = gama;  prev_ww = ww
      end if
    end subroutine apply_metric_update

    subroutine apply_scalar_update
      if (use_picard) then
        scalar_method = 'Picard'
        sphi = (1.0_wp - w_mix) * sphi + w_mix * target_sphi
        prev_sphi = sphi; n_anderson_scalar = 0
      else if (use_chebys) then
        scalar_method = 'Chebys'
        x_k_sphi = sphi
        sphi = prev_sphi + cheb_w * (target_sphi - prev_sphi)
        prev_sphi = x_k_sphi; n_anderson_scalar = 0
      else
        scalar_method = 'Anders'
        n_anderson_scalar = n_anderson_scalar + 1
        call anderson_accel_optimized(sphi, target_sphi, hist_sphi, n_anderson_scalar, M_HIST)
        prev_sphi = sphi
      end if
    end subroutine apply_scalar_update

    subroutine reset_picard_state
      n_picard_stall = 0; picard_best_dif = huge(1.0_wp)
    end subroutine reset_picard_state

    subroutine update_spectral_radius(keep, blend, lo, hi)
      real(wp), intent(in) :: keep, blend, lo, hi
      real(wp) :: rho_obs
      if (prev_dif_local > 0.0_wp .and. dif > 0.0_wp) then
        rho_obs      = min(hi, max(lo, dif / prev_dif_local))
        rho_spec_est = keep * rho_spec_est + blend * rho_obs
      end if
      prev_dif_local = dif
    end subroutine update_spectral_radius

    subroutine adapt_picard_weight
      if (dif <= 1.e-3_wp) return
      if (dif < picard_best_dif) then
        picard_best_dif = dif; n_picard_stall = 0
        w_mix = min(w_mix + 5.e-2_wp, W_PICARD)
      else
        n_picard_stall = n_picard_stall + 1
        if (n_picard_stall >= PICARD_STALL_LIMIT) then
          w_mix = max(W_MIX_MIN, W_MIX_DECAY * w_mix); n_picard_stall = 0
        end if
      end if
    end subroutine adapt_picard_weight

    subroutine check_divergence
      if (abs(rho(2,1))>100.0_wp .or. abs(gama(2,1))>300.0_wp .or. abs(ww(2,1))>100.0_wp &
          .or. abs(sphi(2,1))>10.0_wp) then
        write(*,"(i5,4es18.9)") n_of_it, rho(2,1), gama(2,1), ww(2,1), sphi(2,1)
        error stop "something diverged"
      end if
    end subroutine check_divergence

    subroutine enforce_sphi_floor
      real(wp) :: sphi_floor_decay(SDIV)
      integer :: s, m
      if (abs(mphi_r) < epsilon(mphi_r)) return
      sphi_floor_decay(1) = 1.0_wp
      sphi_floor_decay(2:SDIV) = exp(root_mphi_re * ( &
        s_gp(1:SDIV-1) / (1.0_wp - s_gp(1:SDIV-1)) - &
        s_gp(2:SDIV)   / (1.0_wp - s_gp(2:SDIV)) ))
      do m = 1, MDIV
        do s = 1, (SDIV+1)/2
          if (sphi(s,m) < 1.e-4_wp) sphi(s,m) = 0.0_wp
        end do
        do s = (SDIV+1)/2, SDIV
          if (sphi(s,m) < 0.0_wp) sphi(s,m) = sphi(s-1,m) * sphi_floor_decay(s)
        end do
      end do
    end subroutine enforce_sphi_floor

  end subroutine relaxation

end module spin_relaxation_mod
