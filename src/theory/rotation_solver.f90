module rotation_uniform
  use analysis_mod, only: mass_radius
  use precision_mod, only: wp
  use para_mod, only: active_theory, THEORY_GR, &
                      rho, gama, alpha, ww, omg, sphi, &
                      r_e, r_ratio, h_center, SDIV, MDIV, &
                      B_coup, mphi_r, F_j, &
                      sphi_c, sphi_m, Omega_c, Omega_e, &
                      Fmax_h, n_of_relaxation_steps, timing, solver_type
  use spin_workspace, only: dif, &
    target_rho, target_gama, target_ww, target_sphi, &
    metric_method, scalar_method, &
    dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, &
    dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
    d2g_ss_cache, d2g_mm_cache, e_rsm_cache, &
    sgp_term_2d_cache, sin_theta_2d_cache, sgp_2d_cache, &
    D2_metric_rho, D2_metric_omega, &
    allocate_workspace, deallocate_workspace
  use spin_integration, only: get_all_targets, update_alpha_potential, output_helper
  use spin_relaxation, only: relaxation
  use spin_updates, only: reset_uryu_peak_cache, &
    update_equatorial_radius, update_angular_velocity, update_eos_and_velocity
  implicit none
contains

subroutine rotation_solver
  implicit none
  integer :: n_of_it
  real(wp) :: r_e_old, r_e_new, r_e_new_sq
  real(wp) :: gama_pole_h, gama_center_h, gama_equator_h
  real(wp) :: rho_pole_h, rho_center_h, rho_equator_h, ww_equator_h
  real(wp) :: sphi_pole_h, sphi_center_h, sphi_equator_h
  real(wp) :: root_mphi_re, sqrt_B_coup
  real(wp) :: t0, t1, dt_alpha, dt_relaxation
  logical :: zero_scalar_mode

  ! ---------------------------------------------------------------
  ! Setup phase
  ! ---------------------------------------------------------------
  call reset_uryu_peak_cache()
  zero_scalar_mode = merge(.true., .false., active_theory == THEORY_GR)
  sqrt_B_coup = sqrt(B_coup)
  dif = 1.e0_wp
  n_of_it = 0
  r_e_new = r_e
  r_e_new_sq = r_e_new**2

  if ( maxval(sphi*sqrt_B_coup) < 1.e-3_wp ) sphi = sphi * 10.e0_wp
  if (zero_scalar_mode) sphi = 0.e0_wp
  if ( any(isnan(sphi)) ) stop "NaN found in sphi"

  call allocate_workspace
  if (trim(solver_type) == "uryu") call reset_uryu_peak_cache()
  !write(*,*) h_center, r_ratio
  ! ---------------------------------------------------------------
  ! Main iteration loop
  ! ---------------------------------------------------------------
  iteration: block
    real(wp), allocatable :: rho_prev_iter(:,:), gama_prev_iter(:,:)
    real(wp), allocatable :: ww_prev_iter(:,:), sphi_prev_iter(:,:)
    real(wp) :: drho_prev, dgama_prev, dww_prev, dsphi_prev, dre_prev

    allocate(rho_prev_iter(SDIV,MDIV), source=rho)
    allocate(gama_prev_iter(SDIV,MDIV), source=gama)
    allocate(ww_prev_iter(SDIV,MDIV), source=ww)
    allocate(sphi_prev_iter(SDIV,MDIV), source=sphi)
    drho_prev = -1.e0_wp; dgama_prev = -1.e0_wp; dww_prev = -1.e0_wp
    dsphi_prev = -1.e0_wp; dre_prev = -1.e0_wp

    do while( dif > 1.e-7_wp .or. n_of_it < 2 )
      if (zero_scalar_mode) sphi = 0.e0_wp
      sphi_m = maxval( sphi(:,1) * sqrt_B_coup )
      call rescale_metric(r_e_new_sq)

      ! --- Update equatorial radius ---
      r_e_old = r_e_new
      call update_equatorial_radius( r_e_old, &
        r_e_new, dif, &
        sphi_pole_h, gama_pole_h, rho_pole_h, &
        gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h, &
        sphi_center_h, gama_center_h, rho_center_h )
      r_e_new_sq = r_e_new**2

      ! --- Solve for metric and scalar field targets ---
      call update_angular_velocity(r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, &
        sphi_pole_h, sphi_equator_h, ww_equator_h)

      ! Also rescale back metric potentials here (except for omega / ww)
      call update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, &
                                   sgp_term_2d_cache, sin_theta_2d_cache, sgp_2d_cache)
      root_mphi_re = sqrt(mphi_r * r_e_new_sq)

      call get_all_targets(r_e_new, root_mphi_re, &
        target_rho, target_gama, target_ww, target_sphi)

      ! --- Relaxation iteration ---
      if (timing) call cpu_time(t0)
      call relaxation(target_rho, target_gama, target_ww, target_sphi, root_mphi_re, n_of_it, dif)
      if (timing) then
        call cpu_time(t1); dt_relaxation = t1 - t0; call cpu_time(t0)
      end if

      ! --- Update metric potential (alpha) ---
      if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
        call impose_rigid_rotation()
      else
        call update_alpha_potential(r_e_new, dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, &
            dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_rsm_cache)
      endif
      if (timing) then
        call cpu_time(t1); dt_alpha = t1 - t0; call cpu_time(t0)
      end if
      if (timing) write(*,'(A,7(1X,ES12.5))') "Relaxation + Alpha: ", dt_relaxation, dt_alpha

      ! Rescale back omega / ww
      omg= omg / r_e_new; ww = ww / r_e_new
      
      ! --- Diagnostics and contraction ratios ---
      diagnostics: block
        real(wp) :: drho_norm, dgama_norm, dww_norm, dsphi_norm, dre_norm
        real(wp) :: q_rho, q_gama, q_ww, q_sphi, q_re

        drho_norm  = field_update_norm(rho,  rho_prev_iter)
        dgama_norm = field_update_norm(gama, gama_prev_iter)
        dww_norm   = field_update_norm(ww,   ww_prev_iter)
        dsphi_norm = field_update_norm(sphi, sphi_prev_iter)
        dre_norm   = abs(r_e_new - r_e_old)

        q_rho  = contraction_ratio(drho_norm,  drho_prev)
        q_gama = contraction_ratio(dgama_norm, dgama_prev)
        q_ww   = contraction_ratio(dww_norm,   dww_prev)
        q_sphi = contraction_ratio(dsphi_norm, dsphi_prev)
        q_re   = contraction_ratio(dre_norm,   dre_prev)

        if ( n_of_it > 50 .and. mod(n_of_it,50)==0 ) then
          write(*,'(A,i4,A,es10.3,A,2es12.4,A,1X,A,1X,A,A,5es12.4)') &
            'it= ', n_of_it, ', dif:', dif, &
            ' sphi:', sphi_center_h*r_e_old*sqrt_B_coup, sphi_m, &
            ' ', trim(metric_method), trim(scalar_method), &
            '  |', q_rho, q_gama, q_re, q_ww, q_sphi
        endif

        drho_prev  = drho_norm
        dgama_prev = dgama_norm
        dww_prev   = dww_norm
        dsphi_prev = dsphi_norm
        dre_prev   = dre_norm
      end block diagnostics

      ! --- Update previous-iteration state ---
      rho_prev_iter  = rho
      gama_prev_iter = gama
      ww_prev_iter   = ww
      sphi_prev_iter = sphi

      n_of_it = n_of_it + 1
      if ( n_of_it > 2000 ) stop "Probably won't converge"
    enddo

    deallocate(rho_prev_iter, gama_prev_iter, ww_prev_iter, sphi_prev_iter)
  end block iteration
  !write(*,*) n_of_it
  n_of_relaxation_steps = n_of_relaxation_steps + n_of_it

  ! ---------------------------------------------------------------
  ! Finalization and output
  ! ---------------------------------------------------------------
  Omega_c  = Omega_c / r_e_new
  Omega_e  = Omega_e / r_e_new
  r_e      = r_e_new
  Fmax_h   = maxval(F_j(:,1))
  if (zero_scalar_mode) then
    sphi = 0.e0_wp
    sphi_c = 0.e0_wp
    sphi_m = 0.e0_wp
  else
    sphi_c = sphi(1,1) * sqrt_B_coup
    sphi_m = maxval( sphi(:,1) * sqrt_B_coup )
  end if
  call mass_radius()

  call output_helper(D2_metric_rho, D2_metric_omega)
  
  call deallocate_workspace

contains

  pure real(wp) function field_update_norm(a, b) result(val)
    real(wp), intent(in) :: a(:,:), b(:,:)
    val = sqrt(sum((a - b)**2) / real(size(a), wp))
  end function field_update_norm

  pure real(wp) function contraction_ratio(curr, prev) result(val)
    real(wp), intent(in) :: curr, prev
    if (prev > 0.e0_wp) then
      val = curr / prev
    else
      val = -1.e0_wp
    end if
  end function contraction_ratio

  subroutine rescale_metric(factor)
    implicit none
    real(wp), intent(in) :: factor
    real(wp) :: inv_factor, sqrt_factor
    inv_factor = 1.e0_wp / factor
    sqrt_factor = sqrt(factor)

    rho   = rho   * inv_factor
    gama  = gama  * inv_factor
    alpha = alpha * inv_factor
    ww    = ww    * sqrt_factor
    omg   = omg   * sqrt_factor
    sphi  = sphi  / sqrt_factor
  end subroutine rescale_metric

  subroutine impose_rigid_rotation()
    rho  = spread(rho(:,1),  dim=2, ncopies=MDIV)
    sphi = spread(sphi(:,1), dim=2, ncopies=MDIV)
    gama = spread(gama(:,1), dim=2, ncopies=MDIV)
    ww   = spread(ww(:,1),   dim=2, ncopies=MDIV)
    omg  = spread(omg(:,1),  dim=2, ncopies=MDIV)
    alpha= spread( (gama(:,1) - rho(:,1)) * 0.5e0_wp, dim=2, ncopies=MDIV )
  end subroutine impose_rigid_rotation

end subroutine rotation_solver
end module rotation_uniform
