module spin_helper
  use para_mod
  use brent_mod, only : find_omege_e, zbrent_rot
  use toolkit_mod, only : deriv_sm, deriv_s, deriv_m, besseli, besselk, interp, interp_log_h_to_p, interp_log_p_to_e
  use simpson_mod
  use nag_compat_mod, only : d01gaf
  implicit none
  real(wp) :: mphi_tran = 1.e-11_wp
  real(wp) :: dif
  ! --- Workspace arrays, moved to module scope ---
  real(wp), allocatable, target :: e_gsm_cache(:,:), e_rsm_cache(:,:), e2alpha_r2_cache(:,:)
  real(wp), allocatable, target :: Acoup_cache(:,:), Acoup4_cache(:,:)
  real(wp), allocatable, target :: dg_s_cache(:,:), dg_m_cache(:,:), d2g_ss_cache(:,:), d2g_mm_cache(:,:)
  real(wp), allocatable, target :: dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
  real(wp), allocatable, target :: ds_s_cache(:,:), ds_m_cache(:,:) 
  real(wp), allocatable, target :: mr_cache(:), besseli_cache(:,:), besselk_cache(:,:), wfac_cache(:)
  real(wp), allocatable :: radial_quad_weights(:)
  real(wp), allocatable :: angular_quad_weights(:)
  real(wp), allocatable, target :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
  real(wp), allocatable, target :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
  real(wp), allocatable, target :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
  real(wp), allocatable, target :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)

contains
  pure function deriv_s_vec(f) result(df_ds)
    use para_mod, only : SDIV, MDIV, DS
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_ds
    integer :: m
    df_ds(1,:) = (f(2,:) - f(1,:)) / DS
    df_ds(SDIV,:) = (f(SDIV,:) - f(SDIV-1,:)) / DS
    do m = 1, MDIV
      df_ds(2:SDIV-1,m) = (f(3:SDIV,m) - f(1:SDIV-2,m)) / (2.e0_wp * DS)
    end do
  end function deriv_s_vec

  pure function deriv_m_vec(f) result(df_dm)
    use para_mod, only : SDIV, MDIV, DM
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_dm
    integer :: s
    df_dm(:,1) = (f(:,2) - f(:,1)) / DM
    df_dm(:,MDIV) = (f(:,MDIV) - f(:,MDIV-1)) / DM
    do s = 1, SDIV
      df_dm(s,2:MDIV-1) = (f(s,3:MDIV) - f(s,1:MDIV-2)) / (2.e0_wp * DM)
    end do
  end function deriv_m_vec

  function deriv_sm_vec(f) result(df_dsm)
    use para_mod, only : SDIV, MDIV, DS, DM
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_dsm
    integer :: s, m

    ! Central difference for the interior
    df_dsm(2:SDIV-1, 2:MDIV-1) = (f(3:SDIV, 3:MDIV) - f(1:SDIV-2, 3:MDIV) &
            - f(3:SDIV, 1:MDIV-2) + f(1:SDIV-2, 1:MDIV-2)) / (4.e0_wp * DS * DM)

    ! Forward/backward differences for boundaries
    do s = 1, SDIV
        df_dsm(s, 1) = (deriv_s(f, s, 2) - deriv_s(f, s, 1)) / DM
        df_dsm(s, MDIV) = (deriv_s(f, s, MDIV) - deriv_s(f, s, MDIV-1)) / DM
    end do

    do m = 1, MDIV
        df_dsm(1, m) = (deriv_m(f, 2, m) - deriv_m(f, 1, m)) / DS
        df_dsm(SDIV, m) = (deriv_m(f, SDIV, m) - deriv_m(f, SDIV-1, m)) / DS
    end do
  end function deriv_sm_vec

  subroutine update_equatorial_radius(r_e_old, r_e_new, dif, sphi_pole_h, gama_pole_h, rho_pole_h, &
                                      gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h, &
                                      sphi_center_h, gama_center_h, rho_center_h)
    real(wp), intent(in)    :: r_e_old
    real(wp), intent(out)   :: r_e_new, dif
    real(wp), intent(out)   :: sphi_pole_h, gama_pole_h, rho_pole_h
    real(wp), intent(out)   :: gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h
    real(wp), intent(out)   :: sphi_center_h, gama_center_h, rho_center_h
    real(wp) :: s_p, r_e_new_sq, grgr
    real(wp), dimension(SDIV) :: gama_mu_1, gama_mu_0, rho_mu_1, rho_mu_0, ww_mu_0, sphi_mu_0, sphi_mu_1
    integer :: s

    rho_mu_0  = rho(:,1);    gama_mu_0 = gama(:,1);    sphi_mu_0 = sphi(:,1);    ww_mu_0 = ww(:,1)
    rho_mu_1  = rho(:,MDIV); gama_mu_1 = gama(:,MDIV); sphi_mu_1 = sphi(:,MDIV)
    
    s_p = r_ratio**(1.e0_wp/dble(s_pwr)) / ( 1.e0_wp + r_ratio**(1.e0_wp/dble(s_pwr)) ) 
    call interp(s_gp, sphi_mu_1, SDIV, s_p, sphi_pole_h   )
    call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole_h   )
    call interp(s_gp, rho_mu_1,  SDIV, s_p, rho_pole_h    )
    call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator_h)
    call interp(s_gp, rho_mu_0,  SDIV, s_e, rho_equator_h )
    call interp(s_gp, ww_mu_0,   SDIV, s_e, ww_equator_h  )
    call interp(s_gp, sphi_mu_0, SDIV, s_e, sphi_equator_h)
    sphi_center_h = sphi(1,1)
    gama_center_h = gama(1,1)
    rho_center_h  = rho(1,1)
    
    grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h
    if (active_theory /= THEORY_GR) then
      grgr = grgr + B_coup / 2.e0_wp * ( sphi_center_h**2 - sphi_pole_h**2 )
    end if

    r_e_new_sq = ( 2.e0_wp * ( h_center - enthalpy_min ) ) / grgr
    r_e_new = sqrt( r_e_new_sq )
    dif = abs(r_e_old - r_e_new) / r_e_new

    if (r_e_new / r_e_old > 2 .or. isnan(r_e_new) ) then 
      write(*,*) "r_e_old :", r_e_old
      write(*,*) "r_e_new :", r_e_new
      write(*,*) "grgr    :", grgr
      write(*,*) "r_e_new_sq :", r_e_new_sq
      write(*,*) "sphi_m :", sphi_m
      stop 'r_e cannot be found.'
    endif
  end subroutine update_equatorial_radius

  subroutine update_angular_velocity(r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, &
                                     sphi_pole_h, sphi_equator_h, ww_equator_h)
    use rotation_law_mod
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h
    real(wp), intent(in) :: sphi_pole_h, sphi_equator_h, ww_equator_h
    real(wp) :: metric_diff, term_in_Omega_h
    real(wp), parameter :: TOLERANCE = 1.0e-3_wp
    integer :: s, m
    if (abs(r_ratio - 1.0e0_wp) < TOLERANCE) then
      Omega_c = 0.0e0_wp; Omega_e = 0.0e0_wp; Omg = 0.0e0_wp
      return
    end if
    select case(trim(solver_type))
    case("uniform")
      call uniform_rotation()
    case("const_j")
      call const_j_rotation()
    case("uryu")
      call uryu_rotation()
    case default
      stop "Unknown solver type"
    end select
  contains
    subroutine uniform_rotation()
      metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                  + B_coup / 2.e0_wp * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.e0_wp - exp( r_e_new**2 * metric_diff )
      if (term_in_Omega_h >= 0.e0_wp) then
          Omega_c = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else
          write(*,"('Solving for axis ratio: ', f12.5)") r_ratio
          write(*,"(10A15)") "gama_pole", "rho_pole", "gama_equator", "rho_equator", "sphi_pole", "sphi_equator"
          write(*,"(10es15.3)") gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, sphi_pole_h, sphi_equator_h
          stop "Omega can't be found; Line 99 of spin helper"
      end if
      Omg = Omega_c
      Omega_e = Omega_c
    end subroutine uniform_rotation
    subroutine const_j_rotation()
      real(wp) :: guess, rsm, wwsm, mum, sgp
      real(wp), parameter :: tolerance = 1.e-5_wp
      metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                  + B_coup / 2.e0_wp * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.e0_wp - exp( r_e_new**2 * metric_diff )
      if (term_in_Omega_h >= 0.e0_wp) then
        Omega_e = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else 
        stop "L172 in const_j_rotation"
      endif
      guess = Omega_e * 0.8e0_wp
      call find_omege_e(guess, r_e_new, rho_equator_h, gama_equator_h, &
                       ww_equator_h, rho_pole_h, gama_pole_h, tolerance, Omega_e, diff_rotation_const_j)

      term_in_Omega_h = abs(Omega_e - ww_equator_h) * exp(-2.e0_wp * r_e_new**2 * rho_equator_h)
      Omega_c = Omega_e + term_in_Omega_h / (1.e0_wp - term_in_Omega_h * abs(Omega_e - ww_equator_h)) / A_diff**2
      
      Omg(1,:) = Omega_c
      Omg(1:2*SDIV/3,MDIV) = Omega_c
      do s = 2, 2*SDIV/3
        do m = 1, MDIV-1
          rsm = rho(s,m)
          wwsm = ww(s,m)
          mum = mu(m)
          sgp = s_gp(s)
          call zbrent_rot(Omg(s-1,m) * 8.e-1_wp, r_e_new, rsm, wwsm, sgp, mum, 1.e-5_wp, omg(s,m), rotation_law_const_j)
          F_j(s,m) = (omg(s,m) - wwsm) * sgp**2 * (1.e0_wp - mum**2) &
                / ((1.e0_wp - sgp)**2 * exp(2.e0_wp * r_e_new**2 * rsm) - (omg(s,m) - wwsm)**2 * sgp**2 * (1.e0_wp - mum**2))
        end do
      end do
    end subroutine const_j_rotation
    subroutine uryu_rotation()
      real(wp) :: diff_Fmax, guess, Fa, rsm, wwsm, sgp, mum, omg_max_h
      real(wp), dimension(SDIV) :: omg_mu_0
      integer :: imax
      real(wp), parameter :: tolerance = 1.e-5_wp
      metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                  + B_coup / 2.e0_wp * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.e0_wp - exp( r_e_new**2 * metric_diff )
      if (term_in_Omega_h >= 0.e0_wp) then
        Omega_e = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else
        stop "L205 in uryu"
      endif
      Fmax_h = 2.e-2_wp

      diff_Fmax = 1.e0_wp
      Fmax_h    = Fmax_h / 2.e0_wp
      do while( abs(diff_Fmax) > 1.e-7_wp)
        guess = Omega_e
        call find_omege_e(guess, r_e_new,rho_equator_h,gama_equator_h,ww_equator_h, &
                        rho_pole_h,gama_pole_h, tolerance, Fa, diff_rotation_uryu) ! Fmax_h used here
        Omega_e = fa
        F_equator_h  = (Omega_e - ww_equator_h) / ( exp(2.e0_wp*r_e_new**2*rho_equator_h) - (Omega_e-ww_equator_h)**2 )
        if ( F_equator_h < 0.e0_wp ) stop "negative F_equator_h; L120 in uryu"

        Omega_c = Omega_e / lambda2
        omg_mu_0(1) = Omega_c
        mum = 0.e0_wp
        do s = 2, SDIV*2/3
            rsm = rho(s,1) ! hat
            wwsm= ww (s,1) ! hat
            sgp = s_gp(s)
            guess  = omg_mu_0(s-1)
            call zbrent_rot( guess, r_e_new, rsm, wwsm, sgp, mum, 1.e-5_wp, omg_mu_0(s), rotation_law_uryu)
        enddo
        imax      = maxloc( omg_mu_0, 1 )
        omg_max_h = omg_mu_0(imax)
        diff_Fmax = ( lambda1 - omg_max_h / Omega_c )
        !write(*,*) "Shooting Fmax:", Fmax_h, diff_Fmax; stop 4
        Fmax_h    = Fmax_h - diff_Fmax * 1.e-2_wp
      enddo
      !write(*,*) "Shooting Fmax:", Fmax_h, diff_Fmax, Omega_e, Omega_c
      
      Omg(1,:) = Omega_c
      Omg(1:3*SDIV/4,MDIV) = Omega_c
      do s = 2, SDIV*3/4
        do m = 1, MDIV-1
          rsm = rho(s,m) ! hat
          wwsm= ww (s,m) ! hat
          mum = mu(m)
          sgp = s_gp(s)
          guess  = Omg(s-1,m)
          call zbrent_rot( guess, r_e_new, rsm, wwsm, sgp, mum, 1.e-5_wp, omg(s,m), rotation_law_uryu)
          F_j(s,m) = (omg(s,m) - wwsm) * sgp**2 * (1.e0_wp - mum**2) &
                / ((1.e0_wp - sgp)**2 * exp(2.e0_wp * r_e_new**2 * rsm) - (omg(s,m) - wwsm)**2 * sgp**2 * (1.e0_wp - mum**2))
        enddo
      enddo
    end subroutine uryu_rotation
  end subroutine update_angular_velocity

  subroutine update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h)
    use rotation_law_mod, only: intF
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(wp), dimension(SDIV,MDIV) :: sgp_term, sin_theta_2d, sgp_2d
    logical, dimension(SDIV,MDIV) :: valid
    integer :: s, m

    sgp_term = spread(s_gp / (1.e0_wp - s_gp), 2, MDIV)
    sin_theta_2d = spread(sin_theta, 1, SDIV)
    velocity_sq = merge(0.e0_wp, ((Omg - ww) * sgp_term * sin_theta_2d * exp(-rho * r_e_new**2))**2, r_ratio == 1.e0_wp)

    where (velocity_sq > 1.e0_wp) velocity_sq = 0.e0_wp

    enthalpy = enthalpy_min + 0.5e0_wp * ( &
          r_e_new**2 * ( gama_pole_h + rho_pole_h - gama - rho + ( sphi**2 - sphi_pole_h**2 ) * B_coup / 2.e0_wp ) &
          - log( max(1.e-300_wp, 1.e0_wp-velocity_sq) )  )

    if ( trim(solver_type) == "const_j" ) then
      enthalpy = enthalpy + 0.5e0_wp * A_diff**2 * (Omg - Omega_c)**2
    elseif ( trim(solver_type) == "uryu" ) then
      do m = 1, MDIV
        do s = 1, SDIV
          enthalpy(s,m) = enthalpy(s,m) - intF(omg(s,m), F_j(s,m))
        end do
      end do
    endif

    sgp_2d = spread(s_gp, 2, MDIV)
    valid = (enthalpy > enthalpy_min) .and. (sgp_2d <= s_e)

    where (.not. valid)
      enthalpy = enthalpy_min
      pressure = 0.e0_wp
      energy   = 0.e0_wp
    elsewhere
      pressure = exp( interp_log_h_to_p( log(enthalpy) ) )
      energy   = exp( interp_log_p_to_e( log(pressure) ) )
    end where
    ! Rescale back metric potentials (except omega)
    rho   = rho   * r_e_new**2
    gama  = gama  * r_e_new**2
    alpha = alpha * r_e_new**2
    sphi  = sphi  * r_e_new
  end subroutine update_eos_and_velocity

  subroutine precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, wfac_cache, besseli_cache, besselk_cache, &
                                                dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                                e2alpha_r2_cache, Acoup_cache, Acoup4_cache)
    real(wp), intent(in)  :: r_e_new
    real(wp), intent(in)  :: root_mphi_re
    real(wp), intent(out) :: mr_cache(:), wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    real(wp), intent(out) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(out) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(out) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(wp), intent(out) :: e2alpha_r2_cache(:,:), Acoup_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m, n
    real(wp) :: s1(SDIV), m1(MDIV)

    mr_cache(:)    = root_mphi_re * s_gp(:) / (1.e0_wp - s_gp(:))
    wfac_cache(:)  = 1.e0_wp / (1.e0_wp - s_gp(:))**2

    do n = 0, LMAX
      do s = 1, SDIV
        besseli_cache(n+1,s) = besseli(2*n, mr_cache(s))
        besselk_cache(n+1,s) = besselk(2*n, mr_cache(s))
      end do
    end do

    dg_s_cache = deriv_s_vec(gama)
    dg_m_cache = deriv_m_vec(gama)
    dr_s_cache = deriv_s_vec(rho)
    dr_m_cache = deriv_m_vec(rho)
    dww_s_cache = deriv_s_vec(ww)
    dww_m_cache = deriv_m_vec(ww)
    ds_s_cache = deriv_s_vec(sphi)
    ds_m_cache = deriv_m_vec(sphi)

    s1 = s_gp * (1.e0_wp - s_gp)
    m1 = 1.e0_wp - mu**2
    d2g_ss_cache = spread(s1, 2, MDIV) * deriv_s_vec(dg_s_cache) + spread(1.e0_wp - 2.e0_wp * s_gp, 2, MDIV) * dg_s_cache
    d2g_mm_cache = spread(m1, 1, SDIV) * deriv_m_vec(dg_m_cache) - 2.e0_wp * spread(mu, 1, SDIV) * dg_m_cache
    e_gsm_cache      = exp(0.5e0_wp * gama)
    e_rsm_cache      = exp(-rho)
    e2alpha_r2_cache = exp(2.e0_wp * alpha) * r_e_new**2
    Acoup_cache      = exp(-sphi**2 * B_coup / 4.e0_wp)
    Acoup4_cache     = Acoup_cache**4
  end subroutine precompute_derivatives_and_bessels

  subroutine build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, &
                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                e2alpha_r2_cache, Acoup4_cache)
    real(wp), intent(in)  :: r_e_new
    real(wp), intent(out) :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(wp), intent(in)  :: dr_s_cache(:,:), dr_m_cache(:,:), dg_s_cache(:,:), dg_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(in)  :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(in)  :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(wp), intent(in)  :: e2alpha_r2_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m
    real(wp), dimension(SDIV,MDIV) :: esm, psm, Vphi, scal_p
    real(wp), dimension(SDIV,MDIV) :: s1, s2, m1, mum

    s1 = spread(s_gp * (1.e0_wp - s_gp), 2, MDIV)
    s2 = spread((s_gp / (1.e0_wp - s_gp))**2, 2, MDIV)
    mum = spread(mu, 1, SDIV)
    m1 = 1.e0_wp - mum**2

    esm  = energy   * Acoup4_cache
    psm  = pressure * Acoup4_cache
    Vphi = sphi**2 * mphi_r * 0.5e0_wp * e2alpha_r2_cache
    scal_p = sphi * e_gsm_cache

    S_metric_rho = e_gsm_cache * ( &
        16.e0_wp * pi * e2alpha_r2_cache / 2.e0_wp * (esm + psm) * s2 * (1.e0_wp+velocity_sq) / (1.e0_wp-velocity_sq) &
      + s2 * m1 * e_rsm_cache**2 * ( (s1*dww_s_cache)**2 + m1*dww_m_cache**2 ) &
      + s1 * dg_s_cache - mum * dg_m_cache &
      + rho / 2.e0_wp * ( ( 16.e0_wp * pi * e2alpha_r2_cache * psm - 4.e0_wp * Vphi ) * s2 &
        - s1 * dg_s_cache * ( s1 / 2.e0_wp * dg_s_cache + 1.e0_wp ) &
        - dg_m_cache * ( m1 / 2.e0_wp * dg_m_cache - mum ) ) )

    S_metric_gama = e_gsm_cache * ( ( 16.e0_wp * pi * e2alpha_r2_cache * psm - 4.e0_wp * Vphi ) * s2 &
      + gama / 2.e0_wp * ( ( 16.e0_wp * pi * e2alpha_r2_cache * psm - 4.e0_wp * Vphi ) * s2 &
        - ( s1 * dg_s_cache )**2 / 2.e0_wp - m1 * dg_m_cache**2 / 2.e0_wp ) )

    S_metric_omega = e_gsm_cache * e_rsm_cache * ( &
      - 16.e0_wp * pi * e2alpha_r2_cache * ( Omg - ww ) * ( esm + psm ) * s2 / (1.e0_wp-velocity_sq) &
      + ww * ( -0.5e0_wp * 16.e0_wp * pi * e2alpha_r2_cache * s2 * &
        ( ( ( 1.e0_wp + velocity_sq ) * esm + 2.e0_wp * velocity_sq*psm ) / ( 1.e0_wp - velocity_sq ) ) &
        - s1 * ( 2.e0_wp * dr_s_cache + 0.5e0_wp * dg_s_cache ) &
        + mum * ( 2.e0_wp * dr_m_cache + 0.5e0_wp * dg_m_cache) &
        + 0.25e0_wp * s1**2 * ( 4.e0_wp * dr_s_cache**2 - dg_s_cache**2) &
        + 0.25e0_wp * m1 * ( 4.e0_wp * dr_m_cache**2 - dg_m_cache**2 ) &
        - m1 * e_rsm_cache**2 * ( spread(s_gp,2,MDIV)**4 * dww_s_cache**2 + s2 * m1 * dww_m_cache**2) &
        - 2.e0_wp * Vphi * s2 ) )

    if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp ) then
      S_metric_sphi = - r_e_new**2 * s2 * scal_p * mphi_r &
        + s2 * scal_p * e2alpha_r2_cache * ( -2.e0_wp * pi * B_coup * ( esm - 3.e0_wp * psm ) + mphi_r ) &
        + scal_p * ( s1 * spread(1.e0_wp-s_gp,2,MDIV) * dg_s_cache       &
        + s1 * s1 * ( d2g_ss_cache * 0.5e0_wp + dg_s_cache**2 * 0.25e0_wp ) &
        + m1      * ( d2g_mm_cache * 0.5e0_wp + dg_m_cache**2 * 0.25e0_wp ) &
        - mum * dg_m_cache )
      !> v1
      !S_metric_sphi = - r_e_new**2 * s2 * scal_p * mphi_r &
      !    + r_e_new**2 * s2 * scal_p * exp(2.e0_wp*alpha) * ( -2.e0_wp * pi * B_coup * (esm-3.e0_wp*psm) + mphi_r) &
      !    + scal_p * ( - 2.e0_wp * s1 * spread(s_gp,2,MDIV) * dg_s_cache &
      !    + s1**2 * ( d2g_ss_cache + dg_s_cache**2 / 4.e0_wp ) + s1 * dg_s_cache &
      !    + m1 * d2g_mm_cache / 2.e0_wp + dg_m_cache**2 * r_e_new**2 * s2 / 4.e0_wp &
      !    - 2.e0_wp * mum * dg_m_cache )
    else
      S_metric_sphi = -s1**2 * dg_s_cache * ds_s_cache - m1 * dg_m_cache * ds_m_cache &
        + sphi * ( -2.e0_wp * pi * B_coup * (esm - 3.e0_wp*psm) ) * s2 * e2alpha_r2_cache
    endif
  end subroutine build_source_terms

  subroutine angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                 D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    real(wp), intent(in)  :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(wp), intent(out) :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    integer :: n, k
    real(wp) :: sum_rho, sum_sphi, sum_gama, sum_omega
    real(wp), dimension(MDIV) :: Int_m

    if (allocated(angular_quad_weights)) then
      if (size(angular_quad_weights) /= MDIV) deallocate(angular_quad_weights)
    end if
    if (.not. allocated(angular_quad_weights)) allocate(angular_quad_weights(MDIV))
    call compute_effective_d01gaf_weights(mu, angular_quad_weights)

    n = 0
    do k = 1, SDIV
      Int_m(:) = P_2n(:,n+1) * S_metric_rho(k,:)
      sum_rho = dot_product(angular_quad_weights, Int_m)

      Int_m(:) = P_2n(:,n+1) * S_metric_sphi(k,:)
      sum_sphi = dot_product(angular_quad_weights, Int_m)

      D1_metric_rho  (n+1,k) = sum_rho
      D1_metric_sphi (n+1,k) = sum_sphi
      D1_metric_gama (n+1,k) = 0.e0_wp
      D1_metric_omega(n+1,k) = 0.e0_wp
    enddo

    do n = 1, LMAX
      do k = 1, SDIV
        Int_m(:) = P_2n(:,n+1) * S_metric_rho(k,:)
        sum_rho = dot_product(angular_quad_weights, Int_m)

        Int_m(:) = P_2n(:,n+1) * S_metric_sphi(k,:)
        sum_sphi = dot_product(angular_quad_weights, Int_m)

        Int_m(:) = sin_2n_1_theta(:,n) * S_metric_gama(k,:)
        sum_gama = dot_product(angular_quad_weights, Int_m)

        Int_m(:) = sin_theta(:) * P1_2n_1(:,n+1) * S_metric_omega(k,:)
        sum_omega = dot_product(angular_quad_weights, Int_m)

        D1_metric_rho  (n+1,k) = sum_rho
        D1_metric_sphi (n+1,k) = sum_sphi
        D1_metric_gama (n+1,k) = sum_gama
        D1_metric_omega(n+1,k) = sum_omega
      enddo
    enddo
  end subroutine angular_integration

  subroutine radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
                                D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, &
                                root_mphi_re, wfac_cache, besseli_cache, besselk_cache)
    real(wp), intent(in)  :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    real(wp), intent(out) :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
    real(wp), intent(in)  :: root_mphi_re
    real(wp), intent(in)  :: wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    integer, parameter :: radial_tile = 32
    integer :: n
    real(wp) :: weighted_source(SDIV)

    if (allocated(radial_quad_weights)) then
      if (size(radial_quad_weights) /= SDIV) deallocate(radial_quad_weights)
    end if
    if (.not. allocated(radial_quad_weights)) allocate(radial_quad_weights(SDIV))
    call compute_effective_d01gaf_weights(s_gp, radial_quad_weights)

    D2_metric_gama (:,1) = 0.e0_wp
    D2_metric_omega(:,1) = 0.e0_wp

    do n = 0, LMAX
      weighted_source = radial_quad_weights * D1_metric_rho(n+1,:)
      call integrate_kernel_tiled(f_rho, n+1, weighted_source, D2_metric_rho(:,n+1))

      if (mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp) then
        weighted_source = radial_quad_weights * wfac_cache * D1_metric_sphi(n+1,:) * root_mphi_re
        call integrate_massive_sphi(n+1, weighted_source, D2_metric_sphi(:,n+1))
      else
        weighted_source = radial_quad_weights * D1_metric_sphi(n+1,:)
        call integrate_kernel_tiled(f_rho, n+1, weighted_source, D2_metric_sphi(:,n+1))
      end if
    end do

    do n = 1, LMAX
      weighted_source = radial_quad_weights * D1_metric_gama(n+1,:)
      call integrate_kernel_tiled(f_gama, n+1, weighted_source, D2_metric_gama(:,n+1))

      weighted_source = radial_quad_weights * D1_metric_omega(n+1,:)
      call integrate_omega_tiled(n+1, weighted_source, D2_metric_omega(:,n+1))
    end do

  contains
    subroutine integrate_kernel_tiled(kernel, n_idx, source_weights, out_values)
      real(wp), intent(in) :: kernel(SDIV,LMAX+1,SDIV)
      integer, intent(in) :: n_idx
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: j, s_lo, s_hi, tile_len
      real(wp) :: acc(radial_tile)

      do s_lo = 1, SDIV, radial_tile
        s_hi = min(s_lo + radial_tile - 1, SDIV)
        tile_len = s_hi - s_lo + 1
        acc(1:tile_len) = 0.e0_wp
        do j = 1, SDIV
          acc(1:tile_len) = acc(1:tile_len) + kernel(s_lo:s_hi, n_idx, j) * source_weights(j)
        end do
        out_values(s_lo:s_hi) = acc(1:tile_len)
      end do
    end subroutine integrate_kernel_tiled

    subroutine integrate_omega_tiled(n_idx, source_weights, out_values)
      integer, intent(in) :: n_idx
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: j, s_lo, s_hi, split_hi, split_lo, tile_len
      real(wp) :: acc(radial_tile)

      do s_lo = 1, SDIV, radial_tile
        s_hi = min(s_lo + radial_tile - 1, SDIV)
        tile_len = s_hi - s_lo + 1
        acc(1:tile_len) = 0.e0_wp
        do j = 1, SDIV
          split_hi = min(s_hi, j)
          if (split_hi >= s_lo) then
            acc(1:split_hi-s_lo+1) = acc(1:split_hi-s_lo+1) + &
              f_gama(s_lo:split_hi, n_idx, j) * source_weights(j)
          end if

          split_lo = max(s_lo, j + 1)
          if (split_lo <= s_hi) then
            acc(split_lo-s_lo+1:tile_len) = acc(split_lo-s_lo+1:tile_len) + &
              f_rho(split_lo:s_hi, n_idx, j) * source_weights(j)
          end if
        end do
        out_values(s_lo:s_hi) = acc(1:tile_len)
      end do
    end subroutine integrate_omega_tiled

    subroutine integrate_massive_sphi(n_idx, source_weights, out_values)
      integer, intent(in) :: n_idx
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      real(wp) :: left_prefix(SDIV), right_suffix(SDIV)
      real(wp) :: source_left(SDIV), source_right(SDIV)
      integer :: s

      source_left = source_weights * besseli_cache(n_idx,:)
      source_right = source_weights * besselk_cache(n_idx,:)

      left_prefix(1) = source_left(1)
      do s = 2, SDIV
        left_prefix(s) = left_prefix(s-1) + source_left(s)
      end do

      right_suffix(SDIV) = source_right(SDIV)
      do s = SDIV-1, 1, -1
        right_suffix(s) = right_suffix(s+1) + source_right(s)
      end do

      out_values(1) = besseli_cache(n_idx,1) * right_suffix(1)
      do s = 2, SDIV
        out_values(s) = besselk_cache(n_idx,s) * left_prefix(s-1) + &
                        besseli_cache(n_idx,s) * right_suffix(s)
      end do
    end subroutine integrate_massive_sphi
  end subroutine radial_integration

  subroutine compute_effective_d01gaf_weights(x, w)
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: w(:)
    real(wp), allocatable :: basis(:)
    real(wp) :: ans, er
    integer :: i, ifail, n

    n = size(x)
    if (size(w) /= n) stop "compute_effective_d01gaf_weights: size mismatch"

    allocate(basis(n))
    basis = 0.e0_wp

    do i = 1, n
      basis(i) = 1.e0_wp
      call d01gaf(x, basis, n, ans, er, ifail)
      if (ifail /= 0) stop "compute_effective_d01gaf_weights: d01gaf failed"
      w(i) = ans
      basis(i) = 0.e0_wp
    end do

    deallocate(basis)
  end subroutine compute_effective_d01gaf_weights

  subroutine sum_coefficients_and_get_targets(target_rho, target_gama, target_ww, target_sphi, &
                                              D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    real(wp), intent(out) :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
    real(wp), intent(in)  :: D2_metric_rho(SDIV,LMAX+1), D2_metric_gama(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1), D2_metric_sphi(SDIV,LMAX+1)
    integer :: n
    real(wp), dimension(SDIV,MDIV) :: exp_mhalf_gsm, exp_rsm_mhalf_gsm
    real(wp), dimension(SDIV,MDIV) :: sum_rho, sum_sphi, sum_gama, sum_omega
    real(wp), dimension(MDIV) :: sin_theta_inv

    exp_mhalf_gsm = exp(-0.5e0_wp * gama)
    exp_rsm_mhalf_gsm = exp(rho - 0.5e0_wp * gama)

    ! Monopole term (n=0)
    sum_rho  = -exp_mhalf_gsm * matmul(reshape(D2_metric_rho(:,1), [SDIV, 1]), reshape(P_2n(:,1), [1, MDIV]))
    sum_gama = 0.e0_wp
    sum_omega = 0.e0_wp

    ! Precompute 1/sin(theta) to avoid division in the loop
    sin_theta_inv = 0.e0_wp; where (sin_theta > 1.e-12_wp) sin_theta_inv = 1.e0_wp / sin_theta

    ! Higher multipoles (n=1 to LMAX)
    do n = 1, LMAX
      sum_rho  = sum_rho  - exp_mhalf_gsm * matmul(reshape(D2_metric_rho(:,n+1), [SDIV, 1]), reshape(P_2n(:,n+1), [1, MDIV]))
      
      ! Vectorized update for sum_gama and sum_omega
      sum_gama(:,1:MDIV-1) = sum_gama(:,1:MDIV-1) - (2.e0_wp/pi) * exp_mhalf_gsm(:,1:MDIV-1) * &
            (matmul(reshape(D2_metric_gama(:,n+1), [SDIV, 1]), reshape(sin_2n_1_theta(1:MDIV-1,n) * sin_theta_inv(1:MDIV-1), [1, MDIV-1])) / (2.e0_wp*n-1.e0_wp))
      sum_gama(:,MDIV) = sum_gama(:,MDIV) - (2.e0_wp/pi) * exp_mhalf_gsm(:,MDIV) * D2_metric_gama(:,n+1)

      sum_omega(:,1:MDIV-1) = sum_omega(:,1:MDIV-1) - exp_rsm_mhalf_gsm(:,1:MDIV-1) * &
            (matmul(reshape(D2_metric_omega(:,n+1), [SDIV, 1]), reshape(P1_2n_1(1:MDIV-1,n+1) * sin_theta_inv(1:MDIV-1), [1, MDIV-1])) / (2.e0_wp*n*(2.e0_wp*n-1.e0_wp)))
      sum_omega(:,MDIV) = sum_omega(:,MDIV) + exp_rsm_mhalf_gsm(:,MDIV) * D2_metric_omega(:,n+1) / 2.e0_wp
    end do

    if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp ) then
      sum_sphi = -exp_mhalf_gsm * matmul(reshape(D2_metric_sphi(:,1), [SDIV, 1]), reshape(P_2n(:,1), [1, MDIV]))
      do n = 1, LMAX
        sum_sphi = sum_sphi - exp_mhalf_gsm * (2.e0_wp*dble(n)+1.e0_wp) &
                 * matmul(reshape(D2_metric_sphi(:,n+1), [SDIV, 1]), reshape(P_2n(:,n+1), [1, MDIV]))
      enddo
    else
      sum_sphi = - matmul(reshape(D2_metric_sphi(:,1), [SDIV, 1]), reshape(P_2n(:,1), [1, MDIV]))
      do n = 1, LMAX
        sum_sphi = sum_sphi - matmul(reshape(D2_metric_sphi(:,n+1), [SDIV, 1]), reshape(P_2n(:,n+1), [1, MDIV]))
      enddo
    endif

    target_rho  = sum_rho
    target_gama = sum_gama
    target_ww   = sum_omega
    target_sphi = sum_sphi
  end subroutine sum_coefficients_and_get_targets

  subroutine update_alpha_potential(r_e_new, dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                    ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_rsm_cache)
    real(wp), intent(in) :: r_e_new
    real(wp), intent(in) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(in) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(in) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_rsm_cache(:,:)
    integer :: s, m
    real(wp) :: sgp, s1, sgp_ratio
    real(wp), dimension(SDIV,MDIV) :: da_dm, d_gama_sm_all
    real(wp), dimension(MDIV) :: m1_vec, d_gama_s_row, d_gama_m_row, d_rho_s_row, d_rho_m_row
    real(wp), dimension(MDIV) :: d_sphi_s_row, d_sphi_m_row, d_gama_sm_row, d_ww_s_row, d_ww_m_row, d_gama_ss_row, d_gama_mm_row
    real(wp), dimension(MDIV) :: temp1_row, temp2_row, temp3_row, temp4_row, temp5_row, temp6_row, temp7_row, temp8_row, temp9_row
    real(wp), dimension(SDIV) :: s1_cache, sgp_ratio_cache
    real(wp) :: adj_const(SDIV)

    alpha(:,:) = 0.e0_wp
    if (r_ratio == 1.e0_wp) then
      da_dm(:,:) = 0.0e0_wp
    else
      s1_cache = s_gp * (1.e0_wp-s_gp)
      sgp_ratio_cache = s_gp / (1.e0_wp-s_gp)

      da_dm(1,:) = 0.0e0_wp
      d_gama_sm_all = deriv_sm_vec(gama)
      m1_vec = 1.e0_wp - mu**2
      do s = 2, SDIV
        sgp = s_gp(s)
        s1  = s1_cache(s)
        sgp_ratio = sgp_ratio_cache(s)
        d_gama_s_row  = dg_s_cache(s,:)
        d_gama_m_row  = dg_m_cache(s,:)
        d_rho_s_row   = dr_s_cache(s,:)
        d_rho_m_row   = dr_m_cache(s,:)
        d_sphi_s_row  = ds_s_cache(s,:)
        d_sphi_m_row  = ds_m_cache(s,:)
        d_gama_sm_row = d_gama_sm_all(s,:)
        d_ww_s_row    = dww_s_cache(s,:)
        d_ww_m_row    = dww_m_cache(s,:)
        d_gama_ss_row = d2g_ss_cache(s,:)
        d_gama_mm_row = d2g_mm_cache(s,:)

        temp1_row = 2.e0_wp * sgp**2 * sgp_ratio * m1_vec * d_ww_s_row * d_ww_m_row * (1.e0_wp + s1 * d_gama_s_row) &
          - ( (sgp**2 * d_ww_s_row)**2 - (sgp * d_ww_m_row * sgp_ratio)**2 * m1_vec ) * (-mu + m1_vec * d_gama_m_row)
        temp2_row = 1.e0_wp / ( m1_vec * (1.e0_wp + s1 * d_gama_s_row)**2 + (-mu + m1_vec * d_gama_m_row)**2 )
        temp3_row = s1 * d_gama_ss_row + (s1 * d_gama_s_row)**2
        temp4_row = d_gama_m_row * (-mu + m1_vec * d_gama_m_row)
        temp5_row = ( (s1 * (d_rho_s_row + d_gama_s_row))**2 - m1_vec * (d_rho_m_row + d_gama_m_row)**2 ) * (-mu + m1_vec * d_gama_m_row)
        temp6_row = s1 * m1_vec * (  (d_rho_s_row + d_gama_s_row) * (d_rho_m_row + d_gama_m_row) / 2.e0_wp + d_gama_sm_row + d_gama_s_row * d_gama_m_row  ) * (1.e0_wp + s1 * d_gama_s_row)
        temp7_row = s1 * mu * d_gama_s_row * ( 1.e0_wp + s1 * d_gama_s_row )
        temp8_row = m1_vec * (e_rsm_cache(s,:)**2)
        temp9_row = -temp2_row * (-mu + m1_vec * d_gama_m_row) * ( (s1 * d_sphi_s_row)**2 - m1_vec * d_sphi_m_row**2 ) &
              - m1_vec * s1 * ( 1.e0_wp + s1 * d_gama_s_row ) * 2.e0_wp * d_sphi_m_row * d_sphi_s_row

        da_dm(s,:) = - (d_rho_m_row + d_gama_m_row) / 2.e0_wp &
          - temp2_row * ( (temp3_row - d_gama_mm_row - temp4_row) * (-mu + m1_vec * d_gama_m_row) / 2.e0_wp &
          + temp5_row / 4.e0_wp - temp6_row  + temp7_row + temp8_row * temp1_row / 4.e0_wp ) + temp9_row
      end do
    end if

    alpha(:, 1) = 0.0e0_wp

    do m = 1, MDIV-1
        alpha(:,m+1) = alpha(:,m) + dm * ( da_dm(:,m+1) + da_dm(:,m) ) * 0.5e0_wp
    enddo
    
    alpha(SDIV,:) = 0.e0_wp
    adj_const = alpha(:,MDIV) - ( gama(:,MDIV) - rho(:,MDIV) )/2.e0_wp
    alpha = alpha - spread(adj_const, DIM=2, NCOPIES=MDIV)
    if (any(alpha .ge. 300.0)) then
      write(*,*) "Error: Alpha fails in at least one row."
      stop "alpha fails"
    end if
    omg= omg / r_e_new
    ww = ww / r_e_new
  end subroutine update_alpha_potential

  subroutine get_all_targets(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, root_mphi_re, &
                             out_target_rho, out_target_gama, out_target_ww, out_target_sphi)
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(wp), intent(in) :: root_mphi_re
    real(wp), intent(out) :: out_target_rho(SDIV,MDIV), out_target_gama(SDIV,MDIV), out_target_ww(SDIV,MDIV), out_target_sphi(SDIV,MDIV)
    real(wp) :: t0, t1, dt_precompute, dt_build, dt_angular, dt_radial, dt_sum
    integer, save :: target_call_count = 0

    if (timing) then
      target_call_count = target_call_count + 1
      dt_precompute = 0.e0_wp
      dt_build = 0.e0_wp
      dt_angular = 0.e0_wp
      dt_radial = 0.e0_wp
      dt_sum = 0.e0_wp
      call cpu_time(t0)
    end if
    call precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, wfac_cache, besseli_cache, besselk_cache, &
         dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
         d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup_cache, Acoup4_cache)
    if (timing) then
      call cpu_time(t1)
      dt_precompute = t1 - t0
      call cpu_time(t0)
    end if

    call build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
         dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
         d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup4_cache)
    if (timing) then
      call cpu_time(t1)
      dt_build = t1 - t0
      call cpu_time(t0)
    end if

    call angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
         D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    if (timing) then
      call cpu_time(t1)
      dt_angular = t1 - t0
      call cpu_time(t0)
    end if

    call radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
         D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, root_mphi_re, wfac_cache, &
         besseli_cache, besselk_cache)
    if (timing) then
      call cpu_time(t1)
      dt_radial = t1 - t0
      call cpu_time(t0)
    end if

    call sum_coefficients_and_get_targets(out_target_rho, out_target_gama, out_target_ww, out_target_sphi, &
         D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    if (timing) then
      call cpu_time(t1)
      dt_sum = t1 - t0

      write(*,'(A,I0,A,6(1X,ES12.5))') 'get_all_targets call ', target_call_count, ':', &
        dt_precompute, dt_build, dt_angular, dt_radial, dt_sum, &
        dt_precompute + dt_build + dt_angular + dt_radial + dt_sum
      stop 930
    end if
  end subroutine get_all_targets

  subroutine relaxation(r_e_new, target_rho, target_gama, target_ww, target_sphi, root_mphi_re, n_of_it)
    use hybrid_relaxation, only: hybrid_update
    use anderson_optimized, only: anderson_accel_optimized
    real(wp), intent(in) :: r_e_new
    real(wp), intent(in) :: target_rho(SDIV,MDIV), target_gama(SDIV,MDIV)
    real(wp), intent(in) :: target_ww(SDIV,MDIV),  target_sphi(SDIV,MDIV)
    real(wp), intent(in) :: root_mphi_re
    integer, intent(in) :: n_of_it
    integer :: s, m
    integer, parameter :: m_hist = 3
    real(wp), allocatable, save :: hist_f_rho(:,:,:)
    real(wp), allocatable, save :: hist_f_gama(:,:,:)
    real(wp), allocatable, save :: hist_f_ww(:,:,:)
    real(wp), allocatable, save :: hist_f_sphi(:,:,:)
    real(wp), allocatable, save :: hist_x_rho(:,:,:)
    real(wp), allocatable, save :: hist_x_gama(:,:,:)
    real(wp), allocatable, save :: hist_x_ww(:,:,:)
    real(wp), allocatable, save :: hist_x_sphi(:,:,:)

    select case (trim(relaxation_scheme))
    case ('anderson')
        if (.not. allocated(hist_f_rho)) then
          allocate(hist_f_rho(SDIV,MDIV,m_hist))
          allocate(hist_f_gama(SDIV,MDIV,m_hist))
          allocate(hist_f_ww(SDIV,MDIV,m_hist))
          allocate(hist_f_sphi(SDIV,MDIV,m_hist))
          allocate(hist_x_rho(SDIV,MDIV,m_hist))
          allocate(hist_x_gama(SDIV,MDIV,m_hist))
          allocate(hist_x_ww(SDIV,MDIV,m_hist))
          allocate(hist_x_sphi(SDIV,MDIV,m_hist))
        end if
        call anderson_accel_optimized(rho,  target_rho,  hist_f_rho,  hist_x_rho,  n_of_it, m_hist, .true.)
        call anderson_accel_optimized(gama, target_gama, hist_f_gama, hist_x_gama, n_of_it, m_hist, .true.)
        call anderson_accel_optimized(ww,   target_ww,   hist_f_ww,   hist_x_ww,   n_of_it, m_hist, .true.)
        call anderson_accel_optimized(sphi, target_sphi, hist_f_sphi, hist_x_sphi, n_of_it, m_hist, .true.)
    case ('hybrid')
        if (.not. allocated(hist_f_rho)) then
          allocate(hist_f_rho(SDIV,MDIV,m_hist))
          allocate(hist_f_gama(SDIV,MDIV,m_hist))
          allocate(hist_f_ww(SDIV,MDIV,m_hist))
          allocate(hist_f_sphi(SDIV,MDIV,m_hist))
          allocate(hist_x_rho(SDIV,MDIV,m_hist))
          allocate(hist_x_gama(SDIV,MDIV,m_hist))
          allocate(hist_x_ww(SDIV,MDIV,m_hist))
          allocate(hist_x_sphi(SDIV,MDIV,m_hist))
        end if
        call hybrid_update(rho,  target_rho,  hist_f_rho,  hist_x_rho,  n_of_it, m_hist)
        call hybrid_update(gama, target_gama, hist_f_gama, hist_x_gama, n_of_it, m_hist)
        call hybrid_update(ww,   target_ww,   hist_f_ww,   hist_x_ww,   n_of_it, m_hist)
        call hybrid_update(sphi, target_sphi, hist_f_sphi, hist_x_sphi, n_of_it, m_hist)
    case default
        stop "Unknown relaxation scheme specified"
    end select

    where(sphi .ne. sphi) sphi = 0.e0_wp

    if ( mphi_r == 0.e0_wp ) return
    
    do m=1,MDIV
      do s=2,SDIV
      if (sphi(s,m) < 0.e0_wp) then
          sphi(s,m) = sphi(s-1,m) * exp(root_mphi_re * (s_gp(s-1)/(1.e0_wp-s_gp(s-1)) - s_gp(s)/(1.e0_wp-s_gp(s))))
        end if
      end do
    end do
  end subroutine relaxation

  subroutine allocate_workspace
    allocate(e_gsm_cache(SDIV,MDIV), e_rsm_cache(SDIV,MDIV), e2alpha_r2_cache(SDIV,MDIV))
    allocate(Acoup_cache(SDIV,MDIV), Acoup4_cache(SDIV,MDIV))
    allocate(dg_s_cache(SDIV,MDIV), dg_m_cache(SDIV,MDIV), d2g_ss_cache(SDIV,MDIV), d2g_mm_cache(SDIV,MDIV))
    allocate(dr_s_cache(SDIV,MDIV), dr_m_cache(SDIV,MDIV), dww_s_cache(SDIV,MDIV), dww_m_cache(SDIV,MDIV))
    allocate(ds_s_cache(SDIV,MDIV), ds_m_cache(SDIV,MDIV)) 
    allocate(mr_cache(SDIV), besseli_cache(LMAX+1,SDIV), besselk_cache(LMAX+1,SDIV), wfac_cache(SDIV))
    allocate(target_rho(SDIV,MDIV), target_gama(SDIV,MDIV), target_ww(SDIV,MDIV), target_sphi(SDIV,MDIV))
    allocate(S_metric_rho(SDIV,MDIV), S_metric_gama(SDIV,MDIV), S_metric_omega(SDIV,MDIV), S_metric_sphi(SDIV,MDIV))
    allocate(D1_metric_rho(LMAX+1,SDIV), D1_metric_gama(LMAX+1,SDIV), D1_metric_omega(LMAX+1,SDIV), D1_metric_sphi(LMAX+1,SDIV))
    allocate(D2_metric_rho(SDIV,LMAX+1), D2_metric_gama(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1), D2_metric_sphi(SDIV,LMAX+1))
  end subroutine allocate_workspace

  subroutine deallocate_workspace
    deallocate(e_gsm_cache, e_rsm_cache, e2alpha_r2_cache)
    deallocate(Acoup_cache, Acoup4_cache)
    deallocate(dg_s_cache, dg_m_cache, d2g_ss_cache, d2g_mm_cache)
    deallocate(dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache)
    deallocate(ds_s_cache, ds_m_cache) 
    deallocate(mr_cache, besseli_cache, besselk_cache, wfac_cache)
    deallocate(target_rho, target_gama, target_ww, target_sphi)
    deallocate(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi)
    deallocate(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    deallocate(D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    if (allocated(radial_quad_weights)) deallocate(radial_quad_weights)
    if (allocated(angular_quad_weights)) deallocate(angular_quad_weights)
  end subroutine deallocate_workspace

  subroutine output_helper(D2_metric_rho, D2_metric_omega)
    real(wp), intent(in) :: D2_metric_rho(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1)
    real(wp) :: r_inf, rho_0
    character(32) :: fil1, fil2, fil3, fil4, fil5, fil6, fil7, fil8
    character(64) :: tail_fmt
    integer :: s, m, unit, ios, sig_digits, field_width
    real(wp), external :: n0_at_e

    r_inf = r_e * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1.e0_wp - s_gp(SDIV - 1) ))**s_pwr
    M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2.e0_wp * r_inf**3 * ( C**2 / G / Mass )**3
    S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2.e0_wp * r_inf**5 * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2.e0_wp * r_inf**5 * ( C**2 / G / Mass )**5

    if (.not. output) return
    open(newunit=unit, file="Cont/moment_tail.dat", status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim("Cont/moment_tail.dat")
      return
    end if
    sig_digits = precision(1.0_wp) - 1
    field_width = sig_digits + 10
    write(tail_fmt,'("(30es",i0,".",i0,"e3)")') field_width, sig_digits
    do s = 1, SDIV-1
      write(unit,tail_fmt) s_gp(s), D2_metric_rho(s,:), D2_metric_omega(s,:)
    end do
    close(unit)

    write(*,"(A)") " ", " Off-loading data ..."
    
    write(fil1,"(f6.2)") ang_mom
    write(fil2,"(f16.3)") mass_0/MSUN
    write(fil3,"(es15.2)") B_coup
    write(fil4,"(es15.2)") sqrt(mphi_r*1.e10_wp/KAPPA)*l_uni
    rho_0 = n0_at_e( energy(1,1) ) * MB
    write(fil5,"(es15.3)") rho_0
    write(fil6,"(f15.3)") sphi_m
    select case(trim(solver_type))
    case("uniform")
      write(fil7,"(A)") ".dat"
    case("const_j")
      write(fil7,"(f5.2)") A_diff
      write(fil7,"(A)") "_A"//trim(adjustl(fil7))//".dat"
    case("uryu")
      write(fil7,"(f5.2)") lambda1
      write(fil8,"(f5.2)") lambda2
      write(fil7,"(A)") "_L1"//trim(adjustl(fil7))//"_L2"//trim(adjustl(fil8))//".dat"
    case default
      stop "Unknown solver type"
    end select

    call write_output_file("./Cont/"//trim(adjustl(eos_file))//&
      "_J"//trim(adjustl(fil1))//&
      "_Mb"//trim(adjustl(fil2))//"_B"//trim(adjustl(fil3))//&
      "_mphi"//trim(adjustl(fil4))//"_rhoc"//trim(adjustl(fil5))//"_sphim"//trim(adjustl(fil6))//trim(adjustl(fil7)) )

    call write_output_file("./Res/res.dat")
  end subroutine output_helper

  subroutine write_output_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios, s, m
    integer :: sig_digits, exp_digits, field_width
    real(wp) :: rho_0_val
    character(len=256) :: header_fmt, data_fmt
    real(wp), external :: n0_at_e

    ! Build format strings from precision/range so real128 has enough exponent digits.
    sig_digits = precision(1.0_wp) - 1
    exp_digits = 3
    if (range(1.0_wp) > 999) exp_digits = 4
    field_width = sig_digits + exp_digits + 10
    write(header_fmt,'("(3i5, 5es",i0,".",i0,"e",i0,")")') field_width, sig_digits, exp_digits
    write(data_fmt,  '("(13es",i0,".",i0,"e",i0,")")') field_width, sig_digits, exp_digits

    ! Use newunit to get a free file unit, preventing conflicts
    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim(filename)
      return
    end if

    ! Write header
    write(unit, fmt=header_fmt) SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.e5_wp, &
            energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)), Omega_c* (C/sqrt(kappa))

    ! Write data grid
    do s = 1, SDIV
      do m = 1, MDIV
        if (enthalpy(s,m) > enthalpy_min) then
          rho_0_val = n0_at_e(energy(s,m)) * MB
        else
          rho_0_val = 0.e0_wp
        end if

        write(unit, fmt=data_fmt) s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), &
          ww(s,m) * (C/sqrt(kappa)), pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), &
          enthalpy(s,m), rho_0_val, velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), &
          sphi(s,m) * sqrt(B_coup)
      end do
    end do

    close(unit)
  end subroutine write_output_file

end module spin_helper
