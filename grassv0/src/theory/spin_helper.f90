module spin_helper
  use para_mod
  use brent_mod, only : find_omege_e, zbrent_rot
  use toolkit_mod, only : deriv_sm, deriv_s, deriv_m, besseli, besselk, interp, interp_log_h_to_p, interp_log_p_to_e
  use simpson_mod
  use nag_compat_mod, only : d01gaf
  implicit none
  real(8) :: mphi_tran = 1.d-11
  real(8) :: dif
  interface
    subroutine DPOSV(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
      character, intent(in) :: UPLO
      integer, intent(in) :: N, NRHS, LDA, LDB
      integer, intent(out) :: INFO
      real(8), intent(inout) :: A(LDA,N), B(LDB,NRHS)
    end subroutine DPOSV

    subroutine DTRSV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
      character, intent(in) :: UPLO, TRANS, DIAG
      integer, intent(in) :: N, LDA, INCX
      real(8), intent(in) :: A(LDA,N)
      real(8), intent(inout) :: X(N)
    end subroutine DTRSV
  end interface

  ! --- Workspace arrays, moved to module scope ---
  real(8), allocatable, target :: e_gsm_cache(:,:), e_rsm_cache(:,:), e2alpha_r2_cache(:,:)
  real(8), allocatable, target :: Acoup_cache(:,:), Acoup4_cache(:,:)
  real(8), allocatable, target :: dg_s_cache(:,:), dg_m_cache(:,:), d2g_ss_cache(:,:), d2g_mm_cache(:,:)
  real(8), allocatable, target :: dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
  real(8), allocatable, target :: ds_s_cache(:,:), ds_m_cache(:,:) 
  real(8), allocatable, target :: mr_cache(:), besseli_cache(:,:), besselk_cache(:,:), wfac_cache(:)
  real(8), allocatable, target :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
  real(8), allocatable, target :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
  real(8), allocatable, target :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
  real(8), allocatable, target :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
  real(8), allocatable, target :: target_field_perturbed(:,:)

contains
  pure function deriv_s_vec(f) result(df_ds)
    use para_mod, only : SDIV, MDIV, DS
    real(8), dimension(SDIV,MDIV), intent(in) :: f
    real(8), dimension(SDIV,MDIV) :: df_ds
    integer :: m
    df_ds(1,:) = (f(2,:) - f(1,:)) / DS
    df_ds(SDIV,:) = (f(SDIV,:) - f(SDIV-1,:)) / DS
    do m = 1, MDIV
      df_ds(2:SDIV-1,m) = (f(3:SDIV,m) - f(1:SDIV-2,m)) / (2.d0 * DS)
    end do
  end function deriv_s_vec

  pure function deriv_m_vec(f) result(df_dm)
    use para_mod, only : SDIV, MDIV, DM
    real(8), dimension(SDIV,MDIV), intent(in) :: f
    real(8), dimension(SDIV,MDIV) :: df_dm
    integer :: s
    df_dm(:,1) = (f(:,2) - f(:,1)) / DM
    df_dm(:,MDIV) = (f(:,MDIV) - f(:,MDIV-1)) / DM
    do s = 1, SDIV
      df_dm(s,2:MDIV-1) = (f(s,3:MDIV) - f(s,1:MDIV-2)) / (2.d0 * DM)
    end do
  end function deriv_m_vec

  function deriv_sm_vec(f) result(df_dsm)
    use para_mod, only : SDIV, MDIV, DS, DM
    real(8), dimension(SDIV,MDIV), intent(in) :: f
    real(8), dimension(SDIV,MDIV) :: df_dsm
    integer :: s, m

    ! Central difference for the interior
    df_dsm(2:SDIV-1, 2:MDIV-1) = (f(3:SDIV, 3:MDIV) - f(1:SDIV-2, 3:MDIV) &
            - f(3:SDIV, 1:MDIV-2) + f(1:SDIV-2, 1:MDIV-2)) / (4.d0 * DS * DM)

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
    real(8), intent(in)    :: r_e_old
    real(8), intent(out)   :: r_e_new, dif
    real(8), intent(out)   :: sphi_pole_h, gama_pole_h, rho_pole_h
    real(8), intent(out)   :: gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h
    real(8), intent(out)   :: sphi_center_h, gama_center_h, rho_center_h
    real(8) :: s_p, r_e_new_sq, grgr
    real(8), dimension(SDIV) :: gama_mu_1, gama_mu_0, rho_mu_1, rho_mu_0, ww_mu_0, sphi_mu_0, sphi_mu_1
    integer :: s

    rho_mu_0  = rho(:,1);    gama_mu_0 = gama(:,1);    sphi_mu_0 = sphi(:,1);    ww_mu_0 = ww(:,1)
    rho_mu_1  = rho(:,MDIV); gama_mu_1 = gama(:,MDIV); sphi_mu_1 = sphi(:,MDIV)
    
    s_p = r_ratio**(1.d0/dble(s_pwr)) / ( 1.d0 + r_ratio**(1.d0/dble(s_pwr)) ) 
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
      grgr = grgr + B_coup / 2.d0 * ( sphi_center_h**2 - sphi_pole_h**2 )
    end if

    r_e_new_sq = ( 2.d0 * ( h_center - enthalpy_min ) ) / grgr
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
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h
    real(8), intent(in) :: sphi_pole_h, sphi_equator_h, ww_equator_h
    real(8) :: metric_diff, term_in_Omega_h
    real(8), parameter :: TOLERANCE = 1.0d-3
    integer :: s, m
    if (abs(r_ratio - 1.0d0) < TOLERANCE) then
      Omega_c = 0.0d0; Omega_e = 0.0d0; Omg = 0.0d0
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
                  + B_coup / 2.d0 * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.d0 - exp( r_e_new**2 * metric_diff )
      if (term_in_Omega_h >= 0.d0) then
          Omega_c = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else
          write(*,"(10A15)") "gama_pole", "rho_pole", "gama_equator", "rho_equator", "sphi_pole", "sphi_equator"
          write(*,"(10es15.3)") gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, sphi_pole_h, sphi_equator_h
          stop "Omega can't be found; Line 99 of spin helper"
      end if
      Omg = Omega_c
      Omega_e = Omega_c
    end subroutine uniform_rotation
    subroutine const_j_rotation()
      real(8) :: guess, rsm, wwsm, mum, sgp
      real(8), parameter :: tolerance = 1.d-5
      metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                  + B_coup / 2.d0 * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.d0 - exp( r_e_new**2 * metric_diff )
      if (term_in_Omega_h >= 0.d0) then
        Omega_e = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else 
        stop "L172 in const_j_rotation"
      endif
      guess = Omega_e * 0.8d0
      call find_omege_e(guess, r_e_new, rho_equator_h, gama_equator_h, &
                       ww_equator_h, rho_pole_h, gama_pole_h, tolerance, Omega_e, diff_rotation_const_j)

      term_in_Omega_h = abs(Omega_e - ww_equator_h) * exp(-2.d0 * r_e_new**2 * rho_equator_h)
      Omega_c = Omega_e + term_in_Omega_h / (1.d0 - term_in_Omega_h * abs(Omega_e - ww_equator_h)) / A_diff**2
      
      Omg(1,:) = Omega_c
      Omg(1:2*SDIV/3,MDIV) = Omega_c
      do s = 2, 2*SDIV/3
        do m = 1, MDIV-1
          rsm = rho(s,m)
          wwsm = ww(s,m)
          mum = mu(m)
          sgp = s_gp(s)
          call zbrent_rot(Omg(s-1,m) * 8.d-1, r_e_new, rsm, wwsm, sgp, mum, 1.d-5, omg(s,m), rotation_law_const_j)
          F_j(s,m) = (omg(s,m) - wwsm) * sgp**2 * (1.d0 - mum**2) &
                / ((1.d0 - sgp)**2 * exp(2.d0 * r_e_new**2 * rsm) - (omg(s,m) - wwsm)**2 * sgp**2 * (1.d0 - mum**2))
        end do
      end do
    end subroutine const_j_rotation
    subroutine uryu_rotation()
      real(8) :: diff_Fmax, guess, Fa, rsm, wwsm, sgp, mum, omg_max_h
      real(8), dimension(SDIV) :: omg_mu_0
      integer :: imax
      real(8), parameter :: tolerance = 1.d-5
      metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                  + B_coup / 2.d0 * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.d0 - exp( r_e_new**2 * metric_diff )
      if (term_in_Omega_h >= 0.d0) then
        Omega_e = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else
        stop "L205 in uryu"
      endif
      Fmax_h = 2.d-2

      diff_Fmax = 1.d0
      Fmax_h    = Fmax_h / 2.d0
      do while( abs(diff_Fmax) > 1.d-7)
        guess = Omega_e
        call find_omege_e(guess, r_e_new,rho_equator_h,gama_equator_h,ww_equator_h, &
                        rho_pole_h,gama_pole_h, tolerance, Fa, diff_rotation_uryu) ! Fmax_h used here
        Omega_e = fa
        F_equator_h  = (Omega_e - ww_equator_h) / ( exp(2.d0*r_e_new**2*rho_equator_h) - (Omega_e-ww_equator_h)**2 )
        if ( F_equator_h < 0.d0 ) stop "negative F_equator_h; L120 in uryu"

        Omega_c = Omega_e / lambda2
        omg_mu_0(1) = Omega_c
        mum = 0.d0
        do s = 2, SDIV*2/3
            rsm = rho(s,1) ! hat
            wwsm= ww (s,1) ! hat
            sgp = s_gp(s)
            guess  = omg_mu_0(s-1)
            call zbrent_rot( guess, r_e_new, rsm, wwsm, sgp, mum, 1.d-5, omg_mu_0(s), rotation_law_uryu)
        enddo
        imax      = maxloc( omg_mu_0, 1 )
        omg_max_h = omg_mu_0(imax)
        diff_Fmax = ( lambda1 - omg_max_h / Omega_c )
        !write(*,*) "Shooting Fmax:", Fmax_h, diff_Fmax; stop 4
        Fmax_h    = Fmax_h - diff_Fmax * 1.d-2
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
          call zbrent_rot( guess, r_e_new, rsm, wwsm, sgp, mum, 1.d-5, omg(s,m), rotation_law_uryu)
          F_j(s,m) = (omg(s,m) - wwsm) * sgp**2 * (1.d0 - mum**2) &
                / ((1.d0 - sgp)**2 * exp(2.d0 * r_e_new**2 * rsm) - (omg(s,m) - wwsm)**2 * sgp**2 * (1.d0 - mum**2))
        enddo
      enddo
    end subroutine uryu_rotation
  end subroutine update_angular_velocity

  subroutine update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h)
    use rotation_law_mod, only: intF
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), dimension(SDIV,MDIV) :: sgp_term, sin_theta_2d, sgp_2d
    logical, dimension(SDIV,MDIV) :: valid

    sgp_term = spread(s_gp / (1.d0 - s_gp), 2, MDIV)
    sin_theta_2d = spread(sin_theta, 1, SDIV)
    velocity_sq = merge(0.d0, ((Omg - ww) * sgp_term * sin_theta_2d * exp(-rho * r_e_new**2))**2, r_ratio == 1.d0)

    where (velocity_sq > 1.d0) velocity_sq = 0.d0

    enthalpy = enthalpy_min + 0.5d0 * ( &
          r_e_new**2 * ( gama_pole_h + rho_pole_h - gama - rho + ( sphi**2 - sphi_pole_h**2 ) * B_coup / 2.d0 ) &
          - log( max(1.d-300, 1.d0-velocity_sq) )  )

    if ( trim(solver_type) == "const_j" ) then
      enthalpy = enthalpy + 0.5d0 * A_diff**2 * (Omg - Omega_c)**2
    elseif ( trim(solver_type) == "uryu" ) then
      enthalpy = enthalpy - intF( omg, F_j )
    endif

    sgp_2d = spread(s_gp, 2, MDIV)
    valid = (enthalpy > enthalpy_min) .and. (sgp_2d <= s_e)

    where (.not. valid)
      enthalpy = enthalpy_min
      pressure = 0.d0
      energy   = 0.d0
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
    real(8), intent(in)  :: r_e_new
    real(8), intent(in)  :: root_mphi_re
    real(8), intent(out) :: mr_cache(:), wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    real(8), intent(out) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(8), intent(out) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(8), intent(out) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(8), intent(out) :: e2alpha_r2_cache(:,:), Acoup_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m, n
    real(8) :: s1(SDIV), m1(MDIV)

    mr_cache(:)    = root_mphi_re * s_gp(:) / (1.d0 - s_gp(:))
    wfac_cache(:)  = 1.d0 / (1.d0 - s_gp(:))**2

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

    s1 = s_gp * (1.d0 - s_gp)
    m1 = 1.d0 - mu**2
    d2g_ss_cache = spread(s1, 2, MDIV) * deriv_s_vec(dg_s_cache) + spread(1.d0 - 2.d0 * s_gp, 2, MDIV) * dg_s_cache
    d2g_mm_cache = spread(m1, 1, SDIV) * deriv_m_vec(dg_m_cache) - 2.d0 * spread(mu, 1, SDIV) * dg_m_cache
    e_gsm_cache      = exp(0.5d0 * gama)
    e_rsm_cache      = exp(-rho)
    e2alpha_r2_cache = exp(2.d0 * alpha) * r_e_new**2
    Acoup_cache      = exp(-sphi**2 * B_coup / 4.d0)
    Acoup4_cache     = Acoup_cache**4
  end subroutine precompute_derivatives_and_bessels

  subroutine build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, &
                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                e2alpha_r2_cache, Acoup4_cache)
    real(8), intent(in)  :: r_e_new
    real(8), intent(out) :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(8), intent(in)  :: dr_s_cache(:,:), dr_m_cache(:,:), dg_s_cache(:,:), dg_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(8), intent(in)  :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(8), intent(in)  :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(8), intent(in)  :: e2alpha_r2_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m
    real(8), dimension(SDIV,MDIV) :: esm, psm, Vphi, scal_p
    real(8), dimension(SDIV,MDIV) :: s1, s2, m1, mum

    s1 = spread(s_gp * (1.d0 - s_gp), 2, MDIV)
    s2 = spread((s_gp / (1.d0 - s_gp))**2, 2, MDIV)
    mum = spread(mu, 1, SDIV)
    m1 = 1.d0 - mum**2

    esm  = energy   * Acoup4_cache
    psm  = pressure * Acoup4_cache
    Vphi = sphi**2 * mphi_r * 0.5d0 * e2alpha_r2_cache
    scal_p = sphi * e_gsm_cache

    S_metric_rho = e_gsm_cache * ( &
        16.d0 * pi * e2alpha_r2_cache / 2.d0 * (esm + psm) * s2 * (1.d0+velocity_sq) / (1.d0-velocity_sq) &
      + s2 * m1 * e_rsm_cache**2 * ( (s1*dww_s_cache)**2 + m1*dww_m_cache**2 ) &
      + s1 * dg_s_cache - mum * dg_m_cache &
      + rho / 2.d0 * ( ( 16.d0 * pi * e2alpha_r2_cache * psm - 4.d0 * Vphi ) * s2 &
        - s1 * dg_s_cache * ( s1 / 2.d0 * dg_s_cache + 1.d0 ) &
        - dg_m_cache * ( m1 / 2.d0 * dg_m_cache - mum ) ) )

    S_metric_gama = e_gsm_cache * ( ( 16.d0 * pi * e2alpha_r2_cache * psm - 4.d0 * Vphi ) * s2 &
      + gama / 2.d0 * ( ( 16.d0 * pi * e2alpha_r2_cache * psm - 4.d0 * Vphi ) * s2 &
        - ( s1 * dg_s_cache )**2 / 2.d0 - m1 * dg_m_cache**2 / 2.d0 ) )

    S_metric_omega = e_gsm_cache * e_rsm_cache * ( &
      - 16.d0 * pi * e2alpha_r2_cache * ( Omg - ww ) * ( esm + psm ) * s2 / (1.d0-velocity_sq) &
      + ww * ( -0.5d0 * 16.d0 * pi * e2alpha_r2_cache * s2 * &
        ( ( ( 1.d0 + velocity_sq ) * esm + 2.d0 * velocity_sq*psm ) / ( 1.d0 - velocity_sq ) ) &
        - s1 * ( 2.d0 * dr_s_cache + 0.5d0 * dg_s_cache ) &
        + mum * ( 2.d0 * dr_m_cache + 0.5d0 * dg_m_cache) &
        + 0.25d0 * s1**2 * ( 4.d0 * dr_s_cache**2 - dg_s_cache**2) &
        + 0.25d0 * m1 * ( 4.d0 * dr_m_cache**2 - dg_m_cache**2 ) &
        - m1 * e_rsm_cache**2 * ( spread(s_gp,2,MDIV)**4 * dww_s_cache**2 + s2 * m1 * dww_m_cache**2) &
        - 2.d0 * Vphi * s2 ) )

    if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.d10 ) then
      S_metric_sphi = - r_e_new**2 * s2 * scal_p * mphi_r &
        + s2 * scal_p * e2alpha_r2_cache * ( -2.d0 * pi * B_coup * ( esm - 3.d0 * psm ) + mphi_r ) &
        + scal_p * ( s1 * spread(1.d0-s_gp,2,MDIV) * dg_s_cache       &
        + s1 * s1 * ( d2g_ss_cache * 0.5d0 + dg_s_cache**2 * 0.25d0 ) &
        + m1      * ( d2g_mm_cache * 0.5d0 + dg_m_cache**2 * 0.25d0 ) &
        - mum * dg_m_cache )
      !> v1
      !S_metric_sphi = - r_e_new**2 * s2 * scal_p * mphi_r &
      !    + r_e_new**2 * s2 * scal_p * exp(2.d0*alpha) * ( -2.d0 * pi * B_coup * (esm-3.d0*psm) + mphi_r) &
      !    + scal_p * ( - 2.d0 * s1 * spread(s_gp,2,MDIV) * dg_s_cache &
      !    + s1**2 * ( d2g_ss_cache + dg_s_cache**2 / 4.d0 ) + s1 * dg_s_cache &
      !    + m1 * d2g_mm_cache / 2.d0 + dg_m_cache**2 * r_e_new**2 * s2 / 4.d0 &
      !    - 2.d0 * mum * dg_m_cache )
    else
      S_metric_sphi = -s1**2 * dg_s_cache * ds_s_cache - m1 * dg_m_cache * ds_m_cache &
        + sphi * ( -2.d0 * pi * B_coup * (esm - 3.d0*psm) ) * s2 * e2alpha_r2_cache
    endif
  end subroutine build_source_terms

  subroutine angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                 D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    real(8), intent(in)  :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(8), intent(out) :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    integer :: n, k, ifail
    real(8) :: sum_rho, sum_sphi, sum_gama, sum_omega, er2
    real(8), dimension(MDIV) :: Int_m

    n = 0
    do k = 1, SDIV
      Int_m(:) = P_2n(:,n+1) * S_metric_rho(k,:)
      call d01gaf(mu, Int_m, MDIV, sum_rho, er2, ifail)

      Int_m(:) = P_2n(:,n+1) * S_metric_sphi(k,:)
      call d01gaf(mu, Int_m, MDIV, sum_sphi, er2, ifail)

      D1_metric_rho  (n+1,k) = sum_rho
      D1_metric_sphi (n+1,k) = sum_sphi
      D1_metric_gama (n+1,k) = 0.d0
      D1_metric_omega(n+1,k) = 0.d0
    enddo

    do n = 1, LMAX
      do k = 1, SDIV
        Int_m(:) = P_2n(:,n+1) * S_metric_rho(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_rho, er2, ifail)

        Int_m(:) = P_2n(:,n+1) * S_metric_sphi(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_sphi, er2, ifail)

        Int_m(:) = sin_2n_1_theta(:,n) * S_metric_gama(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_gama, er2, ifail)

        Int_m(:) = sin_theta(:) * P1_2n_1(:,n+1) * S_metric_omega(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_omega, er2, ifail)

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
    real(8), intent(in)  :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    real(8), intent(out) :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
    real(8), intent(in)  :: root_mphi_re
    real(8), intent(in)  :: wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    integer :: s, n, ifail
    real(8) :: sum_val, er2
    real(8), dimension(SDIV) :: Int_s

    n = 0
    do s = 1, SDIV
      Int_s = f_rho(s,1,:) * D1_metric_rho(n+1,:)
      call d01gaf(s_gp, Int_s, SDIV, D2_metric_rho(s,n+1), er2, ifail)
    end do
    D2_metric_gama (:,n+1) = 0.d0
    D2_metric_omega(:,n+1) = 0.d0

    if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.d10 ) then
      do n = 0, LMAX
        do s = 1, SDIV
          Int_s(1:s-1 ) = besseli_cache(n+1,1:s-1) * besselk_cache(n+1,s)  * wfac_cache(1:s-1)  * D1_metric_sphi(n+1,1:s-1 )
          Int_s(s:SDIV) = besseli_cache(n+1,s) * besselk_cache(n+1,s:SDIV) * wfac_cache(s:SDIV) * D1_metric_sphi(n+1,s:SDIV)
          Int_s = Int_s * root_mphi_re
          call d01gaf(s_gp, Int_s, SDIV, D2_metric_sphi(s,n+1), er2, ifail)
        end do
      end do
    else
      do n = 0, LMAX
        do s = 1, SDIV
          Int_s = f_rho(s,n+1,:) * D1_metric_sphi(n+1,:)
          call d01gaf(s_gp, Int_s, SDIV, D2_metric_sphi(s,n+1), er2, ifail)
        end do
      end do
    end if

    do n = 1, LMAX
      do s = 1, SDIV
        ! rho
        Int_s = f_rho(s,n+1,:) * D1_metric_rho(n+1,:)
        call d01gaf(s_gp, Int_s, SDIV, D2_metric_rho(s,n+1), er2, ifail)

        ! gama
        Int_s = f_gama(s,n+1,:) * D1_metric_gama(n+1,:)
        call d01gaf(s_gp, Int_s, SDIV, D2_metric_gama(s,n+1), er2, ifail)

        ! omega
        Int_s(1:s-1) = f_rho(s,n+1,1:s-1) * D1_metric_omega(n+1,1:s-1)
        Int_s(s:SDIV) = f_gama(s,n+1,s:SDIV) * D1_metric_omega(n+1,s:SDIV)
        call d01gaf(s_gp, Int_s, SDIV, D2_metric_omega(s,n+1), er2, ifail)
      end do
    end do

  end subroutine radial_integration

  subroutine sum_coefficients_and_get_targets(target_rho, target_gama, target_ww, target_sphi, &
                                              D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    real(8), intent(out) :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
    real(8), intent(in)  :: D2_metric_rho(SDIV,LMAX+1), D2_metric_gama(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1), D2_metric_sphi(SDIV,LMAX+1)
    integer :: n
    real(8), dimension(SDIV,MDIV) :: exp_mhalf_gsm, exp_rsm_mhalf_gsm
    real(8), dimension(SDIV,MDIV) :: sum_rho, sum_sphi, sum_gama, sum_omega
    real(8), dimension(MDIV) :: sin_theta_inv

    exp_mhalf_gsm = exp(-0.5d0 * gama)
    exp_rsm_mhalf_gsm = exp(rho - 0.5d0 * gama)

    ! Monopole term (n=0)
    sum_rho  = -exp_mhalf_gsm * matmul(reshape(D2_metric_rho(:,1), [SDIV, 1]), reshape(P_2n(:,1), [1, MDIV]))
    sum_gama = 0.d0
    sum_omega = 0.d0

    ! Precompute 1/sin(theta) to avoid division in the loop
    sin_theta_inv = 0.d0; where (sin_theta > 1.d-12) sin_theta_inv = 1.d0 / sin_theta

    ! Higher multipoles (n=1 to LMAX)
    do n = 1, LMAX
      sum_rho  = sum_rho  - exp_mhalf_gsm * matmul(reshape(D2_metric_rho(:,n+1), [SDIV, 1]), reshape(P_2n(:,n+1), [1, MDIV]))
      
      ! Vectorized update for sum_gama and sum_omega
      sum_gama(:,1:MDIV-1) = sum_gama(:,1:MDIV-1) - (2.d0/pi) * exp_mhalf_gsm(:,1:MDIV-1) * &
            (matmul(reshape(D2_metric_gama(:,n+1), [SDIV, 1]), reshape(sin_2n_1_theta(1:MDIV-1,n) * sin_theta_inv(1:MDIV-1), [1, MDIV-1])) / (2.d0*n-1.d0))
      sum_gama(:,MDIV) = sum_gama(:,MDIV) - (2.d0/pi) * exp_mhalf_gsm(:,MDIV) * D2_metric_gama(:,n+1)

      sum_omega(:,1:MDIV-1) = sum_omega(:,1:MDIV-1) - exp_rsm_mhalf_gsm(:,1:MDIV-1) * &
            (matmul(reshape(D2_metric_omega(:,n+1), [SDIV, 1]), reshape(P1_2n_1(1:MDIV-1,n+1) * sin_theta_inv(1:MDIV-1), [1, MDIV-1])) / (2.d0*n*(2.d0*n-1.d0)))
      sum_omega(:,MDIV) = sum_omega(:,MDIV) + exp_rsm_mhalf_gsm(:,MDIV) * D2_metric_omega(:,n+1) / 2.d0
    end do

    if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.d10 ) then
      sum_sphi = -exp_mhalf_gsm * matmul(reshape(D2_metric_sphi(:,1), [SDIV, 1]), reshape(P_2n(:,1), [1, MDIV]))
      do n = 1, LMAX
        sum_sphi = sum_sphi - exp_mhalf_gsm * (2.d0*dble(n)+1.d0) &
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
    real(8), intent(in) :: r_e_new
    real(8), intent(in) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(8), intent(in) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(8), intent(in) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_rsm_cache(:,:)
    integer :: s, m
    real(8) :: sgp, s1, sgp_ratio
    real(8), dimension(SDIV,MDIV) :: da_dm, d_gama_sm_all
    real(8), dimension(MDIV) :: m1_vec, d_gama_s_row, d_gama_m_row, d_rho_s_row, d_rho_m_row
    real(8), dimension(MDIV) :: d_sphi_s_row, d_sphi_m_row, d_gama_sm_row, d_ww_s_row, d_ww_m_row, d_gama_ss_row, d_gama_mm_row
    real(8), dimension(MDIV) :: temp1_row, temp2_row, temp3_row, temp4_row, temp5_row, temp6_row, temp7_row, temp8_row, temp9_row
    real(8), dimension(SDIV) :: s1_cache, sgp_ratio_cache
    real(8) :: adj_const(SDIV)

    alpha(:,:) = 0.d0
    if (r_ratio == 1.d0) then
      da_dm(:,:) = 0.0d0
    else
      s1_cache = s_gp * (1.d0-s_gp)
      sgp_ratio_cache = s_gp / (1.d0-s_gp)

      da_dm(1,:) = 0.0d0
      d_gama_sm_all = deriv_sm_vec(gama)
      m1_vec = 1.d0 - mu**2
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

        temp1_row = 2.d0 * sgp**2 * sgp_ratio * m1_vec * d_ww_s_row * d_ww_m_row * (1.d0 + s1 * d_gama_s_row) &
          - ( (sgp**2 * d_ww_s_row)**2 - (sgp * d_ww_m_row * sgp_ratio)**2 * m1_vec ) * (-mu + m1_vec * d_gama_m_row)
        temp2_row = 1.d0 / ( m1_vec * (1.d0 + s1 * d_gama_s_row)**2 + (-mu + m1_vec * d_gama_m_row)**2 )
        temp3_row = s1 * d_gama_ss_row + (s1 * d_gama_s_row)**2
        temp4_row = d_gama_m_row * (-mu + m1_vec * d_gama_m_row)
        temp5_row = ( (s1 * (d_rho_s_row + d_gama_s_row))**2 - m1_vec * (d_rho_m_row + d_gama_m_row)**2 ) * (-mu + m1_vec * d_gama_m_row)
        temp6_row = s1 * m1_vec * (  (d_rho_s_row + d_gama_s_row) * (d_rho_m_row + d_gama_m_row) / 2.d0 + d_gama_sm_row + d_gama_s_row * d_gama_m_row  ) * (1.d0 + s1 * d_gama_s_row)
        temp7_row = s1 * mu * d_gama_s_row * ( 1.d0 + s1 * d_gama_s_row )
        temp8_row = m1_vec * (e_rsm_cache(s,:)**2)
        temp9_row = -temp2_row * (-mu + m1_vec * d_gama_m_row) * ( (s1 * d_sphi_s_row)**2 - m1_vec * d_sphi_m_row**2 ) &
              - m1_vec * s1 * ( 1.d0 + s1 * d_gama_s_row ) * 2.d0 * d_sphi_m_row * d_sphi_s_row

        da_dm(s,:) = - (d_rho_m_row + d_gama_m_row) / 2.d0 &
          - temp2_row * ( (temp3_row - d_gama_mm_row - temp4_row) * (-mu + m1_vec * d_gama_m_row) / 2.d0 &
          + temp5_row / 4.d0 - temp6_row  + temp7_row + temp8_row * temp1_row / 4.d0 ) + temp9_row
      end do
    end if

    alpha(:, 1) = 0.0d0

    do m = 1, MDIV-1
        alpha(:,m+1) = alpha(:,m) + dm * ( da_dm(:,m+1) + da_dm(:,m) ) * 0.5d0
    enddo
    
    alpha(SDIV,:) = 0.d0
    adj_const = alpha(:,MDIV) - ( gama(:,MDIV) - rho(:,MDIV) )/2.d0
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
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re
    real(8), intent(out) :: out_target_rho(SDIV,MDIV), out_target_gama(SDIV,MDIV), out_target_ww(SDIV,MDIV), out_target_sphi(SDIV,MDIV)
   
    call precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, wfac_cache, besseli_cache, besselk_cache, &
         dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
         d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup_cache, Acoup4_cache)

    call build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
         dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
         d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup4_cache)

    call angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
         D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)

    call radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
         D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, root_mphi_re, wfac_cache, &
         besseli_cache, besselk_cache)

    call sum_coefficients_and_get_targets(out_target_rho, out_target_gama, out_target_ww, out_target_sphi, &
         D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    
  end subroutine get_all_targets

  subroutine relaxation(r_e_new, target_rho, target_gama, target_ww, target_sphi, root_mphi_re, n_of_it)
    use hybrid_relaxation, only: hybrid_update
    use anderson_optimized, only: anderson_accel_optimized
    real(8), intent(in) :: r_e_new
    real(8), intent(in) :: target_rho(SDIV,MDIV), target_gama(SDIV,MDIV)
    real(8), intent(in) :: target_ww(SDIV,MDIV),  target_sphi(SDIV,MDIV)
    real(8), intent(in) :: root_mphi_re
    integer, intent(in) :: n_of_it
    integer :: s, m
    real(8) :: s_p, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8) :: rho_mu_1(SDIV), gama_mu_1(SDIV), sphi_mu_1(SDIV)
    integer, parameter :: m_hist = 3
    real(8), allocatable, save :: hist_f_rho(:,:,:)
    real(8), allocatable, save :: hist_f_gama(:,:,:)
    real(8), allocatable, save :: hist_f_ww(:,:,:)
    real(8), allocatable, save :: hist_f_sphi(:,:,:)
    real(8), allocatable, save :: hist_x_rho(:,:,:)
    real(8), allocatable, save :: hist_x_gama(:,:,:)
    real(8), allocatable, save :: hist_x_ww(:,:,:)
    real(8), allocatable, save :: hist_x_sphi(:,:,:)

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
    case ('newton')
        ! Compute pole values by interpolation along the outer boundary (mu = 1)
        rho_mu_1  = rho (:,MDIV)
        gama_mu_1 = gama(:,MDIV)
        sphi_mu_1 = sphi(:,MDIV)
        s_p = r_ratio**(1.d0/dble(s_pwr)) / (1.d0 + r_ratio**(1.d0/dble(s_pwr)))
        call interp(s_gp, rho_mu_1,  SDIV, s_p, rho_pole_h)
        call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole_h)
        call interp(s_gp, sphi_mu_1, SDIV, s_p, sphi_pole_h)

        call newton_krylov('rho',  rho,  n_of_it, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
        call newton_krylov('gama', gama, n_of_it, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
        call newton_krylov('ww',   ww,   n_of_it, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
        call newton_krylov('sphi', sphi, n_of_it, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    case default
        stop "Unknown relaxation scheme specified"
    end select

    where(sphi .ne. sphi) sphi = 0.d0

    if ( mphi_r == 0.d0 ) return
    
    do m=1,MDIV
      do s=2,SDIV
      if (sphi(s,m) < 0.d0) then
          sphi(s,m) = sphi(s-1,m) * exp(root_mphi_re * (s_gp(s-1)/(1.d0-s_gp(s-1)) - s_gp(s)/(1.d0-s_gp(s))))
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
    allocate(target_field_perturbed(SDIV,MDIV))
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
    deallocate(target_field_perturbed)
  end subroutine deallocate_workspace

  subroutine output_helper(D2_metric_rho, D2_metric_omega)
    real(8), intent(in) :: D2_metric_rho(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1)
    real(8) :: r_inf, M2, S3, M4, rho_0
    character(32) :: fil1, fil2, fil3, fil4, fil5, fil6
    integer :: s, m
    external :: mass_radius
    real(8), external :: n0_at_e

    r_inf = r_e * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1.d0 - s_gp(SDIV - 1) ))**s_pwr
    M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2.d0 * r_inf**3 * ( C**2 / G / Mass )**3
    S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2.d0 * r_inf**5 * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2.d0 * r_inf**5 * ( C**2 / G / Mass )**5
    
    write(*,"(A)") " ", " Off-loading data ..."
    
    write(fil1,"(f6.2)") ang_mom
    write(fil2,"(f16.5)") mass_0/MSUN
    write(fil3,"(es15.2)") B_coup
    write(fil4,"(es15.2)") sqrt(mphi_r*1.d10/KAPPA)*l_uni
    rho_0 = n0_at_e( energy(1,1) ) * MB
    write(fil5,"(es15.3)") rho_0
    write(fil6,"(f15.3)") sphi_m

    call write_output_file("./Cont/"//trim(adjustl(eos_file))//"_J_"//trim(adjustl(fil1))//&
      "_Mb"//trim(adjustl(fil2))//"_B"//trim(adjustl(fil3))//&
      "_mphi"//trim(adjustl(fil4))//"_rhoc"//trim(adjustl(fil5))//"_sphim"//trim(adjustl(fil6))//".dat")

    call write_output_file("./Res/res.dat")
  end subroutine output_helper

  subroutine write_output_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios, s, m
    real(8) :: rho_0_val
    character(len=256) :: header_fmt, data_fmt
    real(8), external :: n0_at_e

    ! Define formats for clarity and easy modification
    header_fmt = "(3i5, 5es27.17e3)"
    data_fmt = "(13es27.17e3)"

    ! Use newunit to get a free file unit, preventing conflicts
    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim(filename)
      return
    end if

    ! Write header
    write(unit, fmt=header_fmt) SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.d5, &
            energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)), Omega_c* (C/sqrt(kappa))

    ! Write data grid
    do s = 1, SDIV
      do m = 1, MDIV
        if (enthalpy(s,m) > enthalpy_min) then
          rho_0_val = n0_at_e(energy(s,m)) * MB
        else
          rho_0_val = 0.d0
        end if

        write(unit, fmt=data_fmt) s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), &
          ww(s,m) * (C/sqrt(kappa)), pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), &
          enthalpy(s,m), rho_0_val, velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), &
          sphi(s,m) * sqrt(B_coup)
      end do
    end do

    close(unit)
  end subroutine write_output_file

  ! ---------------------------

  subroutine residual(field_name, x, F_x, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(in) :: x
    real(8), dimension(SDIV,MDIV), intent(out) :: F_x
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re
    real(8), dimension(SDIV,MDIV) :: G_x

    call evaluate_field_target(field_name, x, G_x, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)

    F_x = G_x - x
  end subroutine residual

  subroutine jacobian_vector_product(field_name, x, F_x, v, Jv, eps, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(in) :: x
    real(8), dimension(SDIV,MDIV), intent(in) :: F_x
    real(8), dimension(SDIV*MDIV), intent(in) :: v
    real(8), dimension(SDIV*MDIV), intent(out) :: Jv
    real(8), intent(in) :: eps
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re

    real(8), dimension(SDIV,MDIV) :: F_x_plus_eps_v
    real(8), dimension(SDIV,MDIV) :: x_plus_eps_v
    real(8) :: norm_v, h

    norm_v = sqrt(sum(v**2))
    if (norm_v < 1.d-12) then
        Jv = 0.d0
        return
    end if

    h = eps / norm_v

    ! F(x + h*v)
    x_plus_eps_v = x + reshape(v, [SDIV, MDIV]) * h
    call residual(field_name, x_plus_eps_v, F_x_plus_eps_v, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)

    ! J*v
    Jv = reshape(F_x_plus_eps_v - F_x, [SDIV*MDIV]) / h
  end subroutine jacobian_vector_product

  subroutine newton_krylov(field_name, current_field, iter, r_e_new, root_mphi_re, &
                           gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(inout) :: current_field
    integer, intent(in) :: iter
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re

    integer, parameter :: gmres_m = 3
    integer, parameter :: max_iter = 50
    real(8), parameter :: epsilon = 1.d-6
    real(8), parameter :: step_clip = 0.5d0
    real(8), parameter :: armijo_c1 = 1.d-4
    real(8), parameter :: gmres_tol_min = 1.d-6
    real(8), parameter :: gmres_tol_scale = 1.d-2
    integer, parameter :: conservative_steps = 3

    integer :: N, restart_iter, i, j, k
    real(8), allocatable, save :: residual_vec(:), v(:), w(:), delta_x(:), s(:)
    real(8), allocatable, save :: V_basis(:,:), H(:,:), g(:), cs(:), sn(:)
    real(8), dimension(SDIV,MDIV) :: F_x
    real(8) :: norm_residual, beta, r_norm, step_scale
    real(8) :: phi_old, phi_new, alpha
    real(8) :: local_gmres_tol, gmres_target
    real(8), dimension(SDIV,MDIV) :: target_current
    real(8), dimension(SDIV,MDIV) :: trial_field, trial_residual_reshaped
    real(8) :: best_r_norm

    N = SDIV * MDIV
    if (.not. allocated(residual_vec) .or. size(residual_vec) /= N) then
      if (allocated(residual_vec)) then
        deallocate(residual_vec, v, w, delta_x, s, V_basis, H, g, cs, sn)
      end if
      allocate(residual_vec(N), v(N), w(N), delta_x(N), s(gmres_m))
      allocate(V_basis(N, gmres_m+1))
      allocate(H(gmres_m+1, gmres_m))
      allocate(g(gmres_m+1), cs(gmres_m), sn(gmres_m))
    end if

    if (iter < conservative_steps) then
      call evaluate_field_target(field_name, current_field, target_current, r_e_new, root_mphi_re, &
           gama_pole_h, rho_pole_h, sphi_pole_h)
      current_field = 0.9d0 * current_field + 0.1d0 * target_current
      return
    end if

    ! 1. Compute residual F(x) = G(x) - x
    call residual(field_name, current_field, F_x, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    residual_vec = reshape(F_x, [N])
    norm_residual = sqrt(sum(residual_vec**2))
    local_gmres_tol = max(gmres_tol_min, gmres_tol_scale * norm_residual)
    gmres_target = local_gmres_tol * norm_residual
    phi_old = 0.5d0 * norm_residual**2
    best_r_norm = norm_residual

    if (norm_residual < 1.d-12) return

    ! 2. Solve J * delta_x = -F_x using GMRES
    delta_x = 0.d0
    r_norm = norm_residual
    
    outer_gmres: do restart_iter = 1, max_iter
      k = gmres_m
      V_basis(:,1) = -residual_vec / r_norm
      g = 0.d0
      g(1) = r_norm
      H = 0.d0

      inner_gmres: do j = 1, gmres_m
        v = V_basis(:,j)
        call jacobian_vector_product(field_name, current_field, F_x, v, w, epsilon, &
             r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)

        do i = 1, j
          H(i,j) = dot_product(w, V_basis(:,i))
          w = w - H(i,j) * V_basis(:,i)
        end do

        beta = sqrt(sum(w**2))
        H(j+1,j) = beta
        if (beta < 1.d-14) then
          k = j
          exit inner_gmres
        end if
        V_basis(:,j+1) = w / H(j+1,j)

        do i = 1, j-1
          call apply_givens(cs(i), sn(i), H(i,j), H(i+1,j))
        end do

        call generate_givens(H(j,j), beta, cs(j), sn(j))
        H(j,j) = cs(j)*H(j,j) + sn(j)*beta
        g(j+1) = -sn(j)*g(j)
        g(j)   =  cs(j)*g(j)

        r_norm = abs(g(j+1))
        if (r_norm < gmres_target) then
          k = j
          exit inner_gmres
        end if
      end do inner_gmres

      s(1:k) = g(1:k)
      call DTRSV('U', 'N', 'N', k, H, gmres_m+1, s, 1)
      delta_x = delta_x + matmul(V_basis(:,1:k), s(1:k))

      best_r_norm = min(best_r_norm, r_norm)
      if (r_norm < gmres_target) exit outer_gmres
      ! Early exit if GMRES stalls and offers little improvement
      if (restart_iter > 1 .and. r_norm > 0.95d0 * norm_residual .and. r_norm > 0.99d0 * best_r_norm) then
        exit outer_gmres
      end if
    end do outer_gmres

    ! 3. Line search and update (Armijo backtracking on 0.5*||F||^2)
    step_scale = min(1.d0, step_clip / max(1.d-12, sqrt(sum(delta_x**2)/dble(N))))
    delta_x = step_scale * delta_x
    alpha = 1.d0

    do i = 1, 3
      trial_field = current_field + reshape(delta_x, [SDIV,MDIV])
      call residual(field_name, trial_field, trial_residual_reshaped, r_e_new, &
           root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
     residual_vec = reshape(trial_residual_reshaped, [N])
      phi_new = 0.5d0 * sum(residual_vec**2)
      if (phi_new <= phi_old - armijo_c1 * alpha * norm_residual**2) then
        current_field = trial_field
        return
      else
        delta_x = 0.5d0 * delta_x
        alpha = 0.5d0 * alpha
      end if
    end do

    ! Fallback: simple mixing
    call evaluate_field_target(field_name, current_field, target_current, r_e_new, root_mphi_re, &
         gama_pole_h, rho_pole_h, sphi_pole_h)
    current_field = 0.75d0 * current_field + 0.25d0 * target_current
  end subroutine newton_krylov

  subroutine apply_givens(c, s, v1, v2)
    real(8), intent(in) :: c, s
    real(8), intent(inout) :: v1, v2
    real(8) :: temp
    temp =  c*v1 + s*v2
    v2   = -s*v1 + c*v2
    v1   = temp
  end subroutine apply_givens

  subroutine generate_givens(a, b, c, s)
    real(8), intent(in) :: a, b
    real(8), intent(out) :: c, s
    real(8) :: r
    if (b == 0.d0) then
      c = 1.d0
      s = 0.d0
    else if (abs(b) > abs(a)) then
      r = a / b
      s = 1.d0 / sqrt(1.d0 + r*r)
      c = r * s
    else
      r = b / a
      c = 1.d0 / sqrt(1.d0 + r*r)
      s = r * c
    end if
  end subroutine generate_givens

  subroutine evaluate_field_target(field_name, trial_field, target_out, r_e_new, root_mphi_re, &
                                   gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(in) :: trial_field
    real(8), dimension(SDIV,MDIV), intent(out) :: target_out
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h 
    real(8), intent(in) :: root_mphi_re
    real(8), allocatable, save :: backup(:,:), temp_target_rho(:,:), temp_target_gama(:,:), temp_target_ww(:,:), temp_target_sphi(:,:)

    if (.not. allocated(backup)) then
      allocate(backup(SDIV,MDIV))
      allocate(temp_target_rho(SDIV,MDIV), temp_target_gama(SDIV,MDIV))
      allocate(temp_target_ww(SDIV,MDIV), temp_target_sphi(SDIV,MDIV))
    end if

    select case (trim(field_name))
    case ('rho')
      backup = rho
      rho = trial_field
    case ('gama')
      backup = gama
      gama = trial_field
    case ('ww')
      backup = ww
      ww = trial_field
    case ('sphi')
      backup = sphi
      sphi = trial_field
    case default
      stop 'evaluate_field_target: unknown field name'
    end select

    call get_all_targets(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, root_mphi_re, &
         temp_target_rho, temp_target_gama, temp_target_ww, temp_target_sphi)

    select case (trim(field_name))
    case ('rho')
      target_out = temp_target_rho
      rho = backup
    case ('gama')
      target_out = temp_target_gama
      gama = backup
    case ('ww')
      target_out = temp_target_ww
      ww = backup
    case ('sphi')
      target_out = temp_target_sphi
      sphi = backup
    end select
  end subroutine evaluate_field_target
end module spin_helper
