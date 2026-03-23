module spin_updates
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  private
  public :: update_equatorial_radius, update_angular_velocity, update_eos_and_velocity

contains

  subroutine update_equatorial_radius(r_e_old, r_e_new, dif, sphi_pole_h, gama_pole_h, rho_pole_h, &
                                      gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h, &
                                      sphi_center_h, gama_center_h, rho_center_h)
    use para_mod, only : wp, SDIV, MDIV, s_gp, s_e, s_pwr, r_ratio, &
                         sphi, gama, rho, ww, has_scalar, B_coup, h_center, enthalpy_min
    use toolkit_mod, only : interp
    real(wp), intent(in)    :: r_e_old
    real(wp), intent(out)   :: r_e_new, dif
    real(wp), intent(out)   :: sphi_pole_h, gama_pole_h, rho_pole_h
    real(wp), intent(out)   :: gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h
    real(wp), intent(out)   :: sphi_center_h, gama_center_h, rho_center_h
    real(wp) :: r_e_new_sq, grgr
    real(wp), save :: s_p_cached = -1.0_wp
    real(wp), save :: r_ratio_prev = -1.0_wp

    if (r_ratio /= r_ratio_prev) then
      r_ratio_prev = r_ratio
      s_p_cached = r_ratio**(1.0_wp/dble(s_pwr)) / (1.0_wp + r_ratio**(1.0_wp/dble(s_pwr)))
    end if

    call interp(s_gp, sphi(:,MDIV), SDIV, s_p_cached, sphi_pole_h   )
    call interp(s_gp, gama(:,MDIV), SDIV, s_p_cached, gama_pole_h   )
    call interp(s_gp, rho(:,MDIV),  SDIV, s_p_cached, rho_pole_h    )
    call interp(s_gp, gama(:,1),    SDIV, s_e,        gama_equator_h)
    call interp(s_gp, rho(:,1),     SDIV, s_e,        rho_equator_h )
    call interp(s_gp, ww(:,1),      SDIV, s_e,        ww_equator_h  )
    call interp(s_gp, sphi(:,1),    SDIV, s_e,        sphi_equator_h)
    sphi_center_h = sphi(1,1)
    gama_center_h = gama(1,1)
    rho_center_h  = rho(1,1)

    grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h
    if (has_scalar) then
      grgr = grgr + B_coup / 2.0_wp * ( sphi_center_h**2 - sphi_pole_h**2 )
    end if

    r_e_new_sq = ( 2.0_wp * ( h_center - enthalpy_min ) ) / grgr
    r_e_new = sqrt( r_e_new_sq )
    dif = abs(r_e_old - r_e_new) / r_e_new

    if (r_e_new / r_e_old > 2 .or. isnan(r_e_new) ) then
      write(*,*) "r_e_old :", r_e_old
      write(*,*) "r_e_new :", r_e_new
      write(*,*) "grgr    :", grgr
      write(*,*) "r_e_new_sq :", r_e_new_sq
      stop 'r_e cannot be found.'
    endif
  end subroutine update_equatorial_radius

  subroutine update_angular_velocity(r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, &
                                    sphi_pole_h, sphi_equator_h, ww_equator_h)
    use para_mod, only : wp, SDIV, MDIV, s_gp, mu, rho, ww, r_ratio, &
                         Omega_c, Omega_e, Omg, F_j, has_scalar, B_coup, A_diff, solver_type, &
                         lambda1, lambda2, Fmax_h, F_equator_h
    use rotation_law_mod, only: diff_rotation_const_j, rotation_law_const_j, &
                                diff_rotation_uryu, rotation_law_uryu
    use brent_mod, only : find_omege_e, zbrent_rot
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h
    real(wp), intent(in) :: sphi_pole_h, sphi_equator_h, ww_equator_h
    real(wp) :: metric_diff, term_in_Omega_h, re2
    real(wp), parameter :: TOLERANCE = 1.0e-3_wp
    integer :: s, m

    if (abs(r_ratio - 1.0e0_wp) < TOLERANCE) then
      Omega_c = 0.0e0_wp; Omega_e = 0.0e0_wp; Omg = 0.0e0_wp
      return
    end if

    re2 = r_e_new**2
    metric_diff = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h &
                + B_coup / 2.0_wp * ( sphi_equator_h**2 - sphi_pole_h**2 )
    term_in_Omega_h = 1.0_wp - exp( re2 * metric_diff )
    if (term_in_Omega_h >= 0.0_wp) then
      Omega_e = ww_equator_h + exp(re2 * rho_equator_h) * sqrt(term_in_Omega_h)
    else
      write(*,"('Solving for axis ratio: ', f12.5)") r_ratio
      write(*,"(10A15)") "gama_pole", "rho_pole", "gama_equator", "rho_equator", "sphi_pole", "sphi_equator"
      write(*,"(10es15.3)") gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, sphi_pole_h, sphi_equator_h
      stop "Omega can't be found; Line 99 of spin helper"
    endif

    select case(trim(solver_type))
    case("uniform")
      Omega_c = Omega_e; Omg = Omega_e
    case("const_j")
      call const_j_rotation()
    case("uryu")
      call uryu_rotation()
    case default
      stop "Unknown solver type"
    end select
  contains
    subroutine const_j_rotation()
      real(wp) :: guess, rsm, wwsm, mum, sgp
      real(wp), parameter :: tolerance = 1.e-5_wp
      guess = Omega_e * 0.8e0_wp
      call find_omege_e(guess, r_e_new, rho_equator_h, gama_equator_h, &
                      ww_equator_h, rho_pole_h, gama_pole_h, tolerance, Omega_e, diff_rotation_const_j)

      term_in_Omega_h = abs(Omega_e - ww_equator_h) * exp(-2.0_wp * r_e_new**2 * rho_equator_h)
      Omega_c = Omega_e + term_in_Omega_h / (1.0_wp - term_in_Omega_h * abs(Omega_e - ww_equator_h)) / A_diff**2

      Omg(1,:) = Omega_c
      Omg(1:2*SDIV/3,MDIV) = Omega_c
      do s = 2, 2*SDIV/3
        do m = 1, MDIV-1
          rsm = rho(s,m)
          wwsm = ww(s,m)
          mum = mu(m)
          sgp = s_gp(s)
          call zbrent_rot(Omg(s-1,m) * 8.e-1_wp, r_e_new, rsm, wwsm, sgp, mum, 1.e-5_wp, omg(s,m), rotation_law_const_j)
          F_j(s,m) = (omg(s,m) - wwsm) * sgp**2 * (1.0_wp - mum**2) &
                / ((1.0_wp - sgp)**2 * exp(2.0_wp * r_e_new**2 * rsm) - (omg(s,m) - wwsm)**2 * sgp**2 * (1.0_wp - mum**2))
        end do
      end do
    end subroutine const_j_rotation
    subroutine uryu_rotation()
      real(wp) :: diff_Fmax, guess, Fa, rsm, wwsm, sgp, mum, omg_max_h
      real(wp) :: exp_term_eq
      real(wp), dimension(SDIV) :: omg_mu_0
      integer :: imax
      real(wp), parameter :: tolerance = 1.e-5_wp
      exp_term_eq = exp(2.0_wp * re2 * rho_equator_h)

      diff_Fmax = 1.0_wp
      Fmax_h    = 1.e-2_wp ! Empirial guess; not sure why it works well
      do while(abs(diff_Fmax) > 1.e-7_wp)
        guess = Omega_e
        call find_omege_e(guess, r_e_new,rho_equator_h,gama_equator_h,ww_equator_h, &
                        rho_pole_h,gama_pole_h, tolerance, Fa, diff_rotation_uryu)
        Omega_e = fa
        F_equator_h  = (Omega_e - ww_equator_h) / ( exp_term_eq - (Omega_e-ww_equator_h)**2 )
        if ( F_equator_h < 0.0_wp ) stop "negative F_equator_h; L120 in uryu"

        Omega_c = Omega_e / lambda2
        omg_mu_0(1) = Omega_c
        omg_max_h = Omega_c
        mum = 0.0_wp
        do s = 2, (SDIV-1)/2
            rsm = rho(s,1)
            wwsm= ww (s,1)
            sgp = s_gp(s)
            guess  = omg_mu_0(s-1)
            call zbrent_rot( guess, r_e_new, rsm, wwsm, sgp, mum, 1.e-5_wp, omg_mu_0(s), rotation_law_uryu)
            !write(*,"(es15.6)",advance='no') omg_mu_0(s)
            if (omg_mu_0(s) > omg_max_h) then 
              omg_max_h = omg_mu_0(s)!; write(*,"(A)",advance='no') " <----"
            else
              exit
            end if
            !write(*,*)" "
        enddo
        !stop 803
        diff_Fmax = ( lambda1 - omg_max_h / Omega_c )
        Fmax_h    = Fmax_h - diff_Fmax * 1.e-2_wp
      enddo

      Omg(1,:) = Omega_c
      Omg(1:3*SDIV/4,MDIV) = Omega_c
      do s = 2, SDIV*3/4
        do m = 1, MDIV-1
          rsm = rho(s,m)
          wwsm= ww (s,m)
          mum = mu(m)
          sgp = s_gp(s)
          guess  = Omg(s-1,m)
          call zbrent_rot( guess, r_e_new, rsm, wwsm, sgp, mum, 1.e-5_wp, omg(s,m), rotation_law_uryu)
          F_j(s,m) = (omg(s,m) - wwsm) * sgp**2 * (1.0_wp - mum**2) &
                / ((1.0_wp - sgp)**2 * exp(2.0_wp * re2 * rsm) - (omg(s,m) - wwsm)**2 * sgp**2 * (1.0_wp - mum**2))
        enddo
      enddo
    end subroutine uryu_rotation
  end subroutine update_angular_velocity

  subroutine update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, &
                                     sgp_term_2d_cache_arg, sin_theta_2d_cache_arg, sgp_2d_cache_arg)
    use para_mod, only : wp, SDIV, MDIV, r_ratio, h_center, &
                         rho, gama, alpha, sphi, velocity_sq, enthalpy, pressure, energy, &
                         Omg, ww, Omega_c, omg, F_j, &
                         s_gp, has_scalar, B_coup, A_diff, solver_type, enthalpy_min, s_e
    use rotation_law_mod, only: intF
    use toolkit_mod, only : interp_log_h_to_p, interp_log_p_to_e
    real(wp), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(wp), intent(in) :: sgp_term_2d_cache_arg(:,:), sin_theta_2d_cache_arg(:,:), sgp_2d_cache_arg(:,:)
    integer :: s, m
    real(wp) :: re2, log_p_val, log_e_val, pp, ee

    re2 = r_e_new**2

    if (abs(r_ratio - 1.0_wp) < epsilon(r_ratio)) then
      velocity_sq = 0.0_wp
      enthalpy = enthalpy_min + 0.5e0_wp * re2 * &
               ( gama_pole_h + rho_pole_h - gama - rho + ( sphi**2 - sphi_pole_h**2 ) * B_coup / 2.0_wp )
    else
      velocity_sq = ((Omg - ww) * sgp_term_2d_cache_arg * sin_theta_2d_cache_arg * exp(-rho * re2))**2
      where (velocity_sq > 1.0_wp) velocity_sq = 0.0_wp
      enthalpy = enthalpy_min + 0.5e0_wp * ( &
            re2 * ( gama_pole_h + rho_pole_h - gama - rho + ( sphi**2 - sphi_pole_h**2 ) * B_coup / 2.0_wp ) &
            - log( max(1.e-300_wp, 1.0_wp-velocity_sq) )  )
    end if

    if ( trim(solver_type) == "const_j" ) then
      enthalpy = enthalpy + 0.5e0_wp * A_diff**2 * (Omg - Omega_c)**2
    elseif ( trim(solver_type) == "uryu" ) then
      do m = 1, MDIV
        do s = 1, SDIV
          enthalpy(s,m) = enthalpy(s,m) - intF(omg(s,m), F_j(s,m))
        end do
      end do
    endif

    do m = 1, MDIV
      do s = 1, SDIV
        associate( h => enthalpy(s,m), p => pressure(s,m), &
                   e => energy(s,m),   g => sgp_2d_cache_arg(s,m) )
        if (h > enthalpy_min .and. g <= s_e) then
          ! Direct scalar calls are usually the fastest path for the CPU
          p = exp(interp_log_h_to_p(log(h)))
          e = exp(interp_log_p_to_e(log(p)))
        else
          h = enthalpy_min
          p = 0.0_wp
          e = 0.0_wp
        end if
        end associate
      end do
    end do

    rho   = rho   * re2
    gama  = gama  * re2
    alpha = alpha * re2
    sphi  = sphi  * r_e_new
  end subroutine update_eos_and_velocity

end module spin_updates
