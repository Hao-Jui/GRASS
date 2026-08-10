module rotation_law_mod
  use precision_mod, only: wp
  implicit none
  public
  real(wp), save :: ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h
  real(wp), save :: ctx_rho_pole_h, ctx_gama_pole_h
  real(wp), save :: cached_aa = 0._wp, cached_bb = 0._wp
  real(wp), save :: ctx_exp2re2rho = 1._wp
contains

  subroutine cache_uryu_ab(valid)
    use para_mod, only: wp, F_equator_h, Fmax_h
    logical, intent(out), optional :: valid
    logical :: coefficients_valid

    call uryu_coefficients(F_equator_h, Fmax_h, cached_aa, cached_bb, coefficients_valid)
    if (.not. coefficients_valid) then
      cached_aa = 0._wp
      cached_bb = 0._wp
      if (present(valid)) then
        valid = .false.
        return
      end if
      error stop "cache_uryu_ab: inadmissible Uryu coefficients"
    end if
    if (present(valid)) valid = .true.
  end subroutine cache_uryu_ab

  subroutine diff_rotation_const_j(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    use para_mod, only: wp, A_diff
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out):: fx
    real(wp) :: LHS
    LHS =  re**2 * (g_h + rho_h - g_p - rho_p) + log(1._wp - ((x-w_h)*exp(-re**2*rho_h))**2 )
    fx = LHS * (1._wp - ((x-w_h)*exp(-re**2*rho_h))**2 )**2 &
        - ( (x-w_h) * exp(-2._wp*re**2*rho_h) / A_diff )**2
  end subroutine diff_rotation_const_j

  subroutine rotation_law_const_j(x, fx, re, rho_p, ww_p, sgp, mugp)
    use para_mod, only: wp, A_diff, Omega_c
    real(wp), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(wp), intent(out):: fx
    real(wp) :: sin2m, inv_e2
    sin2m = 1._wp - mugp * mugp
    inv_e2 = 1._wp / ctx_exp2re2rho
    fx = (Omega_c - x) * ( (1._wp-sgp)**2 - ( sgp*(x-ww_p) )**2 * inv_e2 * sin2m ) &
        - ( (x-ww_p) * sgp * sgp * sin2m * inv_e2 ) / A_diff**2
  end subroutine rotation_law_const_j

  subroutine diff_rotation_uryu(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use para_mod, only: wp, lambda2, Fmax_h
    implicit none
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out):: fx
    real(wp) :: F_e, ocre, RHS, aa, bb, delta, exp2rho, log_argument, atanh_argument
    logical :: coefficients_valid

    delta = x - w_h
    exp2rho = exp(2._wp*re**2*rho_h)
    log_argument = 1._wp - (delta*exp(-re**2*rho_h))**2
    if (.not. ieee_is_finite(exp2rho) .or. exp2rho <= delta**2 .or. &
        .not. ieee_is_finite(log_argument) .or. log_argument <= 0._wp .or. &
        .not. ieee_is_finite(lambda2) .or. abs(lambda2) <= tiny(lambda2)) then
      fx = ieee_value(0._wp, ieee_quiet_nan)
      return
    end if
    ocre = x / lambda2
    F_e = delta / (exp2rho - delta**2)

    call uryu_coefficients(F_e, Fmax_h, aa, bb, coefficients_valid)
    if (.not. coefficients_valid) then
      fx = ieee_value(0._wp, ieee_quiet_nan)
      return
    end if

    atanh_argument = aa * F_e * sqrt(2._wp) / (F_e**2 + aa**2)
    if (.not. ieee_is_finite(atanh_argument) .or. abs(atanh_argument) >= 1._wp) then
      fx = ieee_value(0._wp, ieee_quiet_nan)
      return
    end if

    RHS = F_e * x - aa * ocre / 4._wp * &
          (2._wp * aa / bb * atan(F_e**2/aa**2) &
          - sqrt(2._wp) * ( atan(1._wp - F_e*sqrt(2._wp)/aa) - atan(1._wp + F_e*sqrt(2._wp)/aa) ) &
          + sqrt(2._wp) * ATANH(atanh_argument) )

    fx = re**2 * (g_h + rho_h - g_p - rho_p) + log(log_argument) &
        + 2._wp * RHS
    if (.not. ieee_is_finite(fx)) fx = ieee_value(0._wp, ieee_quiet_nan)
  end subroutine diff_rotation_uryu

  subroutine rotation_law_uryu(x, fx, re, rho_p, ww_p, sgp, mugp)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
    use para_mod, only: wp, uryu_p, uryu_q, Omega_c
    implicit none
    real(wp), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(wp), intent(out):: fx
    real(wp) :: Fj, sin2m, denominator

    sin2m = 1._wp - mugp**2
    denominator = ctx_exp2re2rho * (1._wp-sgp)**2 - (x-ww_p)**2 * sgp**2 * sin2m
    if (.not. ieee_is_finite(denominator) .or. denominator <= 0._wp .or. &
        .not. ieee_is_finite(cached_aa) .or. cached_aa <= 0._wp .or. &
        .not. ieee_is_finite(cached_bb) .or. cached_bb <= 0._wp .or. &
        .not. ieee_is_finite(Omega_c) .or. abs(Omega_c) <= tiny(Omega_c) .or. &
        uryu_p <= 0 .or. uryu_q <= 0) then
      fx = ieee_value(0._wp, ieee_quiet_nan)
      return
    end if
    Fj = (x-ww_p) * sgp**2 * sin2m / denominator
    if (.not. ieee_is_finite(Fj) .or. Fj < 0._wp) then
      fx = ieee_value(0._wp, ieee_quiet_nan)
      return
    end if
    fx = x / Omega_c * ( 1._wp + (Fj / cached_aa)**(uryu_p+uryu_q) ) - ( 1._wp + (Fj / cached_bb)**uryu_p )
    if (.not. ieee_is_finite(fx)) fx = ieee_value(0._wp, ieee_quiet_nan)
  end subroutine rotation_law_uryu

  elemental real(wp) function intF(x, F_at_x)
    use para_mod, only: wp, Omega_c
    implicit none
    real(wp), intent(in) :: x, F_at_x

    if ( abs(x) < epsilon(x) .and. abs(F_at_x) < epsilon(F_at_x) ) then
      intF = 0._wp
    else
      intF = F_at_x * x - cached_aa * omega_c / 4._wp * &
          (2._wp * cached_aa / cached_bb * atan(F_at_x**2/cached_aa**2) &
          - sqrt(2._wp) * ( atan(1._wp-F_at_x*sqrt(2._wp)/cached_aa) - atan(1._wp+F_at_x*sqrt(2._wp)/cached_aa) ) &
          + sqrt(2._wp) * ATANH(cached_aa*F_at_x*sqrt(2._wp)/(F_at_x**2 + cached_aa**2)) )
    endif
  end function intF

  pure subroutine uryu_coefficients(F_e, F_m, aa, bb, valid)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use para_mod, only: wp, lambda1, lambda2, uryu_p, uryu_q
    implicit none
    real(wp), intent(in) :: F_e, F_m
    real(wp), intent(out) :: aa, bb
    logical, intent(out) :: valid
    real(wp) :: tmp1, tmp2_aa, tmp2_bb, aa_base, bb_base, minimum_F_e

    aa = 0._wp
    bb = 0._wp
    valid = .false.
    if (.not. ieee_is_finite(F_e) .or. .not. ieee_is_finite(F_m) .or. &
        .not. ieee_is_finite(lambda1) .or. .not. ieee_is_finite(lambda2) .or. &
        F_e <= 0._wp .or. F_m <= 0._wp .or. lambda1 <= 1._wp .or. &
        lambda2 <= 0._wp .or. lambda2 >= 1._wp .or. uryu_p <= 0 .or. uryu_q <= 0) return
    minimum_F_e = F_m*(lambda1/lambda2)**(1._wp/real(uryu_q, wp))
    if (.not. ieee_is_finite(minimum_F_e) .or. F_e <= minimum_F_e) return
    tmp1 = lambda2 * F_e**uryu_q &
        - lambda1 * F_m**uryu_q
    tmp2_aa = (lambda1-1._wp) * F_e**uryu_p &
        - (lambda2-1._wp) * F_m**uryu_p
    tmp2_bb = lambda2 * (lambda1-1._wp) * F_e**(uryu_p+uryu_q) &
        - lambda1 * (lambda2-1._wp) * F_m**(uryu_p+uryu_q)
    if (.not. ieee_is_finite(tmp1) .or. .not. ieee_is_finite(tmp2_aa) .or. &
        .not. ieee_is_finite(tmp2_bb) .or. abs(tmp2_aa) <= tiny(tmp2_aa) .or. &
        abs(tmp2_bb) <= tiny(tmp2_bb)) return

    aa_base = (F_e*F_m)**uryu_p * tmp1 / tmp2_aa
    bb_base = tmp1 / tmp2_bb
    if (.not. ieee_is_finite(aa_base) .or. .not. ieee_is_finite(bb_base) .or. &
        aa_base <= 0._wp .or. bb_base <= 0._wp) return

    aa = aa_base**(1._wp/real(uryu_p+uryu_q, wp))
    bb = F_e * F_m * bb_base**(1._wp/real(uryu_p, wp))
    valid = ieee_is_finite(aa) .and. ieee_is_finite(bb) .and. aa > 0._wp .and. bb > 0._wp
    if (.not. valid) then
      aa = 0._wp
      bb = 0._wp
    end if
  end subroutine uryu_coefficients

  subroutine set_shoot_context(rho_equator_h, gama_equator_h, ww_equator_h, rho_pole_h, gama_pole_h)
    real(wp), intent(in) :: rho_equator_h, gama_equator_h, ww_equator_h, rho_pole_h, gama_pole_h
    ctx_rho_equator_h  = rho_equator_h
    ctx_gama_equator_h = gama_equator_h
    ctx_ww_equator_h   = ww_equator_h
    ctx_rho_pole_h     = rho_pole_h
    ctx_gama_pole_h    = gama_pole_h
  end subroutine set_shoot_context

  subroutine shoot_Fmax(current_Fmax_h, diff_sq)
    use para_mod, only: wp, Fmax_h, Omega_e, Omega_c, r_e, rho, ww, F_equator_h, &
                        SDIV, s_gp, lambda1, lambda2
    use brent_mod, only: find_omega_e_admissible, zbrent_rot_admissible
    real(wp), intent(in)  :: current_Fmax_h
    real(wp), intent(out) :: diff_sq
    real(wp) :: guess, Fa, rsm, wwsm, sgp, mum, omg_max_h, diff_Fmax
    real(wp), dimension(SDIV) :: omg_mu_0
    integer :: imax, s
    real(wp), parameter :: TOLERANCE = 1.e-5

    Fmax_h = current_Fmax_h

    guess = Omega_e
    call find_omega_e_admissible(guess, r_e, ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h, &
                    ctx_rho_pole_h, ctx_gama_pole_h, tolerance, Fa, diff_rotation_uryu)
    Omega_e = Fa
    F_equator_h = (Omega_e - ww(1,1)) / (exp(2._wp*r_e**2*rho(1,1)) - (Omega_e-ww(1,1))**2)
    if (F_equator_h < 0._wp) then
      diff_sq = 1.e3_wp
      return
    endif

    Omega_c = Omega_e / lambda2
    call cache_uryu_ab()
    omg_mu_0(1) = Omega_c
    mum = 0._wp
    do s = 2, SDIV*2/3
      rsm = rho(s,1)
      wwsm= ww (s,1)
      sgp = s_gp(s)
      guess  = omg_mu_0(s-1)
      ctx_exp2re2rho = exp(2._wp * r_e**2 * rsm)
      call zbrent_rot_admissible( guess, r_e, rsm, wwsm, sgp, mum, 1.e-5_wp, omg_mu_0(s), rotation_law_uryu)
    enddo
    imax = maxloc(omg_mu_0(1:SDIV*2/3), 1)
    omg_max_h = omg_mu_0(imax)
    diff_Fmax = (lambda1 - omg_max_h / Omega_c)
    diff_sq = diff_Fmax**2
  end subroutine shoot_Fmax
end module rotation_law_mod
