program test_rotation_law
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_positive_inf, &
                                           ieee_quiet_nan, ieee_value
  use precision_mod, only: wp
  use para_mod, only: A_diff, F_equator_h, Fmax_h, Omega_c, lambda1, lambda2, uryu_p, uryu_q
  use rotation_law_mod, only: cached_aa, cached_bb, ctx_exp2re2rho, &
                              cache_uryu_ab, diff_rotation_const_j, diff_rotation_uryu, intF, &
                              rotation_law_const_j, rotation_law_uryu, set_shoot_context, &
                              uryu_coefficients, ctx_rho_equator_h, ctx_gama_equator_h, &
                              ctx_ww_equator_h, ctx_rho_pole_h, ctx_gama_pole_h
  use brent_mod, only: find_omega_e, zbrent_rot, find_omega_e_admissible, zbrent_rot_admissible
  use test_utils
  implicit none

  real(wp) :: aa, bb, fx, expected, root, threshold
  logical :: valid, found
  real(wp), parameter :: tol = 1.0e-10_wp

  A_diff = 0.5_wp
  Omega_c = 0.6_wp
  call diff_rotation_const_j(0.2_wp, fx, 0.8_wp, 0.1_wp, 0.2_wp, 0.01_wp, 0.05_wp, 0.12_wp)
  expected = 0.8_wp**2*(0.2_wp+0.1_wp-0.12_wp-0.05_wp) &
      + log(1._wp-((0.2_wp-0.01_wp)*exp(-0.8_wp**2*0.1_wp))**2)
  expected = expected*(1._wp-((0.2_wp-0.01_wp)*exp(-0.8_wp**2*0.1_wp))**2)**2 &
      - ((0.2_wp-0.01_wp)*exp(-2._wp*0.8_wp**2*0.1_wp)/A_diff)**2
  call assert_near("const-j equatorial residual", expected, fx, tol)

  ctx_exp2re2rho = 1.2_wp
  call rotation_law_const_j(0.3_wp, fx, 1._wp, 0._wp, 0.02_wp, 0.25_wp, 0.4_wp)
  expected = (Omega_c-0.3_wp)*((1._wp-0.25_wp)**2 &
      -(0.25_wp*(0.3_wp-0.02_wp))**2/ctx_exp2re2rho*(1._wp-0.4_wp**2)) &
      -(0.3_wp-0.02_wp)*0.25_wp**2*(1._wp-0.4_wp**2)/ctx_exp2re2rho/A_diff**2
  call assert_near("const-j local residual", expected, fx, tol)

  lambda1 = 1.5_wp
  lambda2 = 0.3_wp
  uryu_p = 1
  uryu_q = 3
  call uryu_coefficients(2.222222222222222_wp, 1._wp, aa, bb, valid)
  call assert_true("admissible Uryu coefficients", valid)
  call assert_true("positive finite Uryu aa", aa > 0._wp .and. ieee_is_finite(aa))
  call assert_true("positive finite Uryu bb", bb > 0._wp .and. ieee_is_finite(bb))
  call uryu_coefficients(0.1_wp, 1._wp, aa, bb, valid)
  call assert_true("inadmissible Uryu coefficients rejected", .not. valid)
  call uryu_coefficients(0._wp, 1._wp, aa, bb, valid)
  call assert_true("zero Uryu momentum rejected", .not. valid)
  call uryu_coefficients(-1._wp, 1._wp, aa, bb, valid)
  call assert_true("negative Uryu momentum rejected", .not. valid)
  call uryu_coefficients(ieee_value(0._wp, ieee_quiet_nan), 1._wp, aa, bb, valid)
  call assert_true("NaN Uryu momentum rejected", .not. valid)
  call uryu_coefficients(ieee_value(0._wp, ieee_positive_inf), 1._wp, aa, bb, valid)
  call assert_true("infinite Uryu momentum rejected", .not. valid)
  threshold = (lambda1/lambda2)**(1._wp/real(uryu_q, wp))
  call uryu_coefficients(threshold, 1._wp, aa, bb, valid)
  call assert_true("Uryu coefficient boundary rejected", .not. valid)
  call uryu_coefficients(threshold*(1._wp+sqrt(epsilon(threshold))), 1._wp, aa, bb, valid)
  call assert_true("Uryu coefficients valid above boundary", valid)

  F_equator_h = 2.222222222222222_wp
  Fmax_h = 1._wp
  call cache_uryu_ab(valid)
  call assert_true("valid Uryu coefficients cached", valid)
  call assert_true("cached Uryu coefficients positive", cached_aa > 0._wp .and. cached_bb > 0._wp)
  F_equator_h = 0.1_wp
  call cache_uryu_ab(valid)
  call assert_true("invalid Uryu cache request refused", .not. valid)
  call assert_near("invalid cache clears aa", 0._wp, cached_aa, 0._wp)
  call assert_near("invalid cache clears bb", 0._wp, cached_bb, 0._wp)
  call rotation_law_uryu(0.5_wp, fx, 1._wp, 0._wp, 0._wp, 0.2_wp, 0.3_wp)
  call assert_true("invalid Uryu cache marks local residual invalid", ieee_is_nan(fx))

  F_equator_h = 2.222222222222222_wp
  call cache_uryu_ab(valid)
  Omega_c = 0.8_wp/lambda2
  call diff_rotation_uryu(0.8_wp, fx, 1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp)
  call assert_true("admissible Uryu equatorial residual finite", ieee_is_finite(fx))
  call diff_rotation_uryu(0.1_wp, fx, 1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp)
  call assert_true("inadmissible Uryu equatorial residual marked invalid", ieee_is_nan(fx))

  ctx_exp2re2rho = 1._wp
  call rotation_law_uryu(0.5_wp, fx, 1._wp, 0._wp, 0._wp, 0.2_wp, 0.3_wp)
  call assert_true("Uryu local residual finite", ieee_is_finite(fx))
  call assert_near("zero Uryu integral", 0._wp, intF(0._wp, 0._wp), 0._wp)
  call assert_true("nonzero Uryu integral finite", ieee_is_finite(intF(0.5_wp, 0.2_wp)))

  call set_shoot_context(0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp, 0.5_wp)
  call assert_near("shoot context equatorial rho", 0.1_wp, ctx_rho_equator_h, 0._wp)
  call assert_near("shoot context equatorial gamma", 0.2_wp, ctx_gama_equator_h, 0._wp)
  call assert_near("shoot context equatorial omega", 0.3_wp, ctx_ww_equator_h, 0._wp)
  call assert_near("shoot context polar rho", 0.4_wp, ctx_rho_pole_h, 0._wp)
  call assert_near("shoot context polar gamma", 0.5_wp, ctx_gama_pole_h, 0._wp)

  call find_omega_e_admissible(0.2_wp, 1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
                    1.0e-10_wp, root, bounded_residual, found)
  call assert_true("bounded Omega_e root found", found)
  call assert_near("bounded Omega_e root", 0.8_wp, root, 1.0e-8_wp)
  call find_omega_e_admissible(0.2_wp, 1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
                    1.0e-10_wp, root, no_root_residual, found)
  call assert_true("no-bracket Omega_e search refused", .not. found)
  call assert_true("no-bracket Omega_e result invalid", ieee_is_nan(root))
  call find_omega_e_admissible(0.5_wp, 1._wp, 0._wp, 0._wp, 0._wp, 0._wp, 0._wp, &
                    1.0e-10_wp, root, disconnected_residual, found)
  call assert_true("disconnected finite domains do not bracket", .not. found)

  A_diff = 0.5_wp
  call find_omega_e(0.2_wp, 1._wp, 0._wp, 0.3_wp, 0._wp, 0._wp, 0._wp, &
                    1.0e-10_wp, root, diff_rotation_const_j)
  call diff_rotation_const_j(root, fx, 1._wp, 0._wp, 0.3_wp, 0._wp, 0._wp, 0._wp)
  call assert_near("const-j callback Omega_e residual", 0._wp, fx, 1.0e-8_wp)
  Fmax_h = 0.01_wp
  call find_omega_e_admissible(0.2_wp, 1._wp, 0._wp, 0.1_wp, 0._wp, 0._wp, 0._wp, &
                    1.0e-10_wp, root, diff_rotation_uryu, found)
  call assert_true("Uryu callback Omega_e root found", found)
  call diff_rotation_uryu(root, fx, 1._wp, 0._wp, 0.1_wp, 0._wp, 0._wp, 0._wp)
  call assert_near("Uryu callback Omega_e residual", 0._wp, fx, 1.0e-8_wp)

  call zbrent_rot_admissible(0.2_wp, 1._wp, 0._wp, 0._wp, 0.25_wp, 0._wp, &
                  1.0e-10_wp, root, local_residual, found)
  call assert_true("bounded local rotation root found", found)
  call assert_near("bounded local rotation root", 0.6_wp, root, 1.0e-8_wp)
  call zbrent_rot_admissible(0.2_wp, 1._wp, 0._wp, 0._wp, 0.25_wp, 0._wp, &
                  1.0e-10_wp, root, local_no_root_residual, found)
  call assert_true("no-bracket local rotation search refused", .not. found)
  call assert_true("no-bracket local rotation result invalid", ieee_is_nan(root))

  ctx_exp2re2rho = 1._wp
  Omega_c = 0.6_wp
  call zbrent_rot(0.2_wp, 1._wp, 0._wp, 0._wp, 0.25_wp, 0._wp, &
                  1.0e-10_wp, root, rotation_law_const_j)
  call rotation_law_const_j(root, fx, 1._wp, 0._wp, 0._wp, 0.25_wp, 0._wp)
  call assert_near("const-j callback local residual", 0._wp, fx, 1.0e-8_wp)
  Fmax_h = 1._wp
  F_equator_h = 2.222222222222222_wp
  call cache_uryu_ab(valid)
  Omega_c = 0.8_wp/lambda2
  call zbrent_rot_admissible(2.5_wp, 1._wp, 0._wp, 0._wp, 0.25_wp, 0._wp, &
                  1.0e-10_wp, root, rotation_law_uryu, found)
  call assert_true("Uryu callback local root found", found)
  call rotation_law_uryu(root, fx, 1._wp, 0._wp, 0._wp, 0.25_wp, 0._wp)
  call assert_near("Uryu callback local residual", 0._wp, fx, 1.0e-8_wp)

  call test_summary("rotation_law_mod")

contains

  subroutine bounded_residual(x, residual, re, rho_h, g_h, w_h, rho_p, g_p)
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out) :: residual
    if (x <= 0.75_wp) then
      residual = ieee_value(0._wp, ieee_quiet_nan)
    else
      residual = x - 0.8_wp
    end if
  end subroutine bounded_residual

  subroutine no_root_residual(x, residual, re, rho_h, g_h, w_h, rho_p, g_p)
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out) :: residual
    residual = x + 1._wp
  end subroutine no_root_residual

  subroutine disconnected_residual(x, residual, re, rho_h, g_h, w_h, rho_p, g_p)
    real(wp), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(wp), intent(out) :: residual
    if (x < 0.4_wp) then
      residual = -1._wp
    else if (x > 0.6_wp) then
      residual = 1._wp
    else
      residual = ieee_value(0._wp, ieee_quiet_nan)
    end if
  end subroutine disconnected_residual

  subroutine local_residual(x, residual, re, rho_p, ww_p, sgp, mugp)
    real(wp), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(wp), intent(out) :: residual
    residual = x - 0.6_wp
  end subroutine local_residual

  subroutine local_no_root_residual(x, residual, re, rho_p, ww_p, sgp, mugp)
    real(wp), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(wp), intent(out) :: residual
    residual = x + 1._wp
  end subroutine local_no_root_residual

end program test_rotation_law
