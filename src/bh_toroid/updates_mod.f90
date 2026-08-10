module updates_mod
  use precision_mod, only: wp
  use validation_mod, only: validation_result, validation_ok, validation_error, &
      is_finite, VALID_OK, VALID_BAD_RADIAL_ORDER, VALID_BAD_ROTATION, VALID_BAD_SCALE
  implicit none
  private

  type, public :: bh_toroid_equatorial_point
    real(wp) :: rhat = 0.0_wp
    real(wp) :: nu_hat = 0.0_wp
    real(wp) :: gamma_hat = 0.0_wp
    real(wp) :: omega_hat = 0.0_wp
  end type bh_toroid_equatorial_point

  type, public :: bh_toroid_update_constants
    real(wp) :: r_out = 1.0_wp
    real(wp) :: rotation_A = 1.0_wp
    real(wp) :: poly_k = 1.0_wp
    real(wp) :: poly_n = 3.0_wp
    real(wp) :: omega_c = 0.0_wp
    real(wp) :: bernoulli_c = 0.0_wp
    real(wp) :: surface_mismatch = 0.0_wp
  end type bh_toroid_update_constants

  public :: validate_update_inputs, solve_hydro_rotation_constants
  public :: rotation_omega, rotation_law_residual, toroid_velocity
  public :: integrated_euler_enthalpy, euler_boundary_residual
  public :: checked_rotation_omega, checked_integrated_euler_enthalpy
  public :: polytropic_energy_density, update_hydro_rotation_fields
  public :: finite_equatorial_point

contains

  pure function validate_update_inputs(h0_hat, rin_hat, r_out, rotation_A, poly_k, poly_n) result(res)
    real(wp), intent(in) :: h0_hat, rin_hat, r_out, rotation_A, poly_k, poly_n
    type(validation_result) :: res

    if (.not. all(is_finite([h0_hat, rin_hat, r_out, rotation_A, poly_k, poly_n]))) then
      res = validation_error(VALID_BAD_SCALE, "BH toroid update inputs must be finite")
    else if (h0_hat <= 0.0_wp .or. h0_hat >= 1.0_wp .or. rin_hat <= h0_hat .or. rin_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_RADIAL_ORDER, "update domain must satisfy 0 < h0_hat < rin_hat < 1")
    else if (r_out <= 0.0_wp .or. poly_k <= 0.0_wp .or. poly_n <= 0.0_wp) then
      res = validation_error(VALID_BAD_SCALE, "r_out, poly_k, and poly_n must be positive")
    else if (rotation_A <= 0.0_wp) then
      res = validation_error(VALID_BAD_ROTATION, "rotation_A must be positive")
    else
      res = validation_ok()
    end if
  end function validate_update_inputs

  pure function solve_hydro_rotation_constants(h_point, s_point, t_point, h0_hat, rin_hat, &
      r_out, rotation_A, omega_h, poly_k, poly_n) result(constants)
    type(bh_toroid_equatorial_point), intent(in) :: h_point, s_point, t_point
    real(wp), intent(in) :: h0_hat, rin_hat, r_out, rotation_A, omega_h, poly_k, poly_n
    type(bh_toroid_update_constants) :: constants
    type(validation_result) :: validation
    real(wp) :: omega_s, omega_t, surface_s, surface_t

    constants%r_out = r_out
    constants%rotation_A = rotation_A
    constants%poly_k = poly_k
    constants%poly_n = poly_n

    validation = validate_update_inputs(h0_hat, rin_hat, r_out, rotation_A, poly_k, poly_n)
    if (validation%status /= VALID_OK .or. .not. finite_equatorial_point(h_point) .or. &
        .not. finite_equatorial_point(s_point) .or. .not. finite_equatorial_point(t_point)) then
      constants%surface_mismatch = huge(1.0_wp)
      return
    end if

    constants%omega_c = omega_c_from_horizon(h_point, r_out, rotation_A, omega_h)
    omega_s = rotation_omega(constants, s_point%rhat, 1.0_wp, &
        s_point%gamma_hat, s_point%nu_hat, s_point%omega_hat)
    omega_t = rotation_omega(constants, t_point%rhat, 1.0_wp, &
        t_point%gamma_hat, t_point%nu_hat, t_point%omega_hat)
    surface_s = euler_surface_constant(constants, omega_s, s_point%rhat, 1.0_wp, &
        s_point%gamma_hat, s_point%nu_hat, s_point%omega_hat)
    surface_t = euler_surface_constant(constants, omega_t, t_point%rhat, 1.0_wp, &
        t_point%gamma_hat, t_point%nu_hat, t_point%omega_hat)
    constants%bernoulli_c = 0.5_wp * (surface_s + surface_t)
    constants%surface_mismatch = surface_t - surface_s
  end function solve_hydro_rotation_constants

  ! Nishida & Eriguchi (1994), eqs. (3.10)-(3.14): hatted metric potentials
  ! scaled by r_out. Eq. (3.16) is solved for q = Omega - r_out**2 * omega_hat
  ! by bisection inside |v| < 1.
  pure function rotation_omega(constants, rhat, sin_theta, gamma_hat, nu_hat, omega_hat) result(omega)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp) :: omega
    real(wp) :: low, high, mid, fmid, radial_factor
    integer :: iter

    radial_factor = velocity_radius_factor(constants%r_out, rhat, sin_theta, gamma_hat, nu_hat)
    if (radial_factor <= 0.0_wp .or. .not. is_finite(radial_factor)) then
      omega = constants%r_out**2 * omega_hat
      return
    end if

    low = -safe_velocity_limit(radial_factor)
    high = safe_velocity_limit(radial_factor)
    mid = 0.0_wp
    do iter = 1, 128
      mid = 0.5_wp * (low + high)
      fmid = rotation_law_q_residual(constants, mid, omega_hat, radial_factor)
      if (fmid > 0.0_wp) then
        low = mid
      else
        high = mid
      end if
    end do
    omega = constants%r_out**2 * omega_hat + mid
  end function rotation_omega

  function checked_rotation_omega(constants, rhat, sin_theta, gamma_hat, nu_hat, omega_hat, omega) result(res)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp), intent(out) :: omega
    type(validation_result) :: res
    real(wp) :: factor

    omega = 0.0_wp
    if (.not. all(is_finite([rhat, sin_theta, gamma_hat, nu_hat, omega_hat]))) then
      res = validation_error(VALID_BAD_SCALE, "metric and rotation inputs must be finite")
      return
    end if
    factor = velocity_radius_factor(constants%r_out, rhat, sin_theta, gamma_hat, nu_hat)
    if (constants%rotation_A <= 0.0_wp .or. constants%r_out <= 0.0_wp .or. factor <= 0.0_wp .or. &
        .not. is_finite(factor)) then
      res = validation_error(VALID_BAD_ROTATION, "rotation law denominator is invalid")
      return
    end if

    omega = rotation_omega(constants, rhat, sin_theta, gamma_hat, nu_hat, omega_hat)
    if (.not. is_finite(omega) .or. abs(toroid_velocity(omega, constants%r_out, rhat, sin_theta, &
        gamma_hat, nu_hat, omega_hat)) >= 1.0_wp) then
      omega = 0.0_wp
      res = validation_error(VALID_BAD_ROTATION, "rotation law failed to find subluminal Omega")
    else
      res = validation_ok()
    end if
  end function checked_rotation_omega

  pure function rotation_law_residual(constants, omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat) result(residual)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp) :: residual, q, radial_factor

    radial_factor = velocity_radius_factor(constants%r_out, rhat, sin_theta, gamma_hat, nu_hat)
    q = omega - constants%r_out**2 * omega_hat
    residual = rotation_law_q_residual(constants, q, omega_hat, radial_factor)
  end function rotation_law_residual

  pure function toroid_velocity(omega, r_out, rhat, sin_theta, gamma_hat, nu_hat, omega_hat) result(velocity)
    real(wp), intent(in) :: omega, r_out, rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp) :: velocity

    velocity = (omega - r_out**2 * omega_hat) * velocity_radius_factor(r_out, rhat, sin_theta, gamma_hat, nu_hat)
  end function toroid_velocity

  ! Eq. (3.15) rearranged to return (1+N) ln(K e**(1/N) + 1). Positive values
  ! map to a polytropic energy density; non-positive values are outside matter.
  pure function integrated_euler_enthalpy(constants, omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat) result(enthalpy_term)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp) :: enthalpy_term, velocity, v2

    velocity = toroid_velocity(omega, constants%r_out, rhat, sin_theta, gamma_hat, nu_hat, omega_hat)
    v2 = velocity * velocity
    if (v2 >= 1.0_wp) then
      enthalpy_term = -huge(1.0_wp)
    else
      enthalpy_term = constants%bernoulli_c - constants%r_out**2 * nu_hat &
          - 0.5_wp * log(1.0_wp - v2) &
          + 0.5_wp * constants%rotation_A**2 * (omega - constants%omega_c)**2
    end if
  end function integrated_euler_enthalpy

  function checked_integrated_euler_enthalpy(constants, omega, rhat, sin_theta, gamma_hat, nu_hat, &
      omega_hat, enthalpy_term) result(res)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp), intent(out) :: enthalpy_term
    type(validation_result) :: res
    real(wp) :: velocity

    enthalpy_term = 0.0_wp
    if (.not. all(is_finite([omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat]))) then
      res = validation_error(VALID_BAD_SCALE, "Euler inputs must be finite")
      return
    end if
    velocity = toroid_velocity(omega, constants%r_out, rhat, sin_theta, gamma_hat, nu_hat, omega_hat)
    if (.not. is_finite(velocity) .or. abs(velocity) >= 1.0_wp) then
      res = validation_error(VALID_BAD_ROTATION, "Euler velocity must be finite and subluminal")
      return
    end if

    enthalpy_term = integrated_euler_enthalpy(constants, omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat)
    if (.not. is_finite(enthalpy_term)) then
      enthalpy_term = 0.0_wp
      res = validation_error(VALID_BAD_ROTATION, "Euler enthalpy term is not finite")
    else
      res = validation_ok()
    end if
  end function checked_integrated_euler_enthalpy

  ! Eq. (3.15) residual. Surface points S and T use energy_density = 0,
  ! so this directly checks the boundary constant C.
  pure function euler_boundary_residual(constants, energy_density, omega, rhat, sin_theta, &
      gamma_hat, nu_hat, omega_hat) result(residual)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: energy_density, omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp) :: residual, velocity, v2

    velocity = toroid_velocity(omega, constants%r_out, rhat, sin_theta, gamma_hat, nu_hat, omega_hat)
    v2 = velocity * velocity
    if (v2 >= 1.0_wp .or. energy_density < 0.0_wp .or. .not. is_finite(v2)) then
      residual = huge(1.0_wp)
      return
    end if
    residual = polytropic_enthalpy_log(energy_density, constants%poly_k, constants%poly_n) + &
        euler_surface_constant(constants, omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat) - &
        constants%bernoulli_c
  end function euler_boundary_residual

  pure function polytropic_energy_density(enthalpy_term, constants) result(energy_density)
    real(wp), intent(in) :: enthalpy_term
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp) :: energy_density, base

    if (enthalpy_term <= 0.0_wp) then
      energy_density = 0.0_wp
    else
      base = (exp(enthalpy_term / (1.0_wp + constants%poly_n)) - 1.0_wp) / constants%poly_k
      energy_density = max(0.0_wp, base)**constants%poly_n
    end if
  end function polytropic_energy_density

  function update_hydro_rotation_fields(constants, rin_hat, rhat, sin_theta, gamma_hat, nu_hat, omega_hat, &
      omega, velocity, enthalpy_term, energy_density) result(res)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: rin_hat
    real(wp), intent(in) :: rhat(:), sin_theta(:)
    real(wp), intent(in) :: gamma_hat(:,:), nu_hat(:,:), omega_hat(:,:)
    real(wp), intent(out) :: omega(:,:), velocity(:,:), enthalpy_term(:,:), energy_density(:,:)
    type(validation_result) :: res
    integer :: i, j

    omega = 0.0_wp
    velocity = 0.0_wp
    enthalpy_term = 0.0_wp
    energy_density = 0.0_wp

    if (size(gamma_hat, 1) /= size(rhat) .or. size(gamma_hat, 2) /= size(sin_theta) .or. &
        any(shape(nu_hat) /= shape(gamma_hat)) .or. any(shape(omega_hat) /= shape(gamma_hat)) .or. &
        any(shape(omega) /= shape(gamma_hat)) .or. any(shape(velocity) /= shape(gamma_hat)) .or. &
        any(shape(enthalpy_term) /= shape(gamma_hat)) .or. any(shape(energy_density) /= shape(gamma_hat))) then
      res = validation_error(VALID_BAD_SCALE, "field arrays must share shape (size(rhat), size(sin_theta))")
      return
    end if

    if (.not. all(is_finite(rhat)) .or. .not. all(is_finite(sin_theta)) .or. &
        .not. all(is_finite(gamma_hat)) .or. .not. all(is_finite(nu_hat)) .or. &
        .not. all(is_finite(omega_hat))) then
      res = validation_error(VALID_BAD_SCALE, "field inputs must be finite")
      return
    end if

    do j = 1, size(sin_theta)
      do i = 1, size(rhat)
        if (rhat(i) < rin_hat .or. rhat(i) > 1.0_wp) cycle
        res = checked_rotation_omega(constants, rhat(i), sin_theta(j), gamma_hat(i,j), nu_hat(i,j), &
            omega_hat(i,j), omega(i,j))
        if (res%status /= VALID_OK) exit
        velocity(i,j) = toroid_velocity(omega(i,j), constants%r_out, rhat(i), sin_theta(j), &
            gamma_hat(i,j), nu_hat(i,j), omega_hat(i,j))
        res = checked_integrated_euler_enthalpy(constants, omega(i,j), rhat(i), sin_theta(j), &
            gamma_hat(i,j), nu_hat(i,j), omega_hat(i,j), enthalpy_term(i,j))
        if (res%status /= VALID_OK) exit
        enthalpy_term(i,j) = max(0.0_wp, enthalpy_term(i,j))
        energy_density(i,j) = polytropic_energy_density(enthalpy_term(i,j), constants)
      end do
      if (res%status /= VALID_OK) exit
    end do

    if (res%status /= VALID_OK) then
      omega = 0.0_wp
      velocity = 0.0_wp
      enthalpy_term = 0.0_wp
      energy_density = 0.0_wp
    end if
  end function update_hydro_rotation_fields

  pure function omega_c_from_horizon(point, r_out, rotation_A, omega_h) result(omega_c)
    type(bh_toroid_equatorial_point), intent(in) :: point
    real(wp), intent(in) :: r_out, rotation_A, omega_h
    real(wp) :: omega_c, velocity, q, v2, denom

    velocity = toroid_velocity(omega_h, r_out, point%rhat, 1.0_wp, point%gamma_hat, point%nu_hat, point%omega_hat)
    q = omega_h - r_out**2 * point%omega_hat
    v2 = velocity * velocity
    denom = (1.0_wp - v2) * q
    if (rotation_A <= 0.0_wp .or. abs(denom) <= tiny(1.0_wp)) then
      omega_c = omega_h
    else
      omega_c = omega_h + v2 / (denom * rotation_A**2)
    end if
  end function omega_c_from_horizon

  pure function euler_surface_constant(constants, omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat) result(value)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: omega, rhat, sin_theta, gamma_hat, nu_hat, omega_hat
    real(wp) :: value, velocity, v2

    velocity = toroid_velocity(omega, constants%r_out, rhat, sin_theta, gamma_hat, nu_hat, omega_hat)
    v2 = velocity * velocity
    if (v2 >= 1.0_wp .or. .not. is_finite(v2)) then
      value = huge(1.0_wp)
    else
      value = constants%r_out**2 * nu_hat + 0.5_wp * log(1.0_wp - v2) &
          - 0.5_wp * constants%rotation_A**2 * (omega - constants%omega_c)**2
    end if
  end function euler_surface_constant

  pure function rotation_law_q_residual(constants, q, omega_hat, radial_factor) result(residual)
    type(bh_toroid_update_constants), intent(in) :: constants
    real(wp), intent(in) :: q, omega_hat, radial_factor
    real(wp) :: residual, v2, denom, omega_frame

    omega_frame = constants%r_out**2 * omega_hat
    v2 = (q * radial_factor)**2
    denom = (1.0_wp - v2) * q
    if (constants%rotation_A <= 0.0_wp .or. radial_factor <= 0.0_wp .or. v2 >= 1.0_wp) then
      residual = huge(1.0_wp)
    else if (abs(q) <= tiny(1.0_wp)) then
      residual = constants%rotation_A**2 * (constants%omega_c - omega_frame)
    else
      residual = constants%rotation_A**2 * (constants%omega_c - omega_frame - q) - v2 / denom
    end if
  end function rotation_law_q_residual

  pure function velocity_radius_factor(r_out, rhat, sin_theta, gamma_hat, nu_hat) result(factor)
    real(wp), intent(in) :: r_out, rhat, sin_theta, gamma_hat, nu_hat
    real(wp) :: factor

    factor = r_out * rhat * sin_theta * exp(r_out**2 * (gamma_hat - 2.0_wp * nu_hat))
  end function velocity_radius_factor

  pure function safe_velocity_limit(radial_factor) result(limit)
    real(wp), intent(in) :: radial_factor
    real(wp) :: limit

    limit = (1.0_wp - 128.0_wp * epsilon(1.0_wp)) / radial_factor
  end function safe_velocity_limit

  pure function polytropic_enthalpy_log(energy_density, poly_k, poly_n) result(value)
    real(wp), intent(in) :: energy_density, poly_k, poly_n
    real(wp) :: value

    if (energy_density <= 0.0_wp) then
      value = 0.0_wp
    else
      value = (1.0_wp + poly_n) * log(poly_k * energy_density**(1.0_wp / poly_n) + 1.0_wp)
    end if
  end function polytropic_enthalpy_log

  pure function finite_equatorial_point(point) result(ok)
    type(bh_toroid_equatorial_point), intent(in) :: point
    logical :: ok

    ok = all(is_finite([point%rhat, point%nu_hat, point%gamma_hat, point%omega_hat]))
  end function finite_equatorial_point

end module updates_mod
