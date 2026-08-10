program test_bh_toroid_updates
  use precision_mod, only: wp
  use updates_mod
  use validation_mod, only: validation_result, VALID_OK, VALID_BAD_RADIAL_ORDER, &
      VALID_BAD_ROTATION, VALID_BAD_SCALE
  use test_utils
  implicit none

  real(wp), parameter :: tol = 1.e-11_wp
  type(bh_toroid_update_constants) :: constants
  type(bh_toroid_equatorial_point) :: h_point, s_point, t_point
  type(validation_result) :: res
  real(wp) :: omega, velocity, residual, enthalpy_term, energy_density, zero
  real(wp) :: rhat(5), sin_theta(1)
  real(wp) :: gamma_hat(5,1), nu_hat(5,1), omega_hat(5,1)
  real(wp) :: omega_grid(5,1), velocity_grid(5,1), enthalpy_grid(5,1), energy_grid(5,1)

  call assert_status("valid update inputs", &
      validate_update_inputs(0.2_wp, 0.5_wp, 2.0_wp, 3.0_wp, 0.5_wp, 2.0_wp), VALID_OK)
  call assert_status("bad radial order", &
      validate_update_inputs(0.5_wp, 0.5_wp, 2.0_wp, 3.0_wp, 0.5_wp, 2.0_wp), VALID_BAD_RADIAL_ORDER)
  call assert_status("bad scale", &
      validate_update_inputs(0.2_wp, 0.5_wp, 0.0_wp, 3.0_wp, 0.5_wp, 2.0_wp), VALID_BAD_SCALE)
  call assert_status("bad rotation", &
      validate_update_inputs(0.2_wp, 0.5_wp, 2.0_wp, 0.0_wp, 0.5_wp, 2.0_wp), VALID_BAD_ROTATION)

  velocity = toroid_velocity(0.1_wp, 1.0_wp, 0.5_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp)
  call assert_near("velocity eq 3.17", 0.05_wp, velocity, tol)

  constants = bh_toroid_update_constants(r_out=1.0_wp, rotation_A=2.0_wp, poly_k=0.5_wp, poly_n=2.0_wp, &
      omega_c=0.106265664160401_wp, bernoulli_c=0.0_wp, surface_mismatch=0.0_wp)
  omega = rotation_omega(constants, 0.5_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp)
  call assert_near("rotation law omega", 0.1_wp, omega, 1.e-10_wp)
  call assert_status("checked rotation law ok", &
      checked_rotation_omega(constants, 0.5_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, omega), VALID_OK)
  residual = rotation_law_residual(constants, omega, 0.5_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp)
  call assert_near("rotation law residual", 0.0_wp, residual, 1.e-10_wp)
  residual = rotation_law_residual(constants, omega + 0.01_wp, 0.5_wp, 1.0_wp, &
      0.0_wp, 0.0_wp, 0.0_wp)
  call assert_true("rotation law perturbation is detectable", abs(residual) > 1.e-4_wp)
  call assert_status("checked rotation law invalid denominator", &
      checked_rotation_omega(constants, 0.5_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, omega), VALID_BAD_ROTATION)

  constants = bh_toroid_update_constants(r_out=2.0_wp, rotation_A=3.0_wp, poly_k=0.5_wp, poly_n=2.0_wp, &
      omega_c=0.12_wp, bernoulli_c=0.2_wp, surface_mismatch=0.0_wp)
  enthalpy_term = integrated_euler_enthalpy(constants, 0.1_wp, 0.5_wp, 1.0_wp, 0.0_wp, 0.01_wp, 0.0_wp)
  call assert_near("integrated Euler enthalpy", 0.166078976078639_wp, enthalpy_term, 1.e-8_wp)
  call assert_status("checked Euler enthalpy ok", checked_integrated_euler_enthalpy(constants, 0.1_wp, &
      0.5_wp, 1.0_wp, 0.0_wp, 0.01_wp, 0.0_wp, enthalpy_term), VALID_OK)
  energy_density = polytropic_energy_density(enthalpy_term, constants)
  call assert_near("polytropic density inverse", 0.012959853896132_wp, energy_density, 1.e-10_wp)
  residual = euler_boundary_residual(constants, energy_density, 0.1_wp, 0.5_wp, 1.0_wp, 0.0_wp, 0.01_wp, 0.0_wp)
  call assert_near("Euler residual with density", 0.0_wp, residual, 1.e-12_wp)

  h_point = bh_toroid_equatorial_point(rhat=0.2_wp, nu_hat=0.0_wp, gamma_hat=0.0_wp, omega_hat=0.0_wp)
  s_point = bh_toroid_equatorial_point(rhat=0.5_wp, nu_hat=0.02_wp, gamma_hat=0.0_wp, omega_hat=0.0_wp)
  t_point = bh_toroid_equatorial_point(rhat=1.0_wp, nu_hat=0.02_wp, gamma_hat=0.0_wp, omega_hat=0.0_wp)
  constants = solve_hydro_rotation_constants(h_point, s_point, t_point, 0.2_wp, 0.5_wp, &
      2.0_wp, 3.0_wp, 0.0_wp, 0.5_wp, 2.0_wp)
  call assert_near("H-derived omega_c", 0.0_wp, constants%omega_c, tol)
  call assert_near("surface Bernoulli constant", 0.08_wp, constants%bernoulli_c, tol)
  call assert_near("surface mismatch", 0.0_wp, constants%surface_mismatch, tol)
  call assert_near("H rotation residual", 0.0_wp, &
      rotation_law_residual(constants, 0.0_wp, 0.2_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp), tol)
  call assert_near("S boundary residual", 0.0_wp, &
      euler_boundary_residual(constants, 0.0_wp, 0.0_wp, 0.5_wp, 1.0_wp, 0.0_wp, 0.02_wp, 0.0_wp), tol)
  call assert_near("T boundary residual", 0.0_wp, &
      euler_boundary_residual(constants, 0.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.02_wp, 0.0_wp), tol)

  rhat = [0.2_wp, 0.5_wp, 0.75_wp, 1.0_wp, 1.2_wp]
  sin_theta = [1.0_wp]
  gamma_hat = 0.0_wp
  nu_hat(:,1) = [0.0_wp, 0.02_wp, 0.01_wp, 0.02_wp, 0.0_wp]
  omega_hat = 0.0_wp
  res = update_hydro_rotation_fields(constants, 0.5_wp, rhat, sin_theta, gamma_hat, nu_hat, omega_hat, &
      omega_grid, velocity_grid, enthalpy_grid, energy_grid)
  call assert_true("field update valid", res%status == VALID_OK)
  call assert_near("vacuum gap density", 0.0_wp, energy_grid(1,1), tol)
  call assert_near("vacuum gap enthalpy", 0.0_wp, enthalpy_grid(1,1), tol)
  call assert_near("inner surface density", 0.0_wp, energy_grid(2,1), tol)
  call assert_near("inner surface enthalpy", 0.0_wp, enthalpy_grid(2,1), tol)
  call assert_true("interior density positive", energy_grid(3,1) > 0.0_wp)
  call assert_near("outer surface density", 0.0_wp, energy_grid(4,1), tol)
  call assert_near("outer surface enthalpy", 0.0_wp, enthalpy_grid(4,1), tol)
  call assert_near("outside density", 0.0_wp, energy_grid(5,1), tol)
  call assert_near("outside enthalpy", 0.0_wp, enthalpy_grid(5,1), tol)
  call assert_true("field update outputs finite omega", all(abs(omega_grid) < huge(1.0_wp)))
  call assert_true("field update outputs finite velocity", all(abs(velocity_grid) < huge(1.0_wp)))

  zero = 0.0_wp
  gamma_hat(3,1) = zero / zero
  res = update_hydro_rotation_fields(constants, 0.5_wp, rhat, sin_theta, gamma_hat, nu_hat, omega_hat, &
      omega_grid, velocity_grid, enthalpy_grid, energy_grid)
  call assert_true("field update rejects NaN", res%status == VALID_BAD_SCALE)
  call assert_near("NaN rejection clears density", 0.0_wp, sum(energy_grid), tol)

  call test_summary("test_bh_toroid_updates")

contains

  subroutine assert_status(label, result, expected)
    character(*), intent(in) :: label
    type(validation_result), intent(in) :: result
    integer, intent(in) :: expected

    call assert_true(label, result%status == expected)
  end subroutine assert_status
end program test_bh_toroid_updates
