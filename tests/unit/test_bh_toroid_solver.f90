program test_bh_toroid_solver
  use precision_mod, only: wp
  use bh_toroid_solver_mod
  use bh_toroid_updates_mod, only: bh_toroid_equatorial_point
  use bh_toroid_validation_mod, only: validation_result, VALID_OK, VALID_BAD_GRID_SIZE, &
      VALID_BAD_RADIAL_ORDER, VALID_BAD_SCALE
  use test_utils
  implicit none

  real(wp), parameter :: tol = 1.e-11_wp
  type(bh_toroid_solver_config) :: config, bad_config
  type(bh_toroid_solver_state) :: state, state_again
  type(bh_toroid_solver_diagnostics) :: diagnostics, diagnostics_again
  type(bh_toroid_integrals) :: integrals
  type(bh_horizon_quantities) :: horizon
  type(bh_toroid_equatorial_point) :: h_point, s_point, t_point
  type(validation_result) :: res
  real(wp) :: expected_mass, expected_j

  h_point = bh_toroid_equatorial_point(rhat=0.2_wp, nu_hat=0.0_wp, gamma_hat=0.0_wp, omega_hat=0.0_wp)
  s_point = bh_toroid_equatorial_point(rhat=0.5_wp, nu_hat=0.0_wp, gamma_hat=0.0_wp, omega_hat=0.0_wp)
  t_point = bh_toroid_equatorial_point(rhat=1.0_wp, nu_hat=0.0_wp, gamma_hat=0.0_wp, omega_hat=0.0_wp)
  config = bh_toroid_solver_config(n_r=4, n_theta=2, max_iterations=4, h0_hat=0.2_wp, &
      rin_hat=0.5_wp, r_out=2.0_wp, rotation_A=3.0_wp, omega_h=0.0_wp, poly_k=0.5_wp, &
      poly_n=2.0_wp, relaxation_factor=0.25_wp, tolerance=1.e-8_wp)

  call assert_status("valid solver config", validate_bh_toroid_solver_config(config), VALID_OK)
  bad_config = config
  bad_config%n_r = 1
  call assert_status("bad solver grid", validate_bh_toroid_solver_config(bad_config), VALID_BAD_GRID_SIZE)
  bad_config = config
  bad_config%rin_hat = config%h0_hat
  call assert_status("bad solver radial order", validate_bh_toroid_solver_config(bad_config), VALID_BAD_RADIAL_ORDER)
  bad_config = config
  bad_config%relaxation_factor = 0.0_wp
  call assert_status("bad solver relaxation", validate_bh_toroid_solver_config(bad_config), VALID_BAD_SCALE)
  bad_config = config
  bad_config%tolerance = 0.0_wp
  call assert_status("bad solver tolerance", validate_bh_toroid_solver_config(bad_config), VALID_BAD_SCALE)
  bad_config = config
  bad_config%max_iterations = 0
  call assert_status("bad solver max iterations", validate_bh_toroid_solver_config(bad_config), &
      VALID_BAD_SCALE)
  h_point%nu_hat = huge(1.0_wp)
  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state)
  call assert_status("bad equatorial point rejects initialization", res, VALID_BAD_SCALE)
  h_point%nu_hat = 0.0_wp

  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state)
  call assert_status("initialize solver", res, VALID_OK)
  call assert_true("r grid size", size(state%rhat) == config%n_r)
  call assert_true("theta grid size", size(state%sin_theta) == config%n_theta)
  call assert_true("field shape", all(shape(state%gamma_hat) == [config%n_r, config%n_theta]))
  call assert_true("initial fields finite", all(abs(state%gamma_hat) < huge(1.0_wp)))
  call assert_true("initial iteration zero", state%iteration == 0)

  state%energy_density = 0.0_wp
  state%omega = 0.0_wp
  state%energy_density(3,1) = 2.0_wp
  state%omega(3,1) = 0.5_wp
  res = bh_toroid_solver_step(config, state, diagnostics)
  call assert_status("kernel-backed one step", res, VALID_OK)
  call assert_near("horizon omega target via step", 0.0_wp, state%omega_hat(1,1), tol)
  call assert_true("kernel-backed nu update", state%nu_hat(3,1) < 0.0_wp)
  call assert_true("kernel-backed gamma update", state%gamma_hat(3,1) > 0.0_wp)
  call assert_true("relaxation produced bounded delta", diagnostics%metric_delta > 0.0_wp .and. &
      diagnostics%metric_delta < huge(1.0_wp))

  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state)
  call assert_status("reinitialize solver", res, VALID_OK)
  state%energy_density(3,1) = 1.0_wp
  res = bh_toroid_solver_step(config, state, diagnostics)
  call assert_status("one solver step", res, VALID_OK)
  call assert_true("one step finite residual", diagnostics%residual < huge(1.0_wp))
  call assert_true("one step iteration", diagnostics%iteration == 1)

  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state_again)
  call assert_status("reinitialize deterministic solver", res, VALID_OK)
  state_again%energy_density(3,1) = 1.0_wp
  res = bh_toroid_solver_step(config, state_again, diagnostics_again)
  call assert_status("deterministic solver step", res, VALID_OK)
  call assert_near("deterministic residual", diagnostics%residual, diagnostics_again%residual, tol)
  call assert_near("deterministic mass", diagnostics%integrals%mass, diagnostics_again%integrals%mass, tol)

  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state)
  call assert_status("initialize stationary solver", res, VALID_OK)
  res = solve_bh_toroid(config, state, diagnostics)
  call assert_status("stationary solve converges", res, VALID_OK)
  call assert_true("stationary convergence flag", diagnostics%converged)
  call assert_true("lowres fixture finite omega", all(abs(state%omega) < huge(1.0_wp)))
  call assert_true("lowres fixture finite enthalpy", all(abs(state%enthalpy_term) < huge(1.0_wp)))
  call assert_near("lowres inner edge density", 0.0_wp, state%energy_density(2,1), tol)
  call assert_near("lowres outer edge density", 0.0_wp, state%energy_density(4,1), tol)
  call assert_near("lowres outer nu asymptotic flatness", 0.0_wp, state%nu_hat(config%n_r,1), tol)
  call assert_near("lowres outer gamma asymptotic flatness", 0.0_wp, &
      state%gamma_hat(config%n_r,1), tol)
  call assert_near("lowres outer omega asymptotic flatness", 0.0_wp, &
      state%omega_hat(config%n_r,1), tol)
  call assert_near("lowres horizon angular velocity", 0.0_wp, &
      diagnostics%horizon%omega_residual, tol)

  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state)
  call assert_status("initialize nonconvergence fixture", res, VALID_OK)
  state%energy_density(3,1) = 1.0_wp
  bad_config = config
  bad_config%max_iterations = 1
  bad_config%tolerance = tiny(1.0_wp)
  res = solve_bh_toroid(bad_config, state, diagnostics)
  call assert_true("bounded solve reports nonconvergence", res%status /= VALID_OK)

  res = initialize_bh_toroid_solver(config, h_point, s_point, t_point, state)
  call assert_status("initialize integral fixture", res, VALID_OK)
  state%rhat = [0.2_wp, 0.5_wp, 0.75_wp, 1.0_wp]
  state%radial_weights = [0.1_wp, 0.2_wp, 0.3_wp, 0.4_wp]
  state%sin_theta = [1.0_wp, 1.0_wp]
  state%energy_density = 0.0_wp
  state%omega = 0.0_wp
  state%energy_density(3,1) = 2.0_wp
  state%omega(3,1) = 0.5_wp
  res = compute_bh_toroid_integrals(config, state, integrals)
  call assert_status("compute integrals", res, VALID_OK)
  expected_mass = 2.0_wp * config%r_out**3 * 0.3_wp * 0.5_wp * 0.75_wp**2
  expected_j = expected_mass * 0.5_wp * (config%r_out * 0.75_wp)**2
  call assert_near("mass integral fixture", expected_mass, integrals%mass, tol)
  call assert_near("angular momentum fixture", expected_j, integrals%angular_momentum, tol)

  res = compute_bh_horizon_quantities(config, state, horizon)
  call assert_status("compute horizon quantities", res, VALID_OK)
  call assert_near("horizon radius", config%r_out * config%h0_hat, horizon%radius, tol)
  call assert_near("omega equals omega_h on horizon", 0.0_wp, horizon%omega_residual, tol)

  bad_config = config
  bad_config%omega_h = 0.2_wp
  res = initialize_bh_toroid_solver(bad_config, h_point, s_point, t_point, state)
  call assert_status("initialize nonzero horizon fixture", res, VALID_OK)
  state%omega_hat(1,:) = bad_config%omega_h / bad_config%r_out**2
  res = compute_bh_horizon_quantities(bad_config, state, horizon)
  call assert_status("compute nonzero horizon quantities", res, VALID_OK)
  call assert_near("nonzero omega equals omega_h", 0.0_wp, horizon%omega_residual, tol)

  state%energy_density = 0.0_wp
  state%omega = 0.0_wp
  res = bh_toroid_solver_step(config, state, diagnostics)
  call assert_status("finite-field invariant step", res, VALID_OK)
  call assert_true("finite gamma fields", all(abs(state%gamma_hat) < huge(1.0_wp)))
  call assert_true("finite hydro fields", all(abs(state%energy_density) < huge(1.0_wp)))
  call assert_near("inner boundary density", 0.0_wp, state%energy_density(2,1), tol)
  call assert_near("outer boundary density", 0.0_wp, state%energy_density(4,1), tol)

  call test_summary("test_bh_toroid_solver")

contains

  subroutine assert_status(label, result, expected)
    character(*), intent(in) :: label
    type(validation_result), intent(in) :: result
    integer, intent(in) :: expected

    call assert_true(label, result%status == expected)
  end subroutine assert_status
end program test_bh_toroid_solver
