module solver_mod
  use precision_mod, only: wp
  use validation_mod, only: validation_result, validation_ok, validation_error, &
      is_finite, VALID_OK, VALID_BAD_GRID_SIZE, VALID_BAD_RADIAL_ORDER, VALID_BAD_SCALE
  use radial_map_mod, only: build_rhat_trapezoid_grid, classify_rhat, &
      ZONE_HORIZON, ZONE_TORUS
  use green_mod, only: ne_f1_kernel, ne_f2_kernel
  use updates_mod, only: bh_toroid_equatorial_point, bh_toroid_update_constants, &
      solve_hydro_rotation_constants, update_hydro_rotation_fields, finite_equatorial_point
  implicit none
  private

  type, public :: bh_toroid_solver_config
    integer :: n_r = 0
    integer :: n_theta = 0
    integer :: max_iterations = 100
    real(wp) :: h0_hat = 0.0_wp
    real(wp) :: rin_hat = 0.0_wp
    real(wp) :: r_out = 1.0_wp
    real(wp) :: rotation_A = 1.0_wp
    real(wp) :: omega_h = 0.0_wp
    real(wp) :: poly_k = 1.0_wp
    real(wp) :: poly_n = 3.0_wp
    real(wp) :: relaxation_factor = 0.5_wp
    real(wp) :: tolerance = 1.e-8_wp
  end type bh_toroid_solver_config

  type, public :: bh_toroid_integrals
    real(wp) :: mass = 0.0_wp
    real(wp) :: angular_momentum = 0.0_wp
  end type bh_toroid_integrals

  type, public :: bh_horizon_quantities
    real(wp) :: radius = 0.0_wp
    real(wp) :: omega_h = 0.0_wp
    real(wp) :: omega_model = 0.0_wp
    real(wp) :: omega_residual = 0.0_wp
    real(wp) :: mass = 0.0_wp
    real(wp) :: angular_momentum = 0.0_wp
    real(wp) :: area_proxy = 0.0_wp
  end type bh_horizon_quantities

  type, public :: bh_toroid_solver_diagnostics
    integer :: iteration = 0
    logical :: converged = .false.
    real(wp) :: residual = huge(1.0_wp)
    real(wp) :: metric_delta = huge(1.0_wp)
    real(wp) :: hydro_delta = huge(1.0_wp)
    type(bh_toroid_integrals) :: integrals
    type(bh_horizon_quantities) :: horizon
  end type bh_toroid_solver_diagnostics

  type, public :: bh_toroid_solver_state
    integer :: iteration = 0
    logical :: converged = .false.
    type(bh_toroid_update_constants) :: constants
    real(wp), allocatable :: rhat(:)
    real(wp), allocatable :: radial_weights(:)
    real(wp), allocatable :: sin_theta(:)
    real(wp), allocatable :: gamma_hat(:,:)
    real(wp), allocatable :: nu_hat(:,:)
    real(wp), allocatable :: omega_hat(:,:)
    real(wp), allocatable :: target_gamma_hat(:,:)
    real(wp), allocatable :: target_nu_hat(:,:)
    real(wp), allocatable :: target_omega_hat(:,:)
    real(wp), allocatable :: omega(:,:)
    real(wp), allocatable :: velocity(:,:)
    real(wp), allocatable :: enthalpy_term(:,:)
    real(wp), allocatable :: energy_density(:,:)
  end type bh_toroid_solver_state

  public :: validate_bh_toroid_solver_config, initialize_bh_toroid_solver
  public :: bh_toroid_solver_step, solve_bh_toroid
  public :: compute_bh_toroid_integrals, compute_bh_horizon_quantities

contains

  pure function validate_bh_toroid_solver_config(config) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(validation_result) :: res

    if (config%n_r < 2 .or. config%n_theta < 1) then
      res = validation_error(VALID_BAD_GRID_SIZE, "solver grid must have n_r >= 2 and n_theta >= 1")
    else if (.not. all(is_finite([config%h0_hat, config%rin_hat, config%r_out, config%rotation_A, &
        config%omega_h, config%poly_k, config%poly_n, config%relaxation_factor, config%tolerance]))) then
      res = validation_error(VALID_BAD_SCALE, "solver config values must be finite")
    else if (config%h0_hat <= 0.0_wp .or. config%h0_hat >= 1.0_wp .or. &
        config%rin_hat <= config%h0_hat .or. config%rin_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_RADIAL_ORDER, "solver domain must satisfy 0 < h0_hat < rin_hat < 1")
    else if (config%r_out <= 0.0_wp .or. config%rotation_A <= 0.0_wp .or. &
        config%poly_k <= 0.0_wp .or. config%poly_n <= 0.0_wp) then
      res = validation_error(VALID_BAD_SCALE, "solver scales, rotation_A, and polytrope values must be positive")
    else if (config%relaxation_factor <= 0.0_wp .or. config%relaxation_factor > 1.0_wp) then
      res = validation_error(VALID_BAD_SCALE, "relaxation_factor must be in (0, 1]")
    else if (config%tolerance <= 0.0_wp .or. config%max_iterations < 1) then
      res = validation_error(VALID_BAD_SCALE, "tolerance must be positive and max_iterations >= 1")
    else
      res = validation_ok()
    end if
  end function validate_bh_toroid_solver_config

  function initialize_bh_toroid_solver(config, h_point, s_point, t_point, state) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_equatorial_point), intent(in) :: h_point, s_point, t_point
    type(bh_toroid_solver_state), intent(inout) :: state
    type(validation_result) :: res
    integer :: j
    real(wp), parameter :: pi = acos(-1.0_wp)

    res = validate_bh_toroid_solver_config(config)
    if (res%status /= VALID_OK) return
    if (.not. finite_equatorial_point(h_point) .or. .not. finite_equatorial_point(s_point) .or. &
        .not. finite_equatorial_point(t_point)) then
      res = validation_error(VALID_BAD_SCALE, "equatorial points must be finite")
      return
    end if

    call release_state_arrays(state)
    allocate(state%rhat(config%n_r), state%radial_weights(config%n_r), state%sin_theta(config%n_theta))
    allocate(state%gamma_hat(config%n_r, config%n_theta), state%nu_hat(config%n_r, config%n_theta))
    allocate(state%omega_hat(config%n_r, config%n_theta))
    allocate(state%target_gamma_hat(config%n_r, config%n_theta))
    allocate(state%target_nu_hat(config%n_r, config%n_theta))
    allocate(state%target_omega_hat(config%n_r, config%n_theta))
    allocate(state%omega(config%n_r, config%n_theta), state%velocity(config%n_r, config%n_theta))
    allocate(state%enthalpy_term(config%n_r, config%n_theta))
    allocate(state%energy_density(config%n_r, config%n_theta))

    res = build_rhat_trapezoid_grid(config%n_r, config%h0_hat, state%rhat, state%radial_weights)
    if (res%status /= VALID_OK) return

    do j = 1, config%n_theta
      state%sin_theta(j) = sin((real(j, wp) - 0.5_wp) * 0.5_wp * pi / real(config%n_theta, wp))
    end do

    state%constants = solve_hydro_rotation_constants(h_point, s_point, t_point, config%h0_hat, &
        config%rin_hat, config%r_out, config%rotation_A, config%omega_h, config%poly_k, config%poly_n)
    if (.not. valid_update_constants(state%constants)) then
      res = validation_error(VALID_BAD_SCALE, "hydro-rotation constants are invalid")
      return
    end if

    state%iteration = 0
    state%converged = .false.
    state%gamma_hat = 0.0_wp
    state%nu_hat = 0.0_wp
    state%omega_hat = 0.0_wp
    state%target_gamma_hat = 0.0_wp
    state%target_nu_hat = 0.0_wp
    state%target_omega_hat = 0.0_wp
    state%omega = 0.0_wp
    state%velocity = 0.0_wp
    state%enthalpy_term = 0.0_wp
    state%energy_density = 0.0_wp
    res = validation_ok()
  end function initialize_bh_toroid_solver

  function assemble_bh_toroid_targets(config, state) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(inout) :: state
    type(validation_result) :: res
    integer :: i, j, ip, jp, zone
    real(wp) :: source, horizon_target, angular_weight, volume_weight, r, rp, h0

    res = validate_solver_state(config, state)
    if (res%status /= VALID_OK) return

    state%target_gamma_hat = 0.0_wp
    state%target_nu_hat = 0.0_wp
    state%target_omega_hat = 0.0_wp
    horizon_target = config%omega_h / max(config%r_out**2, tiny(1.0_wp))
    h0 = config%r_out * config%h0_hat
    angular_weight = 1.0_wp / real(config%n_theta, wp)

    do j = 1, config%n_theta
      do i = 1, config%n_r
        zone = classify_rhat(state%rhat(i), config%h0_hat, config%rin_hat)
        if (zone == ZONE_HORIZON) then
          state%target_omega_hat(i,j) = horizon_target
          cycle
        end if
        r = config%r_out * state%rhat(i)
        do jp = 1, config%n_theta
          do ip = 1, config%n_r
            if (classify_rhat(state%rhat(ip), config%h0_hat, config%rin_hat) /= ZONE_TORUS) cycle
            source = max(0.0_wp, state%energy_density(ip,jp))
            if (source <= 0.0_wp) cycle
            rp = config%r_out * state%rhat(ip)
            volume_weight = state%radial_weights(ip) * angular_weight * state%sin_theta(jp)
            state%target_nu_hat(i,j) = state%target_nu_hat(i,j) &
                - ne_f2_kernel(0, r, rp, h0) * source * volume_weight
            state%target_gamma_hat(i,j) = state%target_gamma_hat(i,j) &
                + ne_f1_kernel(1, r, rp, h0) * source * volume_weight * state%sin_theta(jp)
            state%target_omega_hat(i,j) = state%target_omega_hat(i,j) &
                + ne_f2_kernel(0, r, rp, h0) * source * state%omega(ip,jp) * volume_weight
          end do
        end do
      end do
    end do
    res = validation_ok()
  end function assemble_bh_toroid_targets

  function relax_bh_toroid_targets(config, state, max_delta) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(inout) :: state
    real(wp), intent(out) :: max_delta
    type(validation_result) :: res

    res = validate_solver_state(config, state)
    if (res%status /= VALID_OK) then
      max_delta = huge(1.0_wp)
      return
    end if

    max_delta = 0.0_wp
    call relax_in_place(state%gamma_hat, state%target_gamma_hat, config%relaxation_factor, max_delta)
    call relax_in_place(state%nu_hat, state%target_nu_hat, config%relaxation_factor, max_delta)
    call relax_in_place(state%omega_hat, state%target_omega_hat, config%relaxation_factor, max_delta)
    res = validation_ok()
  end function relax_bh_toroid_targets

  function bh_toroid_solver_step(config, state, diagnostics) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(inout) :: state
    type(bh_toroid_solver_diagnostics), intent(out) :: diagnostics
    type(validation_result) :: res
    real(wp), allocatable :: old_energy(:,:)
    real(wp) :: metric_delta, hydro_delta

    res = validate_solver_state(config, state)
    if (res%status /= VALID_OK) return

    allocate(old_energy(config%n_r, config%n_theta), source=state%energy_density)
    res = assemble_bh_toroid_targets(config, state)
    if (res%status /= VALID_OK) return

    res = relax_bh_toroid_targets(config, state, metric_delta)
    if (res%status /= VALID_OK) return

    res = update_hydro_rotation_fields(state%constants, config%rin_hat, state%rhat, state%sin_theta, &
        state%gamma_hat, state%nu_hat, state%omega_hat, state%omega, state%velocity, &
        state%enthalpy_term, state%energy_density)
    if (res%status /= VALID_OK) return

    hydro_delta = maxval(abs(state%energy_density - old_energy))
    state%iteration = state%iteration + 1
    state%converged = max(metric_delta, hydro_delta) <= config%tolerance

    diagnostics%iteration = state%iteration
    diagnostics%converged = state%converged
    diagnostics%metric_delta = metric_delta
    diagnostics%hydro_delta = hydro_delta
    diagnostics%residual = max(metric_delta, hydro_delta)
    res = compute_bh_toroid_integrals(config, state, diagnostics%integrals)
    if (res%status /= VALID_OK) return
    res = compute_bh_horizon_quantities(config, state, diagnostics%horizon)
  end function bh_toroid_solver_step

  function solve_bh_toroid(config, state, diagnostics) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(inout) :: state
    type(bh_toroid_solver_diagnostics), intent(out) :: diagnostics
    type(validation_result) :: res
    integer :: iter

    res = validate_solver_state(config, state)
    if (res%status /= VALID_OK) return

    do iter = 1, config%max_iterations
      res = bh_toroid_solver_step(config, state, diagnostics)
      if (res%status /= VALID_OK) return
      if (diagnostics%converged) return
    end do
    res = validation_error(VALID_BAD_SCALE, "BH toroid solver did not converge within max_iterations")
  end function solve_bh_toroid

  function compute_bh_toroid_integrals(config, state, integrals) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(in) :: state
    type(bh_toroid_integrals), intent(out) :: integrals
    type(validation_result) :: res
    integer :: i, j
    real(wp) :: angular_weight, volume_weight, lever_arm_squared

    res = validate_solver_state(config, state)
    if (res%status /= VALID_OK) return

    integrals%mass = 0.0_wp
    integrals%angular_momentum = 0.0_wp
    angular_weight = 1.0_wp / real(config%n_theta, wp)
    do j = 1, config%n_theta
      do i = 1, config%n_r
        if (classify_rhat(state%rhat(i), config%h0_hat, config%rin_hat) /= ZONE_TORUS) cycle
        ! Low-order quadrature shell for the matter terms in Nishida-Eriguchi
        ! eqs. (4.1)-(4.2). Metric factors are explicit state inputs.
        volume_weight = config%r_out**3 * state%radial_weights(i) * angular_weight &
            * state%rhat(i)**2 * state%sin_theta(j)
        lever_arm_squared = (config%r_out * state%rhat(i) * state%sin_theta(j))**2
        integrals%mass = integrals%mass + state%energy_density(i,j) * volume_weight
        integrals%angular_momentum = integrals%angular_momentum + state%energy_density(i,j) &
            * state%omega(i,j) * lever_arm_squared * volume_weight
      end do
    end do
    res = validation_ok()
  end function compute_bh_toroid_integrals

  function compute_bh_horizon_quantities(config, state, horizon) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(in) :: state
    type(bh_horizon_quantities), intent(out) :: horizon
    type(validation_result) :: res
    integer :: j_eq
    real(wp), parameter :: pi = acos(-1.0_wp)

    res = validate_solver_state(config, state)
    if (res%status /= VALID_OK) return

    j_eq = maxloc(state%sin_theta, dim=1)
    horizon%radius = config%r_out * config%h0_hat
    horizon%omega_h = config%omega_h
    horizon%omega_model = config%r_out**2 * state%omega_hat(1,j_eq)
    horizon%omega_residual = horizon%omega_model - config%omega_h
    horizon%area_proxy = 4.0_wp * pi * horizon%radius**2
    res = validation_ok()
  end function compute_bh_horizon_quantities

  function validate_solver_state(config, state) result(res)
    type(bh_toroid_solver_config), intent(in) :: config
    type(bh_toroid_solver_state), intent(in) :: state
    type(validation_result) :: res

    res = validate_bh_toroid_solver_config(config)
    if (res%status /= VALID_OK) return

    if (.not. state_arrays_allocated(state)) then
      res = validation_error(VALID_BAD_GRID_SIZE, "solver state arrays are not allocated")
    else if (.not. state_shapes_match(state, config%n_r, config%n_theta)) then
      res = validation_error(VALID_BAD_GRID_SIZE, "solver state arrays have incompatible shapes")
    else if (.not. valid_update_constants(state%constants)) then
      res = validation_error(VALID_BAD_SCALE, "solver constants must be finite and valid")
    else if (.not. state_arrays_finite(state)) then
      res = validation_error(VALID_BAD_SCALE, "solver state arrays must be finite")
    else
      res = validation_ok()
    end if
  end function validate_solver_state

  subroutine release_state_arrays(state)
    type(bh_toroid_solver_state), intent(inout) :: state

    if (allocated(state%rhat)) deallocate(state%rhat, state%radial_weights, state%sin_theta, &
        state%gamma_hat, state%nu_hat, state%omega_hat, state%target_gamma_hat, state%target_nu_hat, &
        state%target_omega_hat, state%omega, state%velocity, state%enthalpy_term, state%energy_density)
  end subroutine release_state_arrays

  subroutine relax_in_place(field, target_field, factor, max_delta)
    real(wp), intent(inout) :: field(:,:)
    real(wp), intent(in) :: target_field(:,:)
    real(wp), intent(in) :: factor
    real(wp), intent(inout) :: max_delta
    real(wp) :: delta
    integer :: i, j

    do j = 1, size(field, 2)
      do i = 1, size(field, 1)
        delta = factor * (target_field(i,j) - field(i,j))
        field(i,j) = field(i,j) + delta
        max_delta = max(max_delta, abs(delta))
      end do
    end do
  end subroutine relax_in_place

  pure function state_arrays_allocated(state) result(ok)
    type(bh_toroid_solver_state), intent(in) :: state
    logical :: ok

    ok = allocated(state%rhat) .and. allocated(state%radial_weights) .and. &
        allocated(state%sin_theta) .and. allocated(state%gamma_hat) .and. &
        allocated(state%nu_hat) .and. allocated(state%omega_hat) .and. &
        allocated(state%target_gamma_hat) .and. allocated(state%target_nu_hat) .and. &
        allocated(state%target_omega_hat) .and. allocated(state%omega) .and. &
        allocated(state%velocity) .and. allocated(state%enthalpy_term) .and. &
        allocated(state%energy_density)
  end function state_arrays_allocated

  pure function state_shapes_match(state, n_r, n_theta) result(ok)
    type(bh_toroid_solver_state), intent(in) :: state
    integer, intent(in) :: n_r, n_theta
    logical :: ok
    integer :: expected(2)

    expected = [n_r, n_theta]
    ok = size(state%rhat) == n_r .and. size(state%radial_weights) == n_r .and. &
        size(state%sin_theta) == n_theta .and. &
        all(shape(state%gamma_hat) == expected) .and. all(shape(state%nu_hat) == expected) .and. &
        all(shape(state%omega_hat) == expected) .and. all(shape(state%target_gamma_hat) == expected) .and. &
        all(shape(state%target_nu_hat) == expected) .and. all(shape(state%target_omega_hat) == expected) .and. &
        all(shape(state%omega) == expected) .and. all(shape(state%velocity) == expected) .and. &
        all(shape(state%enthalpy_term) == expected) .and. all(shape(state%energy_density) == expected)
  end function state_shapes_match

  pure function state_arrays_finite(state) result(ok)
    type(bh_toroid_solver_state), intent(in) :: state
    logical :: ok

    ok = all(is_finite(state%rhat)) .and. all(is_finite(state%radial_weights)) .and. &
        all(is_finite(state%sin_theta)) .and. all(is_finite(state%gamma_hat)) .and. &
        all(is_finite(state%nu_hat)) .and. all(is_finite(state%omega_hat)) .and. &
        all(is_finite(state%target_gamma_hat)) .and. all(is_finite(state%target_nu_hat)) .and. &
        all(is_finite(state%target_omega_hat)) .and. all(is_finite(state%omega)) .and. &
        all(is_finite(state%velocity)) .and. all(is_finite(state%enthalpy_term)) .and. &
        all(is_finite(state%energy_density))
  end function state_arrays_finite

  pure function valid_update_constants(constants) result(ok)
    type(bh_toroid_update_constants), intent(in) :: constants
    logical :: ok

    ok = all(is_finite([constants%r_out, constants%rotation_A, constants%poly_k, &
        constants%poly_n, constants%omega_c, constants%bernoulli_c, constants%surface_mismatch])) &
        .and. constants%r_out > 0.0_wp .and. constants%rotation_A > 0.0_wp &
        .and. constants%poly_k > 0.0_wp .and. constants%poly_n > 0.0_wp
  end function valid_update_constants

end module solver_mod
