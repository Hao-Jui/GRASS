program test_gr_uryu
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use analysis_mod, only: solution_properties
  use eos_mod, only: loadEos, p_at_e, h_at_p, n0_at_h, e_at_h
  use grid_mod, only: make_grid, GridTrig
  use miscellaneous_mod, only: print_converged_block
  use para_mod, only: wp, initialize_theory, &
                      active_theory, THEORY_GR, &
                      solver_type, lambda1, lambda2, uryu_p, uryu_q, &
                      SDIV, MDIV, &
                      shooting, SHOOT_2D, &
                      run_mode, MODE_DEFAULT, &
                      run_task, OneModel, &
                      angular_collocation, COLLOCATION_LEG, &
                      eos_file, &
                      C, KSCALE, e_center, p_center, h_center, r_ratio, output
  use rotation_solver_mod, only: rotation_solver
  use sphere_mod, only: sphere
  implicit none
  real(wp) :: ee, rho0

  ! -- Config overrides -------------------------------------------------------
  active_theory      = THEORY_GR
  solver_type        = "uryu"
  lambda1            = 1.5_wp
  lambda2            = 0.3_wp
  uryu_p             = 1
  uryu_q             = 3
  eos_file           = "MPA1"
  SDIV               = 401
  MDIV               = 41
  shooting           = SHOOT_2D
  run_mode           = MODE_DEFAULT
  run_task           = OneModel
  angular_collocation = COLLOCATION_LEG

  ! -- Initialization sequence (mirrors src/main.f90) -------------------------
  call initialize_theory()
  call loadEos
  call make_grid
  call GridTrig

  ! Set the model explicitly so para_panel and starting-model defaults cannot
  ! silently change the regression case.
  e_center = 0.8e15_wp * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
  r_ratio  = 0.7_wp

  call sphere
  output = .true.
  call rotation_solver
  call solution_properties
  rho0 = n0_at_h(h_center)
  ee   = e_at_h(h_center)
  if (.not. ieee_is_finite(rho0) .or. .not. ieee_is_finite(ee)) &
    error stop "test_gr_uryu: non-finite converged state"
  call print_converged_block(rho0, ee)
  output = .false.
  write(*, '(A)') "test_gr_uryu: full-model regression output written"
end program test_gr_uryu
