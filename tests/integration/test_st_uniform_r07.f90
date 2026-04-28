program test_st_uniform_r07
  use eos_mod, only: loadEos
  use grid_mod, only: make_grid, GridTrig
  use para_mod, only: wp, initialize_theory, &
                      active_theory, THEORY_ST, &
                      solver_type, &
                      SDIV, MDIV, &
                      shooting, SHOOT_FIX1_HC, &
                      run_mode, MODE_DEFAULT, &
                      run_task, OneModel, &
                      angular_collocation, COLLOCATION_LEG, &
                      eos_file, &
                      B_goal, mphi_goal, &
                      s_gp, s_pwr
  use starting_model_mod, only: initialize_starting_model
  implicit none

  ! -- Config overrides -------------------------------------------------------
  active_theory      = THEORY_ST
  solver_type        = "uniform"
  eos_file           = "MPA1"
  SDIV               = 401
  MDIV               = 41
  shooting           = SHOOT_FIX1_HC
  run_mode           = MODE_DEFAULT
  run_task           = OneModel
  angular_collocation = COLLOCATION_LEG
  B_goal             = 500.0_wp
  mphi_goal          = 0.2_wp

  ! -- Initialization sequence (mirrors src/main.f90) -------------------------
  call initialize_theory()
  call loadEos
  call make_grid
  call GridTrig
  call initialize_starting_model()

  ! initialize_starting_model calls stop "One model solved!" on success.
  error stop "test_st_uniform_r07: solver returned without stopping"
end program test_st_uniform_r07
