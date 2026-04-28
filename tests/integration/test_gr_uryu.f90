program test_gr_uryu
  use eos_mod, only: loadEos
  use grid_mod, only: make_grid, GridTrig
  use para_mod, only: wp, initialize_theory, &
                      active_theory, THEORY_GR, &
                      solver_type, &
                      SDIV, MDIV, &
                      shooting, SHOOT_2D, &
                      run_mode, MODE_DEFAULT, &
                      run_task, OneModel, &
                      angular_collocation, COLLOCATION_LEG, &
                      eos_file, &
                      s_gp, s_pwr
  use starting_model_mod, only: initialize_starting_model
  implicit none

  ! -- Config overrides -------------------------------------------------------
  active_theory      = THEORY_GR
  solver_type        = "uryu"
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
  call initialize_starting_model()

  ! initialize_starting_model calls stop "One model solved!" on success.
  ! Reaching here means an unexpected return — treat as failure.
  error stop "test_gr_uryu: solver returned without stopping"
end program test_gr_uryu
