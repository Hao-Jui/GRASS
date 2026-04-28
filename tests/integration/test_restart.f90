! Restart round-trip test:
!   1. Solve GR+uniform from scratch → writes Res/res.rst via output_helper
!   2. Read Res/res.rst via MODE_REGRID → re-solve → verify convergence
program test_restart
  use eos_mod, only: loadEos
  use grid_mod, only: make_grid, GridTrig
  use para_mod, only: wp, initialize_theory, &
                      active_theory, THEORY_GR, &
                      solver_type, &
                      SDIV, MDIV, &
                      shooting, SHOOT_FIX1_HC, &
                      run_mode, MODE_DEFAULT, MODE_REGRID, &
                      run_task, OneModel, shoot, &
                      angular_collocation, COLLOCATION_LEG, &
                      eos_file, s_gp, s_pwr, output
  use starting_model_mod, only: initialize_starting_model
  implicit none

  ! --- Phase 1: Solve from scratch (writes Res/res.rst) ---
  active_theory      = THEORY_GR
  solver_type        = "uniform"
  eos_file           = "MPA1"
  SDIV               = 201
  MDIV               = 21
  shooting           = SHOOT_FIX1_HC
  run_mode           = MODE_DEFAULT
  run_task           = shoot
  angular_collocation = COLLOCATION_LEG

  call initialize_theory()
  call loadEos
  call make_grid
  call GridTrig

  ! run_task=shoot → single_model returns instead of calling stop
  call initialize_starting_model()

  ! Res/res.rst was written by output_helper inside rotation_solver

  ! --- Phase 2: Restart from Res/res.rst ---
  run_mode = MODE_REGRID
  run_task = OneModel

  call make_grid
  call GridTrig
  call initialize_starting_model()

  ! If we reach here, solver returned without converging
  error stop "test_restart: solver returned without stopping"
end program test_restart
