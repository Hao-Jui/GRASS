module starting_model_mod
  use analysis_mod, only: solution_properties
  use iso_fortran_env, only: error_unit
  use precision_mod, only: wp
  use eos_mod, only: p_at_e, h_at_p, n0_at_h, e_at_h
  use regrid_mod, only: regrid_read
  use sphere_mod, only: sphere
  use para_mod, only: run_mode, MODE_REGRID, C, KSCALE, KAPPA, &
                      SDIV, MDIV, e_center, p_center, h_center, &
                      solver_type, r_ratio, shooting, SHOOT_2D, output, &
                      active_theory, THEORY_GR, mphi_goal, &
                      mphi_burn_threshold, B_coup, B_goal, mphi_r, l_uni
  use rotation_solver_mod,  only: rotation_solver
  use miscellaneous_mod, only: print_converged_block
  use scalar_burning_mod, only: perform_scalar_burn
  implicit none
contains
  subroutine initialize_starting_model()
    external :: restart_read, refine_read
    real(wp) :: target_mphi
    integer :: regrid_status
    character(len=256) :: regrid_error

    select case (run_mode)
    case (MODE_REGRID)
      call regrid_read(SDIV, MDIV, 2, regrid_status, regrid_error)
      if (regrid_status /= 0) then
        write(error_unit, '(A)') trim(regrid_error)
        error stop "initialize_starting_model: regrid_read failed"
      end if
      e_center = .74e15_wp
      e_center = e_center * C * C * KSCALE
      p_center = p_at_e(e_center)
      h_center = h_at_p(p_center)
    case default
      r_ratio  = merge(0.9e0_wp, 1.e0_wp, trim(adjustl(solver_type)) == "uryu")
      e_center = .8e15_wp
      e_center = e_center * C * C * KSCALE
      p_center = p_at_e(e_center)
      h_center = h_at_p(p_center)
      call sphere

      if (active_theory /= THEORY_GR .and. mphi_goal > mphi_burn_threshold) then
        target_mphi = mphi_goal
        call perform_scalar_burn(target_mphi)
      end if
    end select

    if (shooting == SHOOT_2D .and. abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) r_ratio = min(r_ratio, 0.95e0_wp)

    if (active_theory /= THEORY_GR) then
      B_coup  = B_goal
      mphi_r  = (mphi_goal / l_uni)**2 * KAPPA / 1.e10_wp
    end if

    call single_model()
  contains
    subroutine single_model()
      use para_mod, only: run_task, OneModel
      use constrain_mod, only: hamiltonian
      real(wp) :: ee, rho0, hamL2, t0, t1
      r_ratio = 0.6_wp

      output = .true.; call cpu_time(t0)
          call rotation_solver
          call solution_properties
          rho0 = n0_at_h(h_center)
          ee   = e_at_h(h_center)
          call print_converged_block(rho0, ee)
      output = .false.; call cpu_time(t1)

      write(*,*) " "; write(*,"(A, f10.4)") "Elapsed time [s]: ", t1-t0
      call hamiltonian(hamL2)
      write(*,"(A18,es27.16)") "Ham L2:", hamL2; write(*,*) " "
      if (run_task == OneModel) stop "One model solved!"
    end subroutine single_model
  end subroutine initialize_starting_model
end module starting_model_mod
