module starting_model_mod
  use analysis_mod, only: solution_properties
  use precision_mod, only: wp
  use eos_mod, only: p_at_e, h_at_p, n0_at_h, e_at_h
  use regrid_mod, only: regrid_read
  use sphere_mod, only: sphere
  use para_mod, only: run_mode, MODE_REGRID, &
                      SDIV, MDIV, e_center, p_center, h_center, &
                      C, KSCALE, KAPPA, &
                      solver_type, r_ratio, use_shoot_1d, output, &
                      active_theory, THEORY_GR, &
                      mphi_goal, mphi_burn_threshold, &
                      B_coup, B_goal, mphi_r, l_uni
  use rotation_uniform,  only: rotation_solver
  use miscellaneous_mod, only: print_converged_block
  use scalar_burning_mod, only: perform_scalar_burn
  implicit none
contains
  subroutine initialize_starting_model()
    external :: restart_read, refine_read
    real(wp) :: target_mphi

    select case (run_mode)
    case (MODE_REGRID)
      call regrid_read(SDIV, MDIV, 2)
      e_center = e_center * C * C * KSCALE
      p_center = p_at_e(e_center)
      h_center = h_at_p(p_center)
    case default
      r_ratio  = merge(0.9e0_wp, 1.e0_wp, trim(adjustl(solver_type)) == "uryu")
      e_center = 8.e14_wp
      e_center = e_center * C * C * KSCALE
      p_center = p_at_e(e_center)
      h_center = h_at_p(p_center)
      call sphere

      if (active_theory /= THEORY_GR .and. mphi_goal > mphi_burn_threshold) then
        target_mphi = mphi_goal
        call perform_scalar_burn(target_mphi)
      end if
    end select

    if (.not. use_shoot_1d .and. abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) r_ratio = min(r_ratio, 0.9e0_wp)

    if (active_theory /= THEORY_GR) then
      B_coup  = B_goal
      mphi_r  = (mphi_goal / l_uni)**2 * KAPPA / 1.e10_wp
    end if

    call single_model()
  contains
    subroutine single_model()
      real(wp) :: ee, rho0
      !r_ratio = .7e0_wp
      output = .true.
      call rotation_solver
      call solution_properties

      rho0 = n0_at_h(h_center)
      ee   = e_at_h(h_center)
      call print_converged_block(rho0, ee)
      output = .false.
      print *, " "
      !stop "One model solved!"
    end subroutine single_model
  end subroutine initialize_starting_model
end module starting_model_mod
