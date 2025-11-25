subroutine initialize_starting_model(p_at_e, h_at_p)
  use para_mod
  use rotation_dispatch, only: call_rotation_solver
  implicit none
  real(8), external :: p_at_e, h_at_p
  external :: sphere, restart_read, refine_read, regrid_read
  real(8) :: target_mphi

  select case (run_mode)

  case (MODE_REGRID)
    call regrid_read(SDIV, MDIV, 2)
    e_center = e_center * C * C * KSCALE
    p_center = p_at_e(e_center)
    h_center = h_at_p(p_center)
  case default
    r_ratio  = 1.d0
    e_center = 8.6405005E+14
    if (.not. use_shoot_1d) r_ratio = min(r_ratio, 0.9d0)
    e_center = e_center * C * C * KSCALE
    p_center = p_at_e(e_center)
    h_center = h_at_p(p_center)
    call sphere

    if (active_theory /= THEORY_GR .and. mphi_goal > mphi_burn_threshold) then
      target_mphi = mphi_goal
      call perform_scalar_burn(target_mphi)
    end if
  end select
  !call single_model()
  if (active_theory /= THEORY_GR) then
    B_coup  = B_goal
    mphi_r  = (mphi_goal / l_uni)**2 * KAPPA / 1.d10
  end if
contains
  subroutine single_model()
    real(8) :: n0_at_h, e_at_h, ee, rho0
    r_ratio  = .97d0
    e_center = 2.d16
    e_center = e_center * C * C * KSCALE
    p_center = p_at_e(e_center)
    h_center = h_at_p(p_center)
    output = .true.
    call call_rotation_solver
    call mass_radius
    rho0 = n0_at_h(h_center)
    ee   = e_at_h(h_center)
    call print_converged_block(rho0, ee)
    print *, " "
    stop " Finish single model call"
  end subroutine single_model
end subroutine initialize_starting_model

subroutine shoot_v2
  use para_mod
  use shoot_solver_mod
  use shoot_solver_mod_1d
  use rotation_dispatch, only: call_rotation_solver
  use miscellaneous_mod, only: log_kepler_sequence
  use shoot_newton_helpers, only: evaluate_solution, build_jacobian
  use shoot_newton_helpers_1d, only: evaluate_solution_1d, build_jacobian_1d
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it, i_idx, iteration_cap
  real(8) :: er, rho0, ee
  real(8) :: F(2), x(2), delta_x(2), rhs(2), h_new, r_new
  real(8) :: F1d, x1d, delta_x1d, rhs1d
  type(newton_state)    :: solver_state
  type(newton_state_1d) :: solver_state_1d
  logical :: step_ok
  real(8) :: prev_er1d, prev_er2d
  external :: p_at_e, h_at_p, n0_at_h, e_at_h

  call initialize_starting_model(p_at_e, h_at_p)
  if (.not. use_shoot_1d) then
    call init_newton_state(solver_state, 2)
  end if

  write(unit=*, fmt=*) " "

  iteration_cap = 1

  do i_idx = 1, iteration_cap
    it = 1
    er = 1.d99
    if (use_shoot_1d) then
      call reset_newton_state_1d(solver_state_1d)
      prev_er1d = huge(1.d0)
    else
      call reset_newton_state(solver_state)
      prev_er2d = huge(1.d0)
    endif
    output = .false.

    do
      call cpu_time(start)
      n_of_relaxation_steps = 0
      if (use_shoot_1d) then
        call evaluate_solution_1d(h_center, r_ratio, F1d, rho0, ee)
        call to_solver_coord_1d(h_center, x1d)

        if (solver_state_1d%has_jacobian) then
          call broyden_update_1d(solver_state_1d, x1d, F1d)
        else
          call build_jacobian_1d(solver_state_1d, x1d, F1d, h_center, r_ratio, rho0, ee, reuse_base=.true.)
          call to_solver_coord_1d(h_center, x1d)
        endif

        call commit_state_1d(solver_state_1d, x1d, F1d)
        er = abs(F1d)
      else
        call evaluate_solution(h_center, r_ratio, F, rho0, ee, er)
        call to_solver_coords(h_center, r_ratio, x)

        if (solver_state%has_jacobian) then
          call broyden_update(solver_state, x, F)
        else
          call build_jacobian(solver_state, x, F, h_center, r_ratio, rho0, ee, er, reuse_base=.true.)
          write(*,*) "Jacobian built"
          call to_solver_coords(h_center, r_ratio, x)
        endif

        call commit_state(solver_state, x, F)
      end if

      !if (mod(it, 10) == 0) call print_iter_status(it, rho0, ee, er)
      call cpu_time(finish)
      write(*,"(A,ES15.6,5X,A,I5,5X,A,F8.3)") &
            "Error: ", er, &
            "Iterations: ", n_of_relaxation_steps, &
            "Elapsed Time [s]: ", finish-start
     
      if (er <= accuracy .or. output) exit

      if (use_shoot_1d) then
        rhs1d = -F1d
        step_ok = solve_linear_1d(solver_state_1d%J, rhs1d, delta_x1d)
        if (.not. step_ok) then
          call reset_newton_state_1d(solver_state_1d)
          cycle
        endif

        if (solver_state_1d%has_jacobian .and. abs(delta_x1d) < 1.d-6 .and. abs(F1d) > accuracy*1.d2) then
          call reset_newton_state_1d(solver_state_1d)
          cycle
        end if

        call clamp_step_1d(delta_x1d)
        x1d = x1d + delta_x1d
        call from_solver_coord_1d(x1d, h_new)
        h_center = h_new
        solver_state_1d%has_jacobian = .true.
      else
        rhs = -F
        step_ok = solve_linear(solver_state%J, rhs, delta_x)
        if (.not. step_ok) then
          call reset_newton_state(solver_state)
          cycle
        endif

        if (solver_state%has_jacobian .and. maxval(abs(delta_x(1:2))) < 1.d-4 &
            .and. maxval(abs(F)) > accuracy*1.d2) then
          call reset_newton_state(solver_state)
          write(*,*) "Jacobian refreshed"
          cycle
        end if

        call clamp_step(delta_x)
        x = x + delta_x
        call from_solver_coords(x, h_new, r_new)
        h_center = h_new
        r_ratio  = r_new
        solver_state%has_jacobian = .true.
      end if
      
      it = it + 1
      if (it == 1000) stop "Iteration may never converge."
    end do

    output = .true.
    call call_rotation_solver
    call solution_properties

    rho0 = n0_at_h(h_center)
    ee   = e_at_h(h_center)
    write(unit=*, fmt=*) " "
    call print_converged_block(rho0, ee)

    !call log_kepler_sequence()
  end do ! looping models

end subroutine shoot_v2

subroutine print_iter_status(it, rho0, ee, er)
  use para_mod
  implicit none
  integer, intent(in) :: it
  real(8), intent(in) :: rho0, ee, er

  write(*,'(1X,A)') repeat('=', 36)
  write(*,"(A10,I6)")           " iter :", it/10
  if (active_theory /= THEORY_GR) then
    write(*,"(A10,ES18.9)")     " Bcoup:", B_coup
    write(*,"(A10,ES18.9)")     " mphi :", sqrt(mphi_r*1.d10/KAPPA)*l_uni
    write(*,"(A10,ES18.9)")     " sphi :", sphi_c
  end if
  write(*,"(A10,ES18.9,A10)")    "rho_c :", rho0*MB,"g/cm^3"
  write(*,"(A10,ES18.9,A10,ES18.9)") "e_c   :", ee/(C*C*KSCALE),"g/cm^3", ee/(C*C*KSCALE)*rho_uni
  write(*,"(A10,ES18.9)")       "rp/re :", r_ratio
  write(*,"(A10,ES18.9,A4)")    "OMG_c :", omega_c/(2.d0*pi)*(C/sqrt(kappa)),"Hz"
  write(*,"(A10,ES18.9,A4)")    "Mass  :", mass/MSUN,"M_o"
  write(*,"(A10,ES18.9,A4)")    "M_0   :", mass_0/MSUN,"M_o"
  write(*,"(A10,ES18.9)")       "J     :", ang_mom
  write(*,"(A10,ES18.9)")       "chi   :", chi
  write(*,"(A10,ES18.9)")       "T/W   :", T_kin/abs(mass_p - mass + T_kin)
  write(*,"(A10,ES18.9,A4)")    "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
  write(*,"(A10,ES18.9,A4)")    "R_cir :", r_circ/1.d5,"km"
  write(*,"(A10,ES18.9,A4)")    "r_e/M :", (r_circ/1.d5)/(mass/MSUN*1.4769994423016508d0)
  if (.not. active_theory == THEORY_GR) then
    write(*,"(A10,2ES18.9)")    "Om_K/e:", Omega_K/(2.d0*pi), Omega_e/(2.d0*pi)*(C/sqrt(kappa))
  end if
  write(*,"(A10,3ES18.9)")      "er    :", er
  write(*,'(1X,A)') repeat('=', 36)
  write(*, *) " "
end subroutine print_iter_status

subroutine print_converged_block(rho0, ee)
#include "option_macro.h"
  use para_mod
  implicit none
  real(8), intent(in) :: rho0, ee
  integer :: i

  if (active_theory /= THEORY_GR) then
    open(221, file="./Cont/properties_st.dat")
  else
    open(221, file="./Cont/properties.dat")
  end if

  write(unit=*, fmt=*) " ===================================="
  write(unit=*, fmt=*) "              Converged              "
  do i = 1, 2
     write(6+215*(i-1),"(A18,ES18.9,A8,ES18.9)")  &
          "   Central rho =", rho0*MB,"g/cm^3", rho0*MB*rho_uni
     write(6+215*(i-1),"(A18,ES18.9,A10)")        &
          "Central energy =", ee/(C*C*KSCALE), "g/cm^3"
     write(6+215*(i-1),"(A18,ES18.9)")            "   Axial ratio =", r_ratio
     write(6+215*(i-1),"(A18,F18.9,A4)")          " Central Omega =", Omega_c/(2.d0*pi)*(C/sqrt(kappa)), "Hz"
     write(6+215*(i-1),"(A18,F18.9,A4)")          " Equator Omega =", Omega_e/(2.d0*pi)*(C/sqrt(kappa)), "Hz"
     if (active_theory == THEORY_GR) then
       write(6+215*(i-1),"(A18,F18.9,A4)")        "   Kepler Omega =", Omega_K/(2.d0*pi), "Hz"
     end if
     write(6+215*(i-1),"(A18,F18.9,A4)")          "      ADM Mass =", mass/MSUN, "M_o"
     write(6+215*(i-1),"(A18,F18.9,A16,ES18.9,A2)") &
          " Baryonic Mass =", mass_0/MSUN, "M_o ( binding =", (mass - mass_0)/MSUN, " )"
     write(6+215*(i-1),"(A18,F18.9,A10,F7.4,A2)") &
          " Ang. Momentum =", ang_mom, " ( chi =", chi, ")"
     if (active_theory /= THEORY_GR) then
        write(6+215*(i-1),"(A18,ES18.9)")         "    Coupling B =", B_coup
        write(6+215*(i-1),"(A18,ES18.9)")         "   Scalar mass =", sqrt(mphi_r*1.d10/KAPPA)*l_uni
        write(6+215*(i-1),"(A18,ES18.9)")         "     varphi(0) =", sphi_c
        write(6+215*(i-1),"(A18,ES18.9)")         "    varphi_max =", sphi_m
     endif
     write(6+215*(i-1),"(A18,F18.9)")             "        M2/M^3 =", M2
     write(6+215*(i-1),"(A18,F18.9)")             "        S3/M^4 =", S3
     write(6+215*(i-1),"(A18,F18.9)")             "        M4/M^5 =", M4
     write(6+215*(i-1),"(A18,F18.9)")             "           T/W =", T_kin/abs(mass_p - mass + T_kin)
     write(6+215*(i-1),"(A18,F18.9,A4)")          "       Coord R =", r_e*sqrt(KAPPA)/1.d5,"km"
     write(6+215*(i-1),"(A18,F18.9,A4)")          "       Areal R =", r_circ/1.d5,"km"
  end do

  write(unit=*, fmt=*) " "
  write(unit=*, fmt=*) "In code unit:"
  write(*,"(A6,ES18.9,2X,A4,ES18.9,2X,A5,ES18.9,2X,A8,ES18.9)")  &
       "h_c", h_center, "r_e", r_e, "Fmax", Fmax_h, "Omega_e", Omega_e*r_e
  write(unit=*, fmt=*) " ===================================="
  close(221)
  write(unit=*, fmt=*) " "
  write(unit=*, fmt=*) "Completed!"
  call flush(6)
end subroutine print_converged_block


subroutine perform_scalar_burn(target_mphi)
  use para_mod
  use rotation_dispatch, only: call_rotation_solver
  implicit none
  real(8), intent(in) :: target_mphi
  real(8) :: current_mphi
  integer :: burn_iter
  character(100) :: string

  if (active_theory == THEORY_GR) return

  print *, " "
  print *, "scalar burn stage:  (B, mphi, sphi_c, sphi_m)"
  B_coup = B_burn_init
  mphi_r = (mphi_burn_seed / l_uni)**2 * KAPPA / 1.d10
  call call_rotation_solver()
  current_mphi = sqrt(mphi_r*1.d10/KAPPA) * l_uni

  burn_iter = 0
  do while (current_mphi < target_mphi .and. burn_iter < scalar_burn_max_iter)
    call call_rotation_solver()
    B_coup = merge( B_coup * 1.5d0, B_coup * 1.d0, current_mphi > 30.d0)
    if ( sphi_m < 0.4d0 ) then
      mphi_r = mphi_r * 1.1d0
    else
      mphi_r = mphi_r * 3.0d0
    endif
    current_mphi = sqrt(mphi_r*1.d10/KAPPA) * l_uni
    if ( mod(burn_iter,10)== 0 ) write(*,"(i3,A,2(es12.3),A,2(1X,ES18.9))") burn_iter+1, ") ", &
                                  B_coup, current_mphi, "  |", sphi_c, sphi_m
    burn_iter = burn_iter + 1
  end do
  output = .true.; call call_rotation_solver(); output = .false.
  write(*,"(A)") " ", " Solution saved for restart after burning.", " "
  write(string,"(f12.3)") sphi_m
  if (burn_iter >= scalar_burn_max_iter .and. current_mphi < target_mphi) &
  write(unit=*, fmt=*) "scalar burn stage reached iteration limit before hitting target mass."
  write(*,"(A)") " ", merge("*** initial guess is non-scalarized", &
                    "*** Starts with sphi max : "//trim(adjustl(string)), &
                    sphi_m < 1.d-5), " "
end subroutine perform_scalar_burn
