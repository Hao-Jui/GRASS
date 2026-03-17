subroutine shoot_v2
  use para_mod, only: wp, h_center, r_ratio, mass, mass_0, mass_p, chi, chi_goal, &
                      Omega_c, KAPPA, C, KSCALE, MB, MSUN, pi, n_sat, &
                      I_inertia, Love2, M2, M4, S3, T_kin, &
                      mphi_goal, B_goal, eos_file, sound_speed, r_circ, &
                      accuracy, output, use_shoot_1d, start, finish, &
                      n_of_relaxation_steps
  use shoot_solver_mod
  use shoot_solver_mod_1d
  use rotation_uniform,  only: rotation_solver
  use starting_model_mod, only: initialize_starting_model
  use miscellaneous_mod, only: log_kepler_sequence, print_converged_block
  use shoot_newton_helpers, only: evaluate_solution, build_jacobian
  use shoot_newton_helpers_1d, only: evaluate_solution_1d, build_jacobian_1d
  implicit none
  real(wp) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it, i_idx, iteration_cap
  real(wp) :: er, rho0, ee
  real(wp) :: F(2), x(2), delta_x(2), rhs(2), h_new, r_new
  real(wp) :: F1d, x1d, delta_x1d, rhs1d
  type(newton_state)    :: solver_state
  type(newton_state_1d) :: solver_state_1d
  logical :: step_ok, need_cycle
  real(wp) :: prev_er1d, prev_er2d
  external :: p_at_e, h_at_p, n0_at_h, e_at_h

  call initialize_starting_model(p_at_e, h_at_p, n0_at_h, e_at_h)
  ! Always initialize 2D state since the run can switch from 1D to 2D mid-sequence.
  call init_newton_state(solver_state, 2)

  write(unit=*, fmt=*) " "
  chi_goal = 0.e0_wp
  iteration_cap = 1

  do i_idx = 1, iteration_cap
    if (use_shoot_1d) then
      call reset_newton_state_1d(solver_state_1d); prev_er1d = huge(1.e0_wp)
    else
      call reset_newton_state(solver_state); prev_er2d = huge(1.e0_wp)
    endif

    ! ---------------------------------------------------------------
    ! Shooting stellar parameters
    ! ---------------------------------------------------------------
    it = 1
    er = 1.e99_wp
    do
      call cpu_time(start)
      n_of_relaxation_steps = 0; call evaluate_and_update
      call cpu_time(finish)
      write(*,"(A,ES15.6,5X,A,I5,5X,A,F8.3)") &
            "Error: ", er, &
            "Iterations: ", n_of_relaxation_steps, &
            "Elapsed Time [s]: ", finish-start

      if (er <= accuracy .or. output) exit

      call apply_newton_step
      if (need_cycle) cycle
      
      it = it + 1
      if (it == 1000) stop "Required bulk properties cannot be reached."
    end do

    ! ---------------------------------------------------------------
    ! Recording the desired solution
    ! ---------------------------------------------------------------
    output = .true.
    call rotation_solver; call solution_properties
    output = .false.

    rho0 = n0_at_h(h_center)
    ee   = e_at_h(h_center)
    write(unit=*, fmt=*) " "
    write(*,'(A,2es18.9)') "converged rho0, ee:", rho0, ee
    call print_converged_block(rho0, ee)

    call output_seq()
    chi_goal = 0.2e0_wp + dble(i_idx-1) * 0.05e0_wp
    if (use_shoot_1d) then
      ! Avoid entering 2D Newton from the quasi-spherical limit (r_ratio ~ 1),
      ! where the rotation update treats the model as non-rotating.
      r_ratio = min(r_ratio, 0.9e0_wp)
    end if
    use_shoot_1d = .false.
    !call log_kepler_sequence()
  end do ! looping models

contains
  subroutine output_seq()
    real(wp) :: Q_bar, T_over_W
    integer  :: unit
    character(16) :: fil1, fil2
    Q_bar = merge(-1.d0, M2/chi**2, chi < 1.e-30_wp)
    T_over_W=merge(-1.d0, T_kin/abs(Mass_p - Mass + T_kin), chi < 1.e-30_wp)
    write(fil1,"(f10.1)") mphi_goal
    write(fil2,"(es12.1)") B_goal
    open(newunit=unit,file="/Users/horay/Data4Projects/HT/Seq_"//trim(adjustl(eos_file))//".dat",access='append')
    write(unit,"(99es18.9e3)") ee/(C * C * KSCALE), rho0*MB/n_sat, h_center, sound_speed(1), & ! 1-4
            Mass/MSUN, Mass_0/MSUN, r_circ/1e5, & ! 5-7
            0.d0, 0.d0, & ! 8, 9
            I_inertia, Love2, Q_bar, & ! 10-12
            T_over_W, M2, S3, M4, chi, omega_c/(2.e0_wp*pi)*(C/sqrt(kappa)) ! 13-18
    close(unit)
  end subroutine output_seq

  subroutine apply_newton_step
    need_cycle = .false.
    if (use_shoot_1d) then
      rhs1d = -F1d
      step_ok = solve_linear_1d(solver_state_1d%J, rhs1d, delta_x1d)
      if (.not. step_ok) then
        call reset_newton_state_1d(solver_state_1d)
        need_cycle = .true.; return
      endif

      if (solver_state_1d%has_jacobian .and. abs(delta_x1d) < 1.e-6_wp .and. abs(F1d) > accuracy*1.e4_wp) then
        call reset_newton_state_1d(solver_state_1d)
        need_cycle = .true.; return
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
        need_cycle = .true.; return
      endif

      if (solver_state%has_jacobian .and. maxval(abs(delta_x(1:2))) < 1.e-5_wp &
          .and. maxval(abs(F)) > accuracy*1.e4_wp) then
        call reset_newton_state(solver_state)
        write(*,*) "Jacobian refreshed"
        need_cycle = .true.; return
      end if

      call clamp_step(delta_x)
      x = x + delta_x
      call from_solver_coords(x, h_new, r_new)
      h_center = h_new
      r_ratio  = r_new
      solver_state%has_jacobian = .true.
    end if
  end subroutine apply_newton_step

  subroutine evaluate_and_update
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
  end subroutine evaluate_and_update

end subroutine shoot_v2

subroutine print_iter_status(it, rho0, ee, er)
  use para_mod, only: wp, active_theory, THEORY_GR, &
                      B_coup, mphi_r, sphi_c, l_uni, &
                      KAPPA, C, KSCALE, rho_uni, MB, MSUN, pi, &
                      r_ratio, Omega_c, Omega_K, Omega_e, &
                      mass, mass_0, mass_p, ang_mom, chi, T_kin, r_e, r_circ
  implicit none
  integer, intent(in) :: it
  real(wp), intent(in) :: rho0, ee, er

  write(*,'(1X,A)') repeat('=', 36)
  write(*,"(A10,I6)")           " iter :", it/10
  if (active_theory /= THEORY_GR) then
    write(*,"(A10,ES18.9)")     " Bcoup:", B_coup
    write(*,"(A10,ES18.9)")     " mphi :", sqrt(mphi_r*1.e10_wp/KAPPA)*l_uni
    write(*,"(A10,ES18.9)")     " sphi :", sphi_c
  end if
  write(*,"(A10,ES18.9,A10)")    "rho_c :", rho0*MB,"g/cm^3"
  write(*,"(A10,ES18.9,A10,ES18.9)") "e_c   :", ee/(C*C*KSCALE),"g/cm^3", ee/(C*C*KSCALE)*rho_uni
  write(*,"(A10,ES18.9)")       "rp/re :", r_ratio
  write(*,"(A10,ES18.9,A4)")    "OMG_c :", omega_c/(2.e0_wp*pi)*(C/sqrt(kappa)),"Hz"
  write(*,"(A10,ES18.9,A4)")    "Mass  :", mass/MSUN,"M_o"
  write(*,"(A10,ES18.9,A4)")    "M_0   :", mass_0/MSUN,"M_o"
  write(*,"(A10,ES18.9)")       "J     :", ang_mom
  write(*,"(A10,ES18.9)")       "chi   :", chi
  write(*,"(A10,ES18.9)")       "T/W   :", T_kin/abs(mass_p - mass + T_kin)
  write(*,"(A10,ES18.9,A4)")    "R_is  :", r_e*sqrt(KAPPA)/1.e5_wp,"km"
  write(*,"(A10,ES18.9,A4)")    "R_cir :", r_circ/1.e5_wp,"km"
  write(*,"(A10,ES18.9,A4)")    "r_e/M :", (r_circ/1.e5_wp)/(mass/MSUN*1.4769994423016508e0_wp)
  if (.not. active_theory == THEORY_GR) then
    write(*,"(A10,2ES18.9)")    "Om_K/e:", Omega_K/(2.e0_wp*pi), Omega_e/(2.e0_wp*pi)*(C/sqrt(kappa))
  end if
  write(*,"(A10,3ES18.9)")      "er    :", er
  write(*,'(1X,A)') repeat('=', 36)
  write(*, *) " "
end subroutine print_iter_status
