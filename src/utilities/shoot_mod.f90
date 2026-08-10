module shoot_mod
  use para_mod, only: wp
  use shoot_solver_2d_mod, only: newton_state, init_newton_state, reset_newton_state, &
                                 solve_linear, clamp_step, from_solver_coords, to_solver_coords, &
                                 broyden_update, commit_state
  use shoot_solver_1d_types_mod, only: newton_state_1d
  use shoot_solver_1d_hc_mod, only: reset_newton_state_1d, solve_linear_1d, &
    line_search_1d, from_solver_coord_1d, to_solver_coord_1d, &
    clamp_step_1d, broyden_update_1d, commit_state_1d
  use shoot_solver_1d_r_ratio_mod, only: from_solver_coord_rp, to_solver_coord_rp, &
                                         clamp_step_rp, line_search_rp
  implicit none
  private
  public :: shoot_v2

  real(wp) :: er, rho0, ee
  real(wp) :: F(2), x(2), delta_x(2), rhs(2), h_new, r_new
  real(wp) :: F1d, x1d, delta_x1d, rhs1d
  real(wp) :: F_rp, x_rp, delta_x_rp, rhs_rp, r_new_rp
  type(newton_state)    :: solver_state
  type(newton_state_1d) :: solver_state_1d
  type(newton_state_1d) :: solver_state_rp
  logical :: step_ok, need_cycle
  real(wp) :: prev_er1d, prev_er2d, prev_er_rp

contains

  subroutine shoot_v2
    use para_mod, only: shooting, SHOOT_FIX1_HC
    use starting_model_mod, only: initialize_starting_model
    integer :: iteration_cap

    call initialize_starting_model()
    call init_newton_state(solver_state, 2)

    write(unit=*, fmt=*) " "
    iteration_cap = 1

    if (iteration_cap == 1) then
      call shoot_single()
    else
      call shoot_sequence(iteration_cap)
    end if
  end subroutine shoot_v2

  subroutine shoot_sequence(iteration_cap)
    use para_mod, only: h_center, r_ratio, Omega_e, KAPPA, C, pi, &
                        shooting, SHOOT_FIX1_HC, SHOOT_FIX1_RP
    integer, intent(in) :: iteration_cap
    integer :: i_idx

    do i_idx = 1, iteration_cap
      call shoot_single()
      call output_seq()

      back_bending_project: block
        real(wp) :: Oe
        Oe = Omega_e * (C/sqrt(kappa)) / 2.0_wp / pi
        h_center = merge(h_center + 0.002_wp, h_center + 0.005_wp, Oe < 350.0_wp) ! spin it down
        if (shooting == SHOOT_FIX1_HC) then
          r_ratio = min(r_ratio, 0.95_wp)
        end if
        shooting = SHOOT_FIX1_RP
      end block back_bending_project
    end do
  end subroutine shoot_sequence

  subroutine shoot_single()
    use analysis_mod, only: solution_properties
    use eos_mod, only: n0_at_h, e_at_h
    use para_mod, only: h_center, r_ratio, output, shooting, &
                        SHOOT_FIX1_HC, SHOOT_FIX1_RP, SHOOT_2D, &
                        accuracy, n_of_relaxation_steps, start, finish
    use rotation_solver_mod, only: rotation_solver
    use miscellaneous_mod, only: print_converged_block
    integer :: it
    real(wp) :: er_best
    integer  :: n_stall
    integer, parameter :: STALL_LIMIT = 8
    real(wp), parameter :: STALL_TOL = 1.e-6_wp

    select case (shooting)
    case (SHOOT_FIX1_HC)
      call reset_newton_state_1d(solver_state_1d); prev_er1d = huge(1._wp)
    case (SHOOT_FIX1_RP)
      call reset_newton_state_1d(solver_state_rp); prev_er_rp = huge(1._wp)
    case (SHOOT_2D)
      call reset_newton_state(solver_state); prev_er2d = huge(1._wp)
    end select

    it = 1
    er = 1.e99_wp
    er_best = huge(1._wp)
    n_stall = 0
    do
      call cpu_time(start)
      n_of_relaxation_steps = 0; call evaluate_and_update
      call cpu_time(finish)
      write(*,"(A,ES15.6,5X,A,I5,5X,A,F8.3)") &
            "Error: ", er, &
            "Iterations: ", n_of_relaxation_steps, &
            "Elapsed Time [s]: ", finish-start

      if (er <= accuracy .or. output) exit

      if (er < er_best * (1._wp - STALL_TOL)) then
        er_best = er
        n_stall = 0
      else
        n_stall = n_stall + 1
      end if
      if (n_stall >= STALL_LIMIT) then
        write(*,"(A,ES10.3,A)") "  Shooting stalled (er=", er, "), accepting solution."
        stop
      end if

      call apply_newton_step
      if (need_cycle) cycle

      it = it + 1
      if (it == 1000) stop "Required bulk properties cannot be reached."
    end do

    output = .true.
    call rotation_solver; call solution_properties
    output = .false.

    rho0 = n0_at_h(h_center)
    ee   = e_at_h(h_center)
    write(unit=*, fmt=*) " "
    write(*,'(A,2es18.9)') "converged rho0, ee:", rho0, ee
    call print_converged_block(rho0, ee)
  end subroutine shoot_single

  subroutine output_seq()
    use para_mod, only: rho_uni, ang_mom, h_center, r_ratio, &
                        mass, mass_0, mass_p, chi, T_kin, r_circ, &
                        Omega_e, KAPPA, C, KSCALE, MB, MSUN, pi, n_sat, &
                        I_inertia, Love2, M2, M4, S3, &
                        eos_file, sound_speed
    use eos_mod, only: p_at_e
    character(len=1024) :: filename
    real(wp) :: Q_bar, T_over_W, pp, traceT
    integer  :: unit, ios
    pp   = p_at_e(ee)
    traceT  = 3.0_wp*pp*1.80171810e-39_wp/KSCALE - ee*rho_uni/(C * C * KSCALE)
    Q_bar = merge(-1._wp, M2/chi**2, chi < 1.e-30_wp)
    T_over_W=merge(-1._wp, T_kin/abs(Mass_p - Mass + T_kin), chi < 1.e-30_wp)

    write(filename, '(A, A, A, F0.2, A)') &
      "/Users/horay/Data4Projects/HT/Seq_", trim(eos_file), "_M", mass_0/MSUN, ".dat"

    open(newunit=unit, file=trim(filename), access='append', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,I0,A)') 'ERROR: Cannot open output file. IOSTAT = ', ios, trim(filename)
    end if
    write(unit,"(99es18.9e3)") &
            ee/(C * C * KSCALE), rho0*MB/n_sat,    & ! 1-2
            h_center, traceT, sound_speed(1),      & ! 3-5
            Mass/MSUN, Mass_0/MSUN,                & ! 6-7
            I_inertia, Love2, Q_bar,               & ! 8-10
            M2, S3, M4, chi, T_over_W, ang_mom,    & ! 11-16
            Omega_e*(C/sqrt(kappa)), r_ratio,      & ! 17-18
            r_circ/1e5_wp
    close(unit)
  end subroutine output_seq

  subroutine damp_2d_step_near_spherical(x_current, delta_trial, rep_current)
    use shoot_solver_2d_mod, only: rep_map_scale, r_min_ratio, r_eps
    real(wp), intent(in) :: x_current(2), rep_current
    real(wp), intent(inout) :: delta_trial(2)
    real(wp), parameter :: NEAR_SPHERICAL_RATIO = 0.90_wp
    real(wp), parameter :: MAX_GAP_FRACTION = 0.50_wp
    real(wp) :: rep_cap, x2_limit

    if (rep_current < near_spherical_ratio) return

    rep_cap = min(rep_current + max_gap_fraction * (1.0_wp - rep_current), 1.0_wp - r_eps)
    x2_limit = log((1.0_wp - rep_cap) / (rep_cap - r_min_ratio)) / rep_map_scale

    if (x_current(2) + delta_trial(2) < x2_limit) then
      delta_trial(2) = x2_limit - x_current(2)
    end if
  end subroutine damp_2d_step_near_spherical

  subroutine apply_newton_step
    use para_mod, only: h_center, r_ratio, shooting, &
                        SHOOT_FIX1_HC, SHOOT_FIX1_RP, SHOOT_2D, accuracy
    use shoot_solver_2d_helpers_mod, only: evaluate_solution
    use shoot_solver_1d_hc_helpers_mod, only: evaluate_solution_1d
    use shoot_solver_1d_r_ratio_helpers_mod, only: evaluate_solution_rp

    need_cycle = .false.
    select case (shooting)
    case (SHOOT_FIX1_HC)
      rhs1d = -F1d
      step_ok = solve_linear_1d(solver_state_1d%J, rhs1d, delta_x1d)
      if (.not. step_ok) then
        call reset_newton_state_1d(solver_state_1d)
        need_cycle = .true.; return
      endif
      call clamp_step_1d(delta_x1d, er)

      if (er > 1.e-4_wp) then
        call line_search_1d(x1d, F1d, delta_x1d, r_ratio, evaluate_solution_1d, delta_x1d, &
                            J_est=solver_state_1d%J, success=step_ok)
        if (.not. step_ok) then
          call reset_newton_state_1d(solver_state_1d)
          need_cycle = .true.; return
        end if
      end if

      x1d = x1d + delta_x1d
      call from_solver_coord_1d(x1d, h_new)
      h_center = h_new

    case (SHOOT_FIX1_RP)
      rhs_rp = -F_rp
      step_ok = solve_linear_1d(solver_state_rp%J, rhs_rp, delta_x_rp)
      if (.not. step_ok) then
        call reset_newton_state_1d(solver_state_rp)
        need_cycle = .true.; return
      endif
      call clamp_step_rp(delta_x_rp, er)

      if (er > 1.e-4_wp) then
        call line_search_rp(x_rp, F_rp, delta_x_rp, h_center, evaluate_solution_rp, delta_x_rp, &
                              J_est=solver_state_rp%J, success=step_ok)
        if (.not. step_ok) then
          call reset_newton_state_1d(solver_state_rp)
          need_cycle = .true.; return
        end if
      end if

      x_rp = x_rp + delta_x_rp
      call from_solver_coord_rp(x_rp, r_new_rp)
      r_ratio = r_new_rp

    case (SHOOT_2D)
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
      call damp_2d_step_near_spherical(x, delta_x, r_ratio)
      x = x + delta_x
      call from_solver_coords(x, h_new, r_new)
      h_center = h_new
      r_ratio  = r_new
      solver_state%has_jacobian = .true.
    end select
  end subroutine apply_newton_step

  subroutine evaluate_and_update
    use para_mod, only: h_center, r_ratio, shooting, &
                        SHOOT_FIX1_HC, SHOOT_FIX1_RP, SHOOT_2D
    use shoot_solver_2d_helpers_mod, only: evaluate_solution, build_jacobian
    use shoot_solver_1d_hc_helpers_mod, only: evaluate_solution_1d, build_jacobian_1d
    use shoot_solver_1d_r_ratio_helpers_mod, only: evaluate_solution_rp, build_jacobian_rp

    select case (shooting)
    case (SHOOT_FIX1_HC)
      call evaluate_solution_1d(h_center, r_ratio, F1d, rho0, ee)
      call to_solver_coord_1d(h_center, x1d)
      er = abs(F1d)
      if (solver_state_1d%has_jacobian .and. er > prev_er1d) &
        call reset_newton_state_1d(solver_state_1d)
      prev_er1d = er
      if (solver_state_1d%has_jacobian) then
        call broyden_update_1d(solver_state_1d, x1d, F1d)
      else
        call build_jacobian_1d(solver_state_1d, x1d, F1d, h_center, r_ratio, rho0, ee, reuse_base=.true.)
        call to_solver_coord_1d(h_center, x1d)
      end if
      call commit_state_1d(solver_state_1d, x1d, F1d)

    case (SHOOT_FIX1_RP)
      call evaluate_solution_rp(r_ratio, h_center, F_rp, rho0, ee)
      call to_solver_coord_rp(r_ratio, x_rp)
      er = abs(F_rp)
      if (solver_state_rp%has_jacobian .and. er > prev_er_rp) &
        call reset_newton_state_1d(solver_state_rp)
      prev_er_rp = er
      if (solver_state_rp%has_jacobian) then
        call broyden_update_1d(solver_state_rp, x_rp, F_rp)
      else
        call build_jacobian_rp(solver_state_rp, x_rp, F_rp, r_ratio, h_center, rho0, ee, reuse_base=.true.)
        call to_solver_coord_rp(r_ratio, x_rp)
      end if
      call commit_state_1d(solver_state_rp, x_rp, F_rp)

    case (SHOOT_2D)
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
    end select
  end subroutine evaluate_and_update

end module shoot_mod

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
  write(*,"(A10,ES18.9,A4)")    "OMG_c :", omega_c/(2._wp*pi)*(C/sqrt(kappa)),"Hz"
  write(*,"(A10,ES18.9,A4)")    "Mass  :", mass/MSUN,"M_o"
  write(*,"(A10,ES18.9,A4)")    "M_0   :", mass_0/MSUN,"M_o"
  write(*,"(A10,ES18.9)")       "J     :", ang_mom
  write(*,"(A10,ES18.9)")       "chi   :", chi
  write(*,"(A10,ES18.9)")       "T/W   :", T_kin/abs(mass_p - mass + T_kin)
  write(*,"(A10,ES18.9,A4)")    "R_is  :", r_e*sqrt(KAPPA)/1.e5_wp,"km"
  write(*,"(A10,ES18.9,A4)")    "R_cir :", r_circ/1.e5_wp,"km"
  write(*,"(A10,ES18.9,A4)")    "r_e/M :", (r_circ/1.e5_wp)/(mass/MSUN*1.4769994423016508_wp)
  if (.not. active_theory == THEORY_GR) then
    write(*,"(A10,2ES18.9)")    "Om_K/e:", Omega_K/(2._wp*pi), Omega_e/(2._wp*pi)*(C/sqrt(kappa))
  end if
  write(*,"(A10,3ES18.9)")      "er    :", er
  write(*,'(1X,A)') repeat('=', 36)
  write(*, *) " "
end subroutine print_iter_status
