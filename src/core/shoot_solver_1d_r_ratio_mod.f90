
! 1D shooting for FIX1 via r_ratio: Newton in x = logit(r_ratio),
! with finite-difference Jacobian, Broyden updates, and Armijo
! backtracking.  Same target as the 1D hc solver (M_goal or Mb_goal)
! but adjusts r_ratio instead of hc.

module shoot_solver_1d_r_ratio_mod
  use shoot_solver_1d_types_mod, only: newton_state_1d
  implicit none
  real(8), parameter :: RR_CAP = 0.5d0

contains
  subroutine to_solver_coord_rp(rr, x)
    real(8), intent(in)  :: rr
    real(8), intent(out) :: x
    real(8) :: rc
    rc = max(1.d-3, min(rr, 1.d0 - 1.d-3))
    x  = log(rc / (1.d0 - rc))
  end subroutine to_solver_coord_rp

  subroutine from_solver_coord_rp(x, rr)
    real(8), intent(in)  :: x
    real(8), intent(out) :: rr
    rr = 1.d0 / (1.d0 + exp(-x))
  end subroutine from_solver_coord_rp

  subroutine clamp_step_rp(delta, er)
    real(8), intent(inout) :: delta
    real(8), intent(in)    :: er
    real(8) :: cap
    cap = min(RR_CAP, 0.1d0 + 0.4d0 * max(0.d0, 1.d0 - er))
    delta = max(-cap, min(delta, cap))
  end subroutine clamp_step_rp

  subroutine line_search_rp(x_current, F_current, delta_x, hc, evaluate_func, final_delta, J_est, success)
    use shoot_solver_1d_types_mod, only: evaluation_function_1d
    real(8), intent(in)    :: x_current, F_current, delta_x, hc
    real(8), intent(out)   :: final_delta
    procedure(evaluation_function_1d) :: evaluate_func
    real(8), intent(in), optional :: J_est
    logical, intent(out), optional :: success

    integer, parameter :: max_iter = 15
    real(8), parameter :: tau = 0.5d0, c1 = 1.d-4
    real(8) :: alpha, x_trial, F_trial, rr_trial, rho0_tmp, ee_tmp
    real(8) :: rr_base, phi_old, phi_new, slope0
    integer :: i
    logical :: ok

    alpha = 1.d0
    final_delta = delta_x
    phi_old = 0.5d0 * F_current**2
    ok = .false.
    call from_solver_coord_rp(x_current, rr_base)

    if (abs(delta_x) < 1.d-12) then
      final_delta = delta_x
      if (present(success)) success = .true.
      return
    end if

    if (present(J_est)) then
      slope0 = F_current * J_est * delta_x
    else
      slope0 = -abs(F_current * delta_x)
    end if
    if (slope0 >= 0.d0) slope0 = -abs(F_current * delta_x)

    do i = 1, max_iter
      x_trial = x_current + alpha * delta_x
      call from_solver_coord_rp(x_trial, rr_trial)
      call evaluate_func(rr_trial, hc, F_trial, rho0_tmp, ee_tmp)

      phi_new = 0.5d0 * F_trial**2
      if (phi_new <= phi_old + c1 * alpha * slope0) then
        ok = .true.
        exit
      end if
      alpha = alpha * tau
    end do

    final_delta = alpha * delta_x
    call evaluate_func(rr_base, hc, F_trial, rho0_tmp, ee_tmp)
    if (present(success)) success = ok
  end subroutine line_search_rp

end module shoot_solver_1d_r_ratio_mod

module shoot_solver_1d_r_ratio_helpers_mod
  use analysis_mod, only: mass_radius
  use eos_mod, only: n0_at_h, e_at_h
  use para_mod, only: h_center, r_ratio, Mass, Mass_0, MSUN, M_goal, Mb_goal, FIX1
  use shoot_solver_1d_types_mod, only: newton_state_1d, evaluation_function_1d
  use shoot_solver_1d_r_ratio_mod, only: from_solver_coord_rp
  use rotation_uniform, only: rotation_solver
  implicit none
contains
  subroutine evaluate_solution_rp(rr, hc, F, rho0, ee)
    real(8), intent(in)  :: rr, hc
    real(8), intent(out) :: F, rho0, ee

    r_ratio  = rr
    h_center = hc

    call rotation_solver
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    if (trim(FIX1) == "M_goal") then
      F = Mass / MSUN / M_goal - 1.d0
    else if (trim(FIX1) == "Mb_goal") then
      F = Mass_0 / MSUN / Mb_goal - 1.d0
    else
      stop "evaluate_solution_rp: unknown FIX1"
    end if
  end subroutine evaluate_solution_rp

  subroutine build_jacobian_rp(state, x, F, rr, hc, rho0, ee, reuse_base)
    type(newton_state_1d), intent(inout) :: state
    real(8), intent(in)    :: x
    real(8), intent(inout) :: F
    real(8), intent(inout) :: rr
    real(8), intent(in)    :: hc
    real(8), intent(out)   :: rho0, ee
    logical, intent(in), optional :: reuse_base
    real(8) :: delta, xp, Fp, rho_tmp, ee_tmp, rr_p

    delta = max(abs(x), 1.d0) * 1.d-4
    xp = x + delta
    call from_solver_coord_rp(xp, rr_p)
    call evaluate_solution_rp(rr_p, hc, Fp, rho_tmp, ee_tmp)
    state%J = (Fp - F) / delta

    state%has_jacobian = .true.
    state%has_prev     = .false.

    call from_solver_coord_rp(x, rr)
    if (.not. (present(reuse_base) .and. reuse_base)) then
      call evaluate_solution_rp(rr, hc, F, rho0, ee)
    end if
  end subroutine build_jacobian_rp

end module shoot_solver_1d_r_ratio_helpers_mod
