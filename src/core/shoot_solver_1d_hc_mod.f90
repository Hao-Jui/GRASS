
! Active 1D shooting path: textbook Newton in x = log(hc), with
! a finite-difference Jacobian rebuilt at every iteration and an
! Armijo backtracking line search for globalization.

module shoot_solver_1d_types_mod
  use precision_mod, only: wp
  implicit none
  type :: newton_state_1d
    logical :: has_jacobian = .false.
    logical :: has_prev     = .false.
    real(wp) :: J            = 0._wp
    real(wp) :: x_prev       = 0._wp
    real(wp) :: F_prev       = 0._wp
  end type newton_state_1d

  abstract interface
    subroutine evaluation_function_1d(hc, rep, F, rho0, ee)
      import :: wp
      real(wp), intent(in)  :: hc, rep
      real(wp), intent(out) :: F, rho0, ee
    end subroutine evaluation_function_1d
  end interface

end module shoot_solver_1d_types_mod

module shoot_solver_1d_hc_mod
  use precision_mod, only: wp
  use shoot_solver_1d_types_mod, only: newton_state_1d
  implicit none
  ! Adaptive step size cap based on error magnitude
  ! Linear interpolation: cap = CAP_MIN + (CAP_MAX - CAP_MIN) * max(0, 1 - er)
  ! Reference table (er = error, cap = log-space bound, max change = exp(cap)):
  !   er = 1.0  -->  cap = 0.2  -->  max h change = ±22%
  !   er = 0.5  -->  cap = 1.1  -->  max h change = ±3.0×
  !   er = 0.1  -->  cap = 1.82 -->  max h change = ±6.2×
  !   er < 0.01 -->  cap ≈ 2.0  -->  max h change = ±7.4× (full Newton)
  real(wp), parameter :: STEP_CAP_MIN = 0.2_wp   ! tight cap when er ~ 1 (prevent catastrophe)
  real(wp), parameter :: STEP_CAP_MAX = 2.0_wp   ! loose cap near convergence (preserve Newton speed)

contains
  subroutine reset_newton_state_1d(state)
    type(newton_state_1d), intent(inout) :: state
    state%has_jacobian = .false.
    state%has_prev     = .false.
    state%J            = 0._wp
    state%x_prev       = 0._wp
    state%F_prev       = 0._wp
  end subroutine reset_newton_state_1d
  subroutine to_solver_coord_1d(hc, x)
    real(wp), intent(in)  :: hc
    real(wp), intent(out) :: x
    x = log(max(hc, 1.e-12_wp))
  end subroutine to_solver_coord_1d
  subroutine from_solver_coord_1d(x, hc)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: hc
    hc = exp(x)
  end subroutine from_solver_coord_1d
  subroutine clamp_step_1d(delta, er)
    real(wp), intent(inout) :: delta
    real(wp), intent(in)    :: er
    real(wp) :: cap
    ! Adaptive cap: tight when far from solution, loose when close
    cap = STEP_CAP_MIN + (STEP_CAP_MAX - STEP_CAP_MIN) * max(0._wp, 1._wp - er)
    delta = max(-cap, min(delta, cap))
  end subroutine clamp_step_1d

  logical function solve_linear_1d(J, rhs, delta)
    real(wp), intent(in)  :: J, rhs
    real(wp), intent(out) :: delta
    if (abs(J) < 1.e-12_wp) then
      delta = 0._wp
      solve_linear_1d = .false.
    else
      delta = rhs / J
      solve_linear_1d = .true.
    end if
  end function solve_linear_1d
  subroutine broyden_update_1d(state, x, F)
    type(newton_state_1d), intent(inout) :: state
    real(wp), intent(in) :: x, F
    real(wp) :: dx, dF, x_scale, f_scale

    if (.not. state%has_prev) return
    if (.not. state%has_jacobian) return

    dx      = x - state%x_prev
    x_scale = max(1._wp, abs(x), abs(state%x_prev))
    if (abs(dx) <= 1.e2_wp * epsilon(x_scale) * x_scale) return

    dF = F - state%F_prev
    f_scale = max(1._wp, abs(F), abs(state%F_prev))
    if (abs(dF) <= 1.e2_wp * epsilon(f_scale) * f_scale) return

    state%J = state%J + (dF - state%J * dx) / dx
  end subroutine broyden_update_1d
  subroutine commit_state_1d(state, x, F)
    type(newton_state_1d), intent(inout) :: state
    real(wp), intent(in) :: x, F
    state%x_prev = x
    state%F_prev = F
    state%has_prev = .true.
  end subroutine commit_state_1d

  subroutine line_search_1d(x_current, F_current, delta_x, rep, evaluate_func, final_delta, J_est, success)
    ! Armijo backtracking on phi = 0.5*F^2; uses Jacobian estimate when supplied for slope
    use shoot_solver_1d_types_mod, only: evaluation_function_1d
    real(wp), intent(in)    :: x_current, F_current, delta_x, rep
    real(wp), intent(out)   :: final_delta
    procedure(evaluation_function_1d) :: evaluate_func
    real(wp), intent(in), optional :: J_est
    logical, intent(out), optional :: success
    
    integer, parameter :: max_iter = 15
    real(wp), parameter :: TAU = 0.5_wp
    real(wp), parameter :: C1 = 1.e-4
    real(wp)             :: alpha, x_trial, F_trial, hc_trial, rho0_tmp, ee_tmp
    real(wp)             :: hc_base
    real(wp)             :: phi_old, phi_new, slope0
    integer             :: i
    logical             :: ok

    alpha = 1._wp
    final_delta = delta_x
    phi_old = 0.5_wp * F_current**2
    ok = .false.
    call from_solver_coord_1d(x_current, hc_base)
    if (abs(delta_x) < 1.e-12_wp) then
      final_delta = delta_x
      if (present(success)) success = .true.
      return
    end if
    if (present(J_est)) then
      slope0 = F_current * J_est * delta_x
    else
      slope0 = -abs(F_current * delta_x)
    end if
    if (slope0 >= 0._wp) slope0 = -abs(F_current * delta_x)

    do i = 1, max_iter
      x_trial = x_current + alpha * delta_x
      call from_solver_coord_1d(x_trial, hc_trial)
      call evaluate_func(hc_trial, rep, F_trial, rho0_tmp, ee_tmp)

      phi_new = 0.5_wp * F_trial**2
      if (phi_new <= phi_old + c1 * alpha * slope0) then
        ok = .true.
        exit
      end if

      alpha = alpha * tau
    end do
    final_delta = alpha * delta_x
    ! Restore base state to avoid leaving globals at the last trial if caller rejects the step
    call evaluate_func(hc_base, rep, F_trial, rho0_tmp, ee_tmp)
    if (present(success)) success = ok
  end subroutine line_search_1d

end module shoot_solver_1d_hc_mod

module shoot_solver_1d_hc_helpers_mod
  use precision_mod, only: wp
  use analysis_mod, only: mass_radius
  use eos_mod, only: n0_at_h, e_at_h
  use para_mod, only: wp, h_center, r_ratio, Mass, Mass_0, MSUN, M_goal, Mb_goal, FIX1
  use shoot_solver_1d_types_mod, only: newton_state_1d, evaluation_function_1d
  use shoot_solver_1d_hc_mod, only: from_solver_coord_1d
  use rotation_solver_mod, only: rotation_solver
  implicit none
contains  
  subroutine evaluate_solution_1d(hc, rep, F, rho0, ee)
    real(wp), intent(in)  :: hc, rep
    real(wp), intent(out) :: F, rho0, ee
    real(wp) :: devi

    h_center = hc
    r_ratio  = rep

    call rotation_solver
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    if ( trim(FIX1) == "M_goal" ) then
      devi = Mass/MSUN/M_goal - 1._wp
    elseif ( trim(FIX1) == "Mb_goal" ) then
      devi = Mass_0/MSUN/Mb_goal - 1._wp
    else
      stop "evaluate_solution: unknown FIX1"
    endif

    F = devi
  end subroutine evaluate_solution_1d

  subroutine build_jacobian_1d(state, x, F, hc, rep, rho0, ee, reuse_base)
    ! Finite-difference J; optionally reuse the caller’s base residual to save a solve
    type(newton_state_1d), intent(inout) :: state
    real(wp), intent(in)    :: x
    real(wp), intent(inout) :: F
    real(wp), intent(inout) :: hc
    real(wp), intent(in)    :: rep
    real(wp), intent(out)   :: rho0, ee
    logical, intent(in), optional :: reuse_base
    real(wp) :: delta, xp, Fp, rho_tmp, ee_tmp
    real(wp) :: hc_p

    ! Use a larger FD step to avoid flat slopes from numerical noise near root
    delta = max(abs(x), 1._wp) * 1.e-4
    xp = x + delta
    call from_solver_coord_1d(xp, hc_p)
    call evaluate_solution_1d(hc_p, rep, Fp, rho_tmp, ee_tmp)
    state%J = (Fp - F) / delta

    state%has_jacobian = .true.
    state%has_prev     = .false.

    call from_solver_coord_1d(x, hc)
    ! Caller can skip recomputing the base residual if it was just evaluated
    if (.not. (present(reuse_base) .and. reuse_base)) then
      call evaluate_solution_1d(hc, rep, F, rho0, ee)
    end if
  end subroutine build_jacobian_1d

end module shoot_solver_1d_hc_helpers_mod
