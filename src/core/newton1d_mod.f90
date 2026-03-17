
!Subroutine	                                    Purpose
!reset_newton_state_1d	                        Reinitializes solver state
!to_solver_coord_1d / from_solver_coord_1d	    Transforms between physical (hc) and solver (x = log(hc)) coordinates
!clamp_step_1d	                                Trust-region style step limiting (±1.0 max)
!solve_linear_1d	                              Solves J * delta = rhs with singularity check
!broyden_update_1d	                            Rank-1 Broyden Jacobian update (secant method style)
!commit_state_1d	                              Saves current iterate for next Broyden update
!line_search_1d

module newton_types_mod
  implicit none
  type :: newton_state_1d
    logical :: has_jacobian = .false.
    logical :: has_prev     = .false.
    real(8) :: J            = 0.d0
    real(8) :: x_prev       = 0.d0
    real(8) :: F_prev       = 0.d0
  end type newton_state_1d

  abstract interface
    subroutine evaluation_function_1d(hc, rep, F, rho0, ee)
      real(8), intent(in)  :: hc, rep
      real(8), intent(out) :: F, rho0, ee
    end subroutine evaluation_function_1d
  end interface

end module newton_types_mod

module shoot_solver_mod_1d
  use newton_types_mod
  implicit none
  real(8), parameter :: max_step_1d = 1.d0  ! trust-region style clamp on log(h) step; smaller to be less aggressive

contains
  subroutine reset_newton_state_1d(state)
    type(newton_state_1d), intent(inout) :: state
    state%has_jacobian = .false.
    state%has_prev     = .false.
    state%J            = 0.d0
    state%x_prev       = 0.d0
    state%F_prev       = 0.d0
  end subroutine reset_newton_state_1d
  subroutine to_solver_coord_1d(hc, x)
    real(8), intent(in)  :: hc
    real(8), intent(out) :: x
    x = log(max(hc, 1.d-12))
  end subroutine to_solver_coord_1d
  subroutine from_solver_coord_1d(x, hc)
    real(8), intent(in)  :: x
    real(8), intent(out) :: hc
    hc = exp(x)
  end subroutine from_solver_coord_1d
  subroutine clamp_step_1d(delta)
    real(8), intent(inout) :: delta
    delta = max(-max_step_1d, min(delta, max_step_1d))
  end subroutine clamp_step_1d

  logical function solve_linear_1d(J, rhs, delta)
    real(8), intent(in)  :: J, rhs
    real(8), intent(out) :: delta
    if (abs(J) < 1.d-12) then
      delta = 0.d0
      solve_linear_1d = .false.
    else
      delta = rhs / J
      solve_linear_1d = .true.
    end if
  end function solve_linear_1d
  subroutine broyden_update_1d(state, x, F)
    type(newton_state_1d), intent(inout) :: state
    real(8), intent(in) :: x, F
    real(8) :: dx, dF, x_scale, f_scale

    if (.not. state%has_prev) return
    if (.not. state%has_jacobian) return

    dx      = x - state%x_prev
    x_scale = max(1.d0, abs(x), abs(state%x_prev))
    if (abs(dx) <= 1.d2 * epsilon(x_scale) * x_scale) return

    dF = F - state%F_prev
    f_scale = max(1.d0, abs(F), abs(state%F_prev))
    if (abs(dF) <= 1.d2 * epsilon(f_scale) * f_scale) return

    state%J = state%J + (dF - state%J * dx) / dx
  end subroutine broyden_update_1d
  subroutine commit_state_1d(state, x, F)
    type(newton_state_1d), intent(inout) :: state
    real(8), intent(in) :: x, F
    state%x_prev = x
    state%F_prev = F
    state%has_prev = .true.
  end subroutine commit_state_1d

  subroutine line_search_1d(x_current, F_current, delta_x, rep, evaluate_func, final_delta, J_est, success)
    ! Armijo backtracking on phi = 0.5*F^2; uses Jacobian estimate when supplied for slope
    use newton_types_mod
    real(8), intent(in)    :: x_current, F_current, delta_x, rep
    real(8), intent(out)   :: final_delta
    procedure(evaluation_function_1d) :: evaluate_func
    real(8), intent(in), optional :: J_est
    logical, intent(out), optional :: success
    
    integer, parameter :: max_iter = 15
    real(8), parameter :: tau = 0.5d0
    real(8), parameter :: c1 = 1.d-4
    real(8)             :: alpha, x_trial, F_trial, hc_trial, rho0_tmp, ee_tmp
    real(8)             :: hc_base
    real(8)             :: phi_old, phi_new, slope0
    integer             :: i
    logical             :: ok

    alpha = 1.d0
    final_delta = delta_x
    phi_old = 0.5d0 * F_current**2
    ok = .false.
    call from_solver_coord_1d(x_current, hc_base)
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
      call from_solver_coord_1d(x_trial, hc_trial)
      call evaluate_func(hc_trial, rep, F_trial, rho0_tmp, ee_tmp)

      phi_new = 0.5d0 * F_trial**2
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

end module shoot_solver_mod_1d

module shoot_newton_helpers_1d
  use para_mod, only: h_center, r_ratio, Mass, Mass_0, MSUN, M_goal, Mb_goal, FIX1
  use newton_types_mod
  use shoot_solver_mod_1d, only: from_solver_coord_1d
  use rotation_uniform, only: rotation_solver
  implicit none
contains  
  subroutine evaluate_solution_1d(hc, rep, F, rho0, ee)
    real(8), intent(in)  :: hc, rep
    real(8), intent(out) :: F, rho0, ee
    real(8) :: n0_at_h, e_at_h
    external :: n0_at_h, e_at_h
    external :: mass_radius
    real(8) :: devi

    h_center = hc
    r_ratio  = rep

    call rotation_solver
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    if ( trim(FIX1) == "M_goal" ) then
      devi = Mass/MSUN/M_goal - 1.d0
    elseif ( trim(FIX1) == "Mb_goal" ) then
      devi = Mass_0/MSUN/Mb_goal - 1.d0
    else
      stop "evaluate_solution: unknown FIX1"
    endif

    F = devi
  end subroutine evaluate_solution_1d

  subroutine build_jacobian_1d(state, x, F, hc, rep, rho0, ee, reuse_base)
    ! Finite-difference J; optionally reuse the caller’s base residual to save a solve
    type(newton_state_1d), intent(inout) :: state
    real(8), intent(in)    :: x
    real(8), intent(inout) :: F
    real(8), intent(inout) :: hc
    real(8), intent(in)    :: rep
    real(8), intent(out)   :: rho0, ee
    logical, intent(in), optional :: reuse_base
    real(8) :: delta, xp, Fp, rho_tmp, ee_tmp
    real(8) :: hc_p

    ! Use a larger FD step to avoid flat slopes from numerical noise near root
    delta = max(abs(x), 1.d0) * 1.d-4
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

end module shoot_newton_helpers_1d


