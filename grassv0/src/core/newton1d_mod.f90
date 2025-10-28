module shoot_solver_mod_1d
  implicit none
  real(8), parameter :: max_step_1d = 0.5d0

  type :: newton_state_1d
    logical :: has_jacobian = .false.
    logical :: has_prev     = .false.
    real(8) :: J            = 0.d0
    real(8) :: x_prev       = 0.d0
    real(8) :: F_prev       = 0.d0
  end type newton_state_1d

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
    real(8) :: dx, dF

    if (.not. state%has_prev) return
    if (.not. state%has_jacobian) return

    dx = x - state%x_prev
    if (abs(dx) <= 1.d-16) return

    dF = F - state%F_prev
    state%J = state%J + (dF - state%J * dx) / dx
  end subroutine broyden_update_1d
  subroutine commit_state_1d(state, x, F)
    type(newton_state_1d), intent(inout) :: state
    real(8), intent(in) :: x, F
    state%x_prev = x
    state%F_prev = F
    state%has_prev = .true.
  end subroutine commit_state_1d

end module shoot_solver_mod_1d

module shoot_newton_helpers_1d
  use para_mod
  use rotation_dispatch, only: call_rotation_solver
  use shoot_solver_mod_1d
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

    call call_rotation_solver
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    devi = mass_0/MSUN - Mb_goal

    F = devi
  end subroutine evaluate_solution_1d

  subroutine build_jacobian_1d(state, x, F, hc, rep, rho0, ee)
    type(newton_state_1d), intent(inout) :: state
    real(8), intent(in)    :: x
    real(8), intent(inout) :: F
    real(8), intent(inout) :: hc
    real(8), intent(in)    :: rep
    real(8), intent(out)   :: rho0, ee
    real(8) :: delta, xp, Fp, rho_tmp, ee_tmp
    real(8) :: hc_p

    delta = max(0.05d0, 0.2d0*abs(x))
    xp = x + delta
    call from_solver_coord_1d(xp, hc_p)
    call evaluate_solution_1d(hc_p, rep, Fp, rho_tmp, ee_tmp)
    state%J = (Fp - F) / delta

    state%has_jacobian = .true.
    state%has_prev     = .false.

    call from_solver_coord_1d(x, hc)
    write(*,*) hc_p, hc
    call evaluate_solution_1d(hc, rep, F, rho0, ee)
  end subroutine build_jacobian_1d

end module shoot_newton_helpers_1d
