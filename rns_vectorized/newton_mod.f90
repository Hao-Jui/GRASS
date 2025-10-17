module shoot_solver_mod
  implicit none
  real(8), parameter :: r_eps = 1.d-8
  real(8), parameter :: max_step = 0.5d0

  type :: newton_state
    logical :: has_jacobian = .false.
    logical :: has_prev     = .false.
    real(8) :: J(2,2)       = 0.d0
    real(8) :: x_prev(2)    = 0.d0
    real(8) :: F_prev(2)    = 0.d0
  end type newton_state
contains
  subroutine reset_newton_state(state)
    type(newton_state), intent(inout) :: state
    state%has_jacobian = .false.
    state%has_prev     = .false.
    state%J            = 0.d0
    state%x_prev       = 0.d0
    state%F_prev       = 0.d0
  end subroutine reset_newton_state

  subroutine to_solver_coords(hc, rep, x)
    real(8), intent(in)  :: hc, rep
    real(8), intent(out) :: x(2)
    real(8) :: rep_clip

    rep_clip = min(max(rep, r_eps), 1.d0 - r_eps)
    x(1) = log(max(hc, 1.d-12))
    x(2) = log(rep_clip / (1.d0 - rep_clip))
  end subroutine to_solver_coords

  subroutine from_solver_coords(x, hc, rep)
    real(8), intent(in)  :: x(2)
    real(8), intent(out) :: hc, rep
    real(8) :: exp_arg

    hc = exp(x(1))
    exp_arg = exp(-x(2))
    rep = 1.d0 / (1.d0 + exp_arg)
    rep = min(max(rep, r_eps), 1.d0 - r_eps)
  end subroutine from_solver_coords

  subroutine clamp_step(delta)
    real(8), intent(inout) :: delta(2)
    delta = max(-max_step, min(delta, max_step))
  end subroutine clamp_step

  logical function solve_linear(J, rhs, delta)
    real(8), intent(in)  :: J(2,2), rhs(2)
    real(8), intent(out) :: delta(2)
    real(8) :: det

    det = J(1,1) * J(2,2) - J(1,2) * J(2,1)
    if (abs(det) < 1.d-12) then
      delta = 0.d0
      solve_linear = .false.
    else
      delta(1) = ( rhs(1) * J(2,2) - J(1,2) * rhs(2) ) / det
      delta(2) = ( J(1,1) * rhs(2) - rhs(1) * J(2,1) ) / det
      solve_linear = .true.
    endif
  end function solve_linear

  subroutine broyden_update(state, x, F)
    type(newton_state), intent(inout) :: state
    real(8), intent(in) :: x(2), F(2)
    real(8) :: dx(2), dF(2), denom, Jdx(2)

    if (.not. state%has_prev) return
    if (.not. state%has_jacobian) return

    dx = x - state%x_prev
    denom = dot_product(dx, dx)
    if (denom <= 1.d-16) return

    dF  = F - state%F_prev
    Jdx = matmul(state%J, dx)

    state%J = state%J + matmul(reshape(dF - Jdx, (/2,1/)), reshape(dx, (/1,2/))) / denom
  end subroutine broyden_update

  subroutine commit_state(state, x, F)
    type(newton_state), intent(inout) :: state
    real(8), intent(in) :: x(2), F(2)

    state%x_prev = x
    state%F_prev = F
    state%has_prev = .true.
  end subroutine commit_state

end module shoot_solver_mod

module shoot_newton_helpers
  use para_mod
  use shoot_solver_mod
  implicit none
contains
  subroutine evaluate_solution(hc, rep, F, rho0, ee, er)
    real(8), intent(in)  :: hc, rep
    real(8), intent(out) :: F(2), rho0, ee, er
    real(8) :: n0_at_h, e_at_h
    external :: n0_at_h, e_at_h
    external :: call_rotation_solver, mass_radius
    real(8) :: deviA, deviB

    r_ratio  = rep
    h_center = hc

    call call_rotation_solver
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    deviA = Mass_0/MSUN/Mb_goal - 1.d0
    deviB = ( Omega_e - Omega_K / (C/sqrt(kappa)) ) / 5.d0

    F(1) = deviA
    F(2) = deviB
    er   = abs(deviA) + abs(deviB)
  end subroutine evaluate_solution

  subroutine build_jacobian(state, x, F, hc, rep, rho0, ee, er)
    type(newton_state), intent(inout) :: state
    real(8), intent(in)    :: x(2)
    real(8), intent(inout) :: F(2)
    real(8), intent(inout) :: hc, rep
    real(8), intent(out)   :: rho0, ee, er
    real(8) :: delta(2), xp(2), Fp(2), rho_tmp, ee_tmp, er_tmp
    real(8) :: hc_p, rep_p
    integer :: i

    delta = max(0.05d0, 0.2d0*abs(x))
    do i = 1, 2
      xp = x
      xp(i) = xp(i) + delta(i)
      call from_solver_coords(xp, hc_p, rep_p)
      call evaluate_solution(hc_p, rep_p, Fp, rho_tmp, ee_tmp, er_tmp)
      state%J(:,i) = (Fp - F) / delta(i)
    enddo

    state%has_jacobian = .true.
    state%has_prev     = .false.

    call from_solver_coords(x, hc, rep)
    call evaluate_solution(hc, rep, F, rho0, ee, er)
  end subroutine build_jacobian

end module shoot_newton_helpers

