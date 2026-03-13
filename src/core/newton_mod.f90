! uses the Newton-Raphson method (exact Jacobian via finite difference) 
! to get started and then switches to the more computationally efficient 
! Broyden's method (Jacobian updated algebraically) for subsequent steps, 
! which is common when function evaluation is expensive.
module shoot_solver_mod
  use precision_mod, only: wp
  implicit none
  real(wp), parameter :: r_eps = 1.e-8_wp
  real(wp), parameter :: r_min_ratio = 0.35e0_wp
  real(wp), parameter :: max_step = 0.5e0_wp
  real(wp), parameter :: rep_map_scale = 6.e0_wp

  type, public :: newton_state
    logical :: has_jacobian = .false.
    logical :: has_prev     = .false.
    real(wp), allocatable :: J(:,:)
    real(wp), allocatable :: x_prev(:)
    real(wp), allocatable :: F_prev(:)
  end type newton_state

  abstract interface
    subroutine evaluation_function(hc, rep, F, rho0, ee, er)
      import :: wp
      real(wp), intent(in)  :: hc, rep
      real(wp), intent(out) :: F(2), rho0, ee, er
    end subroutine evaluation_function
  end interface

contains
  subroutine init_newton_state(state, n_dim)
    type(newton_state), intent(inout) :: state
    integer, intent(in) :: n_dim
    if (allocated(state%J)) deallocate(state%J, state%x_prev, state%F_prev)
    allocate(state%J(n_dim, n_dim), state%x_prev(n_dim), state%F_prev(n_dim))
    call reset_newton_state(state)
  end subroutine init_newton_state

  subroutine reset_newton_state(state)
    type(newton_state), intent(inout) :: state
    if (.not. allocated(state%J)) return
    state%has_jacobian = .false.
    state%has_prev     = .false.
    state%J            = 0.e0_wp
    state%x_prev       = 0.e0_wp
    state%F_prev       = 0.e0_wp
  end subroutine reset_newton_state

  subroutine to_solver_coords(hc, rep, x)
    real(wp), intent(in)  :: hc, rep
    real(wp), intent(out) :: x(2)
    real(wp) :: rep_clip

    rep_clip = min(max(rep, r_eps), 1.e0_wp - r_eps)
    x(1) = log(max(hc, 1.e-12_wp))
    x(2) = log( (1.e0_wp - rep_clip) / ( rep_clip - r_min_ratio ) ) / rep_map_scale
  end subroutine to_solver_coords

  subroutine from_solver_coords(x, hc, rep)
    real(wp), intent(in)  :: x(2)
    real(wp), intent(out) :: hc, rep
    real(wp) :: exp_arg

    hc = exp(x(1))
    exp_arg = exp( rep_map_scale * x(2) )
    rep = ( exp_arg * r_min_ratio + 1.e0_wp ) / ( exp_arg + 1.e0_wp )
  end subroutine from_solver_coords

  subroutine clamp_step(delta)
    real(wp), intent(inout) :: delta(:)
    ! Only clamp the physical variables, not the path parameter
    integer :: n_clamp
    n_clamp = min(size(delta), 2)
    delta(1:n_clamp) = max(-max_step, min(delta(1:n_clamp), max_step))
  end subroutine clamp_step

  logical function solve_linear(J, rhs, delta)
    real(wp), intent(in)  :: J(:,:), rhs(:)
    real(wp), intent(out) :: delta(:)
    integer :: n, info
    real(wp) :: a, b, c, d, det

    n = size(rhs)
    ! Fast path for the fixed 2x2 systems we solve in the shooting method
    if (n == 2 .and. size(J, 1) == 2 .and. size(J, 2) == 2) then
      a = J(1,1); b = J(1,2)
      c = J(2,1); d = J(2,2)
      det = a * d - b * c
      if (abs(det) < 1.e-24_wp) then
        solve_linear = .false.
        return
      end if
      delta(1) = ( rhs(1) * d - b * rhs(2) ) / det
      delta(2) = ( a * rhs(2) - c * rhs(1) ) / det
      solve_linear = .true.
      return
    end if

    ! Generic fallback keeps allocation on the stack to avoid heap churn
    block
      real(wp) :: J_copy(size(J,1), size(J,2))
      integer :: ipiv(max(1, size(rhs)))
      J_copy = J
      delta = rhs
      call dgesv(n, 1, J_copy, n, ipiv, delta, n, info)
      solve_linear = (info == 0)
    end block
  end function solve_linear

  subroutine broyden_update(state, x, F)
    type(newton_state), intent(inout) :: state
    real(wp), intent(in) :: x(2), F(2)
    real(wp) :: dx(2), dF(2), denom, Jdx(2)

    if (.not. state%has_prev) return
    if (.not. state%has_jacobian) return

    dx = x - state%x_prev
    denom = dot_product(dx, dx)
    if (denom <= 1.e-16_wp) return

    dF  = F - state%F_prev
    Jdx = matmul(state%J, dx)

    state%J = state%J + matmul(reshape(dF - Jdx, (/2,1/)), reshape(dx, (/1,2/))) / denom
  end subroutine broyden_update

  subroutine commit_state(state, x, F)
    type(newton_state), intent(inout) :: state
    real(wp), intent(in) :: x(:), F(:)

    state%x_prev = x
    state%F_prev = F
    state%has_prev = .true.
  end subroutine commit_state

  subroutine line_search(x_current, F_current, delta_x, evaluate_func, final_delta, J_est, success)
    ! Armijo backtracking on phi = 0.5*|F|^2; uses Jacobian estimate when supplied for slope
    procedure(evaluation_function) :: evaluate_func
    real(wp), intent(in)    :: x_current(2), F_current(2), delta_x(2)
    real(wp), intent(out)   :: final_delta(2)
    real(wp), intent(in), optional :: J_est(2,2)
    logical, intent(out), optional :: success
    
    integer, parameter :: max_iter = 10
    real(wp), parameter :: tau = 0.5e0_wp
    real(wp), parameter :: c1 = 1.e-4_wp
    real(wp)             :: alpha
    real(wp)             :: x_trial(2), F_trial(2), hc_trial, rep_trial
    real(wp)             :: rho0_tmp, ee_tmp, er_tmp
    real(wp)             :: phi_old, phi_new, slope0
    integer             :: i
    logical             :: ok

    alpha = 1.e0_wp
    final_delta = delta_x
    phi_old = 0.5e0_wp * dot_product(F_current, F_current)
    ok = .false.

    if (present(J_est)) then
      slope0 = dot_product(F_current, matmul(J_est, delta_x))
    else
      slope0 = -phi_old * 2.e0_wp ! Fallback to a steep descent assumption
    end if
    if (slope0 > -1.e-12_wp) slope0 = -phi_old * 2.e0_wp

    do i = 1, max_iter
      x_trial = x_current + alpha * delta_x
      call from_solver_coords(x_trial, hc_trial, rep_trial)
      call evaluate_func(hc_trial, rep_trial, F_trial, rho0_tmp, ee_tmp, er_tmp)

      phi_new = 0.5e0_wp * dot_product(F_trial, F_trial)
      if (phi_new <= phi_old + c1 * alpha * slope0) then
        ok = .true.
        exit
      end if
      alpha = alpha * tau
    end do
    final_delta = alpha * delta_x
    if (present(success)) success = ok
  end subroutine line_search
end module shoot_solver_mod

module shoot_newton_helpers
  use precision_mod, only: wp
  use para_mod, only: r_ratio, h_center, Mass, MSUN, M_goal, Mass_0, Mb_goal, &
                      J_goal, ang_mom, chi, chi_goal, Omega_c, omc_goal, &
                      Omega_K, C, kappa, Omega_e, FIX1, FIX2
  use shoot_solver_mod
  use rotation_uniform,  only: rotation_solver
  implicit none
contains
  subroutine evaluate_solution(hc, rep, F, rho0, ee, er)
    real(wp), intent(in)  :: hc, rep
    real(wp), intent(out) :: F(2), rho0, ee, er
    real(wp) :: n0_at_h, e_at_h
    external :: n0_at_h, e_at_h
    external :: mass_radius
    real(wp) :: deviA, deviB

    r_ratio  = rep
    h_center = hc

    call rotation_solver
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    select case (trim(FIX1))
    case ('M_goal')
      deviA = Mass/MSUN/M_goal - 1.e0_wp
    case ('Mb_goal')
      deviA = Mass_0/MSUN/Mb_goal - 1.e0_wp
    case default
      stop "evaluate_solution: unknown FIX1"
    end select

    select case (trim(FIX2))
    case ('J_goal')
      deviB = J_goal / ang_mom - 1.e0_wp
    case ('chi_goal')
      deviB = chi / chi_goal - 1.e0_wp
    case ('omc_goal')
      deviB = Omega_c / omc_goal - 1.e0_wp
    case default
      deviB = (Omega_K / (C/sqrt(kappa))) / Omega_e - 1.e0_wp
    end select

    F(1) = deviA
    F(2) = deviB
    er   = abs(deviA) + abs(deviB)
  end subroutine evaluate_solution

  subroutine build_jacobian(state, x, F, hc, rep, rho0, ee, er, reuse_base)
    type(newton_state), intent(inout) :: state
    real(wp), intent(in)    :: x(2)
    real(wp), intent(inout) :: F(2)
    real(wp), intent(inout) :: hc, rep
    real(wp), intent(inout) :: rho0, ee, er
    logical, intent(in), optional :: reuse_base
    real(wp) :: delta(2), xp(2), Fp(2), rho_tmp, ee_tmp, er_tmp
    real(wp) :: hc_p, rep_p
    integer :: i

    delta = max(abs(x), 1.e0_wp) * epsilon(x(1))**(1.e0_wp/3.e0_wp)
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
    ! Caller can skip recomputing the base residual if it was just evaluated
    if (.not. (present(reuse_base) .and. reuse_base)) then
      call evaluate_solution(hc, rep, F, rho0, ee, er)
    end if
  end subroutine build_jacobian

end module shoot_newton_helpers


