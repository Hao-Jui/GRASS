module eos_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use precision_mod, only: wp
  use toolkit_mod, only: same_abscissa, pt_interp_action, n_order, bary_w
  use para_mod, only: log_e, log_p, log_h, log_n0, num_tab, n_PT, phase_transition
  use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
  implicit none
  private
  real(wp), allocatable, save :: e_tab(:), p_tab(:)
  real(wp), allocatable, save :: deriv_coeffs_work(:,:)
  ! Public elemental functions
  public :: e_at_p, p_at_e, p_at_e_dual, n0_at_e, n0_at_h, e_at_h, p_at_h, h_at_p
  ! Public subroutines
  public :: loadEos, pressure_derivative_n, pe_at_h
contains

  ! Unified log-space interpolation: dispatches to phase-aware or regular stencil
  real(wp) function interp_eos(val_in, t_in, t_out)
    real(wp), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out(:)
    real(wp) :: res_log
    if (phase_transition) then
      call interp_scalar_pt(t_in, t_out, num_tab, log(val_in), res_log)
    else
      call interp_scalar(t_in, t_out, num_tab, log(val_in), res_log)
    end if
    interp_eos = exp(res_log)
  end function interp_eos

  ! Efficient simultaneous interpolation of two outputs from same input
  subroutine interp_eos_pair(val_in, t_in, t_out1, t_out2, val_out1, val_out2, idx_hint)
    real(wp), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out1(:), t_out2(:)
    real(wp), intent(out) :: val_out1, val_out2
    real(wp) :: res_log1, res_log2
    integer, intent(inout), optional :: idx_hint
    integer :: idx_used

    if (present(idx_hint)) then
      if (phase_transition) then
        call interp_pair_pt(t_in, t_out1, t_out2, num_tab, log(val_in), res_log1, res_log2, idx_hint, idx_used)
      else
        call interp_pair(t_in, t_out1, t_out2, num_tab, log(val_in), res_log1, res_log2, idx_hint, idx_used)
      end if
      idx_hint = idx_used
    else
      if (phase_transition) then
        call interp_pair_pt(t_in, t_out1, t_out2, num_tab, log(val_in), res_log1, res_log2)
      else
        call interp_pair(t_in, t_out1, t_out2, num_tab, log(val_in), res_log1, res_log2)
      end if
    end if

    val_out1 = exp(res_log1)
    val_out2 = exp(res_log2)
  end subroutine interp_eos_pair

  ! Dual-number version for automatic differentiation
  type(dual) function interp_eos_dual(val_in, t_in, t_out)
    use ad_mod, only: log, exp
    type(dual), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out(:)
    type(dual) :: res_log
    if (phase_transition) then
      call interp_dual_pt(t_in, t_out, num_tab, log(val_in), res_log)
    else
      call interp_dual_fast(t_in, t_out, num_tab, log(val_in), res_log)
    end if
    interp_eos_dual = exp(res_log)
  end function interp_eos_dual

  subroutine loadEos
    use iso_fortran_env, only: output_unit
    use para_mod, only: p_at_PT, n_PT, eos_file, num_tab, log_e, log_p, log_h, log_n0, enthalpy_min, C, KSCALE
    implicit none
    integer :: i, unit, ios
    real(wp) :: dle, dlp
    real(wp), parameter :: GAMMA_PT = 6.e-2_wp  ! polytropic-index ceiling for a phase transition
    real(wp), parameter :: DLE_MIN  = 1.e-2_wp  ! minimum d(ln e) to count as sizable
    integer :: pt_buf(1024), pt_count

    pt_count = 0
    if (allocated(p_at_PT)) deallocate(p_at_PT)
    if (allocated(log_e)) deallocate(log_e, log_p, log_h, log_n0)
    if (allocated(e_tab)) deallocate(e_tab, p_tab)

    open(newunit=unit, file="./eos/" // trim(adjustl(eos_file)) // ".dat", status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "loadEos: failed to open EOS table"

    read(unit, *, iostat=ios) num_tab
    if (ios /= 0) error stop "loadEos: failed to read table size"

    allocate(log_e(num_tab), log_p(num_tab), log_h(num_tab), log_n0(num_tab))
    read(unit, *, iostat=ios) log_e(1), log_p(1), log_h(1), log_n0(1)
    if (ios /= 0) error stop "loadEos: malformed EOS row"
    enthalpy_min = log(log_h(1))
    if (log_h(1) <= 1._wp) error stop "loadEos: first enthalpy h(1) <= 1, and thus log(log(h)) is undefined"

    do i = 2, num_tab
      read(unit, *, iostat=ios) log_e(i), log_p(i), log_h(i), log_n0(i)
      if (ios /= 0) error stop "loadEos: malformed EOS row"
      dle = log(log_e(i) / log_e(i-1))
      dlp = log(log_p(i) / log_p(i-1))
      if (dle > DLE_MIN .and. abs(dlp / dle) < GAMMA_PT) then
        pt_count = pt_count + 1
        pt_buf(pt_count) = i
      end if
    end do
    close(unit)

    n_PT = pt_count
    if (n_PT > 0) then
      allocate(p_at_PT(n_PT))
      p_at_PT(:) = pt_buf(1:n_PT)
    else
      allocate(p_at_PT(0))
    end if

    log_e = log(log_e * C * C * KSCALE)
    log_p = log(log_p * KSCALE)
    log_h = log(log(log_h))
    log_n0 = log(log_n0)
    allocate(e_tab(num_tab), p_tab(num_tab))
    e_tab = exp(log_e)
    p_tab = exp(log_p)

    write(output_unit, *) " "
    write(output_unit, *) "# EOS: ", eos_file
    if (n_PT > 0) then
      do i = 1, n_PT
        write(output_unit, "(A,i0,A,1X,i0)") " p_at_PT(", i, "):", p_at_PT(i)
      end do
      phase_transition = .true.
    end if
    write(output_unit, *) "EOS data is in with log(h_min) =", enthalpy_min

    if (any(ieee_is_nan(log_e)) .or. any(ieee_is_nan(log_p)) .or. any(ieee_is_nan(log_h))) error stop "loadEos: wrong table"
  end subroutine loadEos


  ! ========== High-level API: all wrap interp_eos or interp_eos_dual
  real(wp) function e_at_p(pp)
    real(wp), intent(in) :: pp
    e_at_p = interp_eos(pp, log_p, log_e)
  end function e_at_p

  real(wp) function p_at_e(ee)
    real(wp), intent(in) :: ee
    p_at_e = interp_eos(ee, log_e, log_p)
  end function p_at_e

  function p_at_e_dual(ee) result(res)
    type(dual), intent(in) :: ee
    type(dual) :: res
    res = interp_eos_dual(ee, log_e, log_p)
  end function p_at_e_dual

  real(wp) function n0_at_e(ee)
    real(wp), intent(in) :: ee
    n0_at_e = interp_eos(ee, log_e, log_n0)
  end function n0_at_e

  real(wp) function n0_at_h(hh)
    real(wp), intent(in) :: hh
    n0_at_h = interp_eos(hh, log_h, log_n0)
  end function n0_at_h

  real(wp) function e_at_h(hh)
    real(wp), intent(in) :: hh
    e_at_h = interp_eos(hh, log_h, log_e)
  end function e_at_h

  real(wp) function p_at_h(hh)
    real(wp), intent(in) :: hh
    p_at_h = interp_eos(hh, log_h, log_p)
  end function p_at_h

  real(wp) function h_at_p(pp)
    real(wp), intent(in) :: pp
    h_at_p = interp_eos(pp, log_p, log_h)
  end function h_at_p

  subroutine pe_at_h(hh, pp, ee, idx_hint)
    real(wp), intent(in) :: hh
    real(wp), intent(out) :: pp, ee
    integer, intent(inout), optional :: idx_hint

    call interp_eos_pair(hh, log_h, log_p, log_e, pp, ee, idx_hint)
  end subroutine pe_at_h

  pure integer function nearest_monotone_index(xp, xb, idx_hint) result(idx)
    real(wp), intent(in) :: xp(:), xb
    integer, intent(in), optional :: idx_hint
    integer :: lo, hi, mid, n, i, steps
    integer, parameter :: max_hint_steps = 8

    n = size(xp)
    if (n <= 1) then
      idx = 1
      return
    end if

    if (xb <= xp(1)) then
      idx = 1
      return
    end if
    if (xb >= xp(n)) then
      idx = n
      return
    end if

    if (present(idx_hint)) then
      i = min(n, max(1, idx_hint))
      if (same_abscissa(xb, xp(i))) then
        idx = i
        return
      end if

      steps = 0
      if (xb > xp(i)) then
        do while (i < n .and. xp(i+1) <= xb .and. steps < max_hint_steps)
          i = i + 1
          steps = steps + 1
        end do
        if (i < n .and. xp(i) <= xb .and. xb <= xp(i+1)) then
          if (abs(xb - xp(i)) <= abs(xp(i+1) - xb)) then
            idx = i
          else
            idx = i + 1
          end if
          return
        end if
      else
        do while (i > 1 .and. xp(i-1) >= xb .and. steps < max_hint_steps)
          i = i - 1
          steps = steps + 1
        end do
        if (i > 1 .and. xp(i-1) <= xb .and. xb <= xp(i)) then
          if (abs(xb - xp(i-1)) <= abs(xp(i) - xb)) then
            idx = i - 1
          else
            idx = i
          end if
          return
        end if
      end if
    end if

    lo = 1
    hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (xp(mid) <= xb) then
        lo = mid
      else
        hi = mid
      end if
    end do

    if (abs(xb - xp(lo)) <= abs(xp(hi) - xb)) then
      idx = lo
    else
      idx = hi
    end if
  end function nearest_monotone_index

  ! ========== Pair interpolation (two outputs, one input)
  subroutine interp_pair(xp, yp1, yp2, np, xb, y1, y2, idx_hint, idx_used)
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp1(np), yp2(np), xb
    real(wp), intent(out) :: y1, y2
    integer, intent(in), optional :: idx_hint
    integer, intent(out), optional :: idx_used
    integer :: n_nearest_pt, ir, ii
    real(wp) :: dx, wi, den, num1, num2

    n_nearest_pt = nearest_monotone_index(xp, xb, idx_hint)
    if (present(idx_used)) idx_used = n_nearest_pt
    ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))

    call barycentric_pair_stencil(xp, yp1, yp2, ir, xb, y1, y2)
  end subroutine interp_pair

  subroutine interp_pair_pt(xp, yp1, yp2, np, xb, y1, y2, idx_hint, idx_used)
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp1(np), yp2(np), xb
    real(wp), intent(out) :: y1, y2
    integer, intent(in), optional :: idx_hint
    integer, intent(out), optional :: idx_used
    integer :: n_nearest_pt, action, il, ir

    n_nearest_pt = nearest_monotone_index(xp, xb, idx_hint)
    if (present(idx_used)) idx_used = n_nearest_pt
    if (same_abscissa(xb, xp(n_nearest_pt))) then
      y1 = yp1(n_nearest_pt); y2 = yp2(n_nearest_pt); return
    end if

    call pt_interp_action(xp, np, xb, n_nearest_pt, action, il, ir)
    select case (action)
    case (1)
      call linear_segment_pair(xp, yp1, yp2, il, ir, xb, y1, y2)
    case (2)
      y1 = yp1(il); y2 = yp2(il)
    case default
      call barycentric_pair_stencil(xp, yp1, yp2, il, xb, y1, y2)
    end select
  end subroutine interp_pair_pt

  pure subroutine barycentric_pair_stencil(xp, yp1, yp2, i_center, xb, y1, y2)
    real(wp), intent(in) :: xp(:), yp1(:), yp2(:), xb
    integer, intent(in) :: i_center
    real(wp), intent(out) :: y1, y2
    integer :: ii
    real(wp) :: dx, wi, den, num1, num2

    num1 = 0.0_wp; num2 = 0.0_wp; den = 0.0_wp
    do ii = -n_order, n_order
      dx = xb - xp(i_center + ii)
      if (abs(dx) < epsilon(dx)) then
        y1 = yp1(i_center + ii); y2 = yp2(i_center + ii); return
      end if
      wi = bary_w(ii) / dx
      num1 = num1 + wi * yp1(i_center + ii)
      num2 = num2 + wi * yp2(i_center + ii)
      den = den + wi
    end do
    y1 = num1 / den; y2 = num2 / den
  end subroutine barycentric_pair_stencil

  pure subroutine linear_segment_pair(xp, yp1, yp2, i_left, i_right, xb, y1, y2)
    real(wp), intent(in) :: xp(:), yp1(:), yp2(:), xb
    integer, intent(in) :: i_left, i_right
    real(wp), intent(out) :: y1, y2
    real(wp) :: slope1, slope2

    if (i_left >= i_right) then
      y1 = yp1(i_left); y2 = yp2(i_left)
    elseif (same_abscissa(xb, xp(i_left))) then
      y1 = yp1(i_left); y2 = yp2(i_left)
    elseif (same_abscissa(xb, xp(i_right))) then
      y1 = yp1(i_right); y2 = yp2(i_right)
    elseif (same_abscissa(xp(i_left), xp(i_right))) then
      y1 = yp1(i_left); y2 = yp2(i_left)
    else
      slope1 = (yp1(i_right) - yp1(i_left)) / (xp(i_right) - xp(i_left))
      slope2 = (yp2(i_right) - yp2(i_left)) / (xp(i_right) - xp(i_left))
      y1 = yp1(i_left) + (xb - xp(i_left)) * slope1
      y2 = yp2(i_left) + (xb - xp(i_left)) * slope2
    end if
  end subroutine linear_segment_pair

  ! ========== Scalar interpolation (one output)
  subroutine interp_scalar(xp, yp, np, xb, y, idx_hint, idx_used)
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp(np), xb
    real(wp), intent(out) :: y
    integer, intent(in), optional :: idx_hint
    integer, intent(out), optional :: idx_used
    integer :: n_nearest_pt, ir

    n_nearest_pt = nearest_monotone_index(xp, xb, idx_hint)
    if (present(idx_used)) idx_used = n_nearest_pt
    ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))
    call barycentric_scalar_stencil(xp, yp, ir, xb, y)
  end subroutine interp_scalar

  subroutine interp_scalar_pt(xp, yp, np, xb, y, idx_hint, idx_used)
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp(np), xb
    real(wp), intent(out) :: y
    integer, intent(in), optional :: idx_hint
    integer, intent(out), optional :: idx_used
    integer :: n_nearest_pt, action, il, ir

    n_nearest_pt = nearest_monotone_index(xp, xb, idx_hint)
    if (present(idx_used)) idx_used = n_nearest_pt
    if (same_abscissa(xb, xp(n_nearest_pt))) then
      y = yp(n_nearest_pt); return
    end if

    call pt_interp_action(xp, np, xb, n_nearest_pt, action, il, ir)
    select case (action)
    case (1)
      call linear_segment_scalar(xp, yp, il, ir, xb, y)
    case (2)
      y = yp(il)
    case default
      call barycentric_scalar_stencil(xp, yp, il, xb, y)
    end select
  end subroutine interp_scalar_pt

  pure subroutine barycentric_scalar_stencil(xp, yp, i_center, xb, y)
    real(wp), intent(in) :: xp(:), yp(:), xb
    integer, intent(in) :: i_center
    real(wp), intent(out) :: y
    integer :: ii
    real(wp) :: dx, wi, den, num

    num = 0.0_wp; den = 0.0_wp
    do ii = -n_order, n_order
      dx = xb - xp(i_center + ii)
      if (abs(dx) < epsilon(dx)) then
        y = yp(i_center + ii); return
      end if
      wi = bary_w(ii) / dx
      num = num + wi * yp(i_center + ii)
      den = den + wi
    end do
    y = num / den
  end subroutine barycentric_scalar_stencil

  pure subroutine linear_segment_scalar(xp, yp, i_left, i_right, xb, y)
    real(wp), intent(in) :: xp(:), yp(:), xb
    integer, intent(in) :: i_left, i_right
    real(wp), intent(out) :: y

    if (i_left >= i_right .or. same_abscissa(xb, xp(i_left))) then
      y = yp(i_left)
    elseif (same_abscissa(xb, xp(i_right))) then
      y = yp(i_right)
    elseif (same_abscissa(xp(i_left), xp(i_right))) then
      y = yp(i_left)
    else
      y = yp(i_left) + (xb - xp(i_left)) * (yp(i_right) - yp(i_left)) / (xp(i_right) - xp(i_left))
    end if
  end subroutine linear_segment_scalar

  ! ========== Dual interpolation (one output)
  subroutine interp_dual_fast(xp, yp, np, xb, yb)
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp(np)
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb
    integer :: n_nearest_pt, ir

    n_nearest_pt = nearest_monotone_index(xp, xb%val)
    ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))
    call barycentric_dual_stencil(xp, yp, ir, xb, yb)
  end subroutine interp_dual_fast

  subroutine interp_dual_pt(xp, yp, np, xb, yb)
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp(np)
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb
    integer :: n_nearest_pt, action, il, ir

    n_nearest_pt = nearest_monotone_index(xp, xb%val)
    if (same_abscissa(xb%val, xp(n_nearest_pt))) then
      yb = dual_const(yp(n_nearest_pt))
      return
    end if

    call pt_interp_action(xp, np, xb%val, n_nearest_pt, action, il, ir)
    select case (action)
    case (1)
      call linear_segment_dual(xp, yp, il, ir, xb, yb)
    case (2)
      yb = dual_const(yp(il))
    case default
      call barycentric_dual_stencil(xp, yp, il, xb, yb)
    end select
  end subroutine interp_dual_pt

  pure subroutine barycentric_dual_stencil(xp, yp, i_center, xb, yb)
    real(wp), intent(in) :: xp(:), yp(:)
    integer, intent(in) :: i_center
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb
    integer :: ii
    type(dual) :: dx, wi, den, num

    num = dual_const(0.0_wp)
    den = dual_const(0.0_wp)
    do ii = -n_order, n_order
      dx = xb - dual_const(xp(i_center + ii))
      if (abs(dx%val) < epsilon(dx%val)) then
        yb = dual_const(yp(i_center + ii))
        return
      end if
      wi = dual_const(bary_w(ii)) / dx
      num = num + wi * dual_const(yp(i_center + ii))
      den = den + wi
    end do
    yb = num / den
  end subroutine barycentric_dual_stencil

  pure subroutine linear_segment_dual(xp, yp, i_left, i_right, xb, yb)
    real(wp), intent(in) :: xp(:), yp(:)
    integer, intent(in) :: i_left, i_right
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb

    if (i_left >= i_right .or. same_abscissa(xb%val, xp(i_left))) then
      yb = dual_const(yp(i_left))
    elseif (same_abscissa(xb%val, xp(i_right))) then
      yb = dual_const(yp(i_right))
    elseif (same_abscissa(xp(i_left), xp(i_right))) then
      yb = dual_const(yp(i_left))
    else
      yb = dual_const(yp(i_left)) + (xb - dual_const(xp(i_left))) &
         * dual_const((yp(i_right) - yp(i_left)) / (xp(i_right) - xp(i_left)))
    end if
  end subroutine linear_segment_dual

  ! ========== Derivative computation (separate subsystem) ==========
  ! Evaluates d^n p / d e^n using Fornberg finite-difference weights
  ! on a clipped energy range. Error flags: 1=bad order, 2=no EOS, 3=bad energy,
  ! 4=small stencil, 5=out-of-bounds (clamped)
  subroutine pressure_derivative_n(ee, n, derivative, status, idx_hint)
    implicit none
    real(wp), intent(in) :: ee
    integer, intent(in) :: n
    real(wp), intent(out) :: derivative
    integer, intent(out), optional :: status
    integer, intent(inout), optional :: idx_hint

    integer :: info, idx, half_width, left, right, n_points
    real(wp) :: x0, x0_log, min_e, max_e
    real(wp) :: ee_clamped
    logical :: fatal_error
    real(wp), parameter :: CLAMP_TOL = 1.e-12_wp

    info = 0
    fatal_error = .false.
    derivative = 0._wp

    if (n < 0) then
      info = 1
      fatal_error = .true.
    else if (.not. allocated(log_e) .or. .not. allocated(log_p) .or. &
             .not. allocated(e_tab) .or. .not. allocated(p_tab) .or. num_tab <= 0) then
      info = 2
      fatal_error = .true.
    else if (n >= num_tab) then
      info = 4
      fatal_error = .true.
    else if (ee <= 0._wp) then
      info = 3
      fatal_error = .true.
    end if

    if (fatal_error) then
      if (present(status)) status = info
      derivative = 0._wp
      return
    end if

    min_e = e_tab(1)
    max_e = e_tab(num_tab)

    ee_clamped = max(min_e, min(max_e, ee))
    if (abs(ee_clamped - ee) > clamp_tol * max(1._wp, abs(ee_clamped))) info = 5

    x0 = ee_clamped
    x0_log = log(x0)
    idx = nearest_monotone_index(log_e, x0_log, idx_hint)
    if (present(idx_hint)) idx_hint = idx

    half_width = max(n, 4)
    half_width = min(half_width, num_tab - 1)

    left = idx - half_width
    right = idx + half_width

    if (left < 1) then
      right = min(num_tab, right + (1 - left))
      left = 1
    end if
    if (right > num_tab) then
      left = max(1, left - (right - num_tab))
      right = num_tab
    end if

    n_points = right - left + 1

    if (n_points <= n) then
      info = 4
      fatal_error = .true.
      if (present(status)) status = info
      derivative = 0._wp
      return
    end if
    
    if (.not. allocated(deriv_coeffs_work)) then
      allocate(deriv_coeffs_work(n_points, n+1))
    else if (size(deriv_coeffs_work, 1) < n_points .or. size(deriv_coeffs_work, 2) < n+1) then
      deallocate(deriv_coeffs_work)
      allocate(deriv_coeffs_work(n_points, n+1))
    end if

    call fornberg_weights(x0, e_tab(left:right), n_points, n, deriv_coeffs_work(1:n_points, 1:n+1))
    derivative = dot_product(deriv_coeffs_work(1:n_points, n+1), p_tab(left:right))

    if (abs(ee_clamped - ee) > clamp_tol * max(1._wp, abs(ee_clamped))) &
      info = max(info, 5)

    if (present(status)) status = info
    if (fatal_error) derivative = 0._wp

  contains
    ! helper that constructs finite-difference weights for arbitrary 
    ! derivative order, allowing the outer routine to reuse the existing 
    ! EOS samples without new interpolation logic.
    subroutine fornberg_weights(x0_local, x, m, n_deriv, c)
      implicit none
      integer, intent(in) :: m, n_deriv
      real(wp), intent(in) :: x0_local
      real(wp), intent(in) :: x(m)
      real(wp), intent(out) :: c(m, n_deriv + 1)

      integer :: i, j, k, max_k
      real(wp) :: c1, c2, c3, c4, c5

      c(:, :) = 0._wp
      c(1,1) = 1._wp
      c1 = 1._wp
      c4 = x(1) - x0_local

      do i = 2, m
        c2 = 1._wp
        c5 = c4
        c4 = x(i) - x0_local
        max_k = min(i - 1, n_deriv)
        do j = 1, i - 1
          c3 = x(i) - x(j)
          c2 = c2 * c3
          if (j == i - 1) then
            do k = max_k, 1, -1
              c(i, k+1) = c1 * (k * c(i-1, k) - c5 * c(i-1, k+1)) / c2
            end do
            c(i,1) = -c1 * c5 * c(i-1,1) / c2
          end if
          do k = max_k, 1, -1
            c(j, k+1) = (c4 * c(j, k+1) - k * c(j, k)) / c3
          end do
          c(j,1) = c4 * c(j,1) / c3
        end do
        c1 = c2
      end do
    end subroutine fornberg_weights
  end subroutine pressure_derivative_n

end module eos_mod
