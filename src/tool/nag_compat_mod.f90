module nag_compat_mod
  implicit none

  ! to include: F06QWF, D02NBF
  abstract interface
    subroutine nag_ode_rhs(t, y, dy)
      real(8), intent(in) :: t
      real(8), intent(in) :: y(:)
      real(8), intent(out) :: dy(:)
    end subroutine nag_ode_rhs
  end interface
contains

  subroutine d01gaf(x, y, n, ans, er, ifail)
    implicit none
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n), y(n)
    real(8), intent(out) :: ans, er
    integer, intent(out) :: ifail

    integer :: i, last
    real(8) :: c, d1, d2, d3, h1, h2, h3, h4
    real(8) :: r1, r2, r3, r4, s
    real(8) :: first_step

    ifail = 0
    ans   = 0.d0
    er    = 0.d0

    if (n < 4) then
      ifail = 1
      return
    end if

    if (size(x) /= size(y)) then
      ifail = 2
      return
    end if

    first_step = x(2) - x(1)
    if (first_step == 0.d0) then
      ifail = 3
      return
    end if

    do i = 3, n
      h3 = x(i) - x(i-1)
      if (h3 * first_step <= 0.d0) then
        ifail = 3
        return
      end if
    end do

    ! piecewise cubic / higher-order Newton–Cotes adaptive-like interpolatory method
    h2 = first_step
    d3 = (y(2) - y(1)) / h2
    h3 = x(3) - x(2)
    d1 = (y(3) - y(2)) / h3
    h1 = h2 + h3
    d2 = (d1 - d3) / h1
    h4 = x(4) - x(3)
    r1 = (y(4) - y(3)) / h4
    r2 = (r1 - d1) / (h4 + h3)
    h1 = h1 + h4
    r3 = (r2 - d2) / h1

    ans = h2 * ( y(1) + h2 * ( d3/2.d0 - h2 * ( d2/6.d0 - (h2 + 2.d0*h3) * r3 / 12.d0 ) ) )
    s   = -h2**3 * ( h2*(3.d0*h2 + 5.d0*h4) + 10.d0*h3*h1 ) / 60.d0
    r4  = 0.d0

    last = n - 1
    do i = 3, last
      ans = ans + h3 * ( (y(i) + y(i-1))/2.d0 - h3*h3*(d2 + r2 + (h2 - h4)*r3)/12.d0 )
      c   = h3**3 * ( 2.d0*h3*h3 + 5.d0*( h3*(h4 + h2) + 2.d0*h4*h2 ) ) / 120.d0
      er  = er + (c + s) * r4

      if (i == 3) then
        s = s + 2.d0 * c
      else
        s = c
      end if

      if (i < last) then
        h1 = h2
        h2 = h3
        h3 = h4
        d1 = r1
        d2 = r2
        d3 = r3
        h4 = x(i+2) - x(i+1)
        r1 = (y(i+2) - y(i+1)) / h4
        r4 = h4 + h3
        r2 = (r1 - d1) / r4
        r4 = r4 + h2
        r3 = (r2 - d2) / r4
        r4 = r4 + h1
        r4 = (r3 - d3) / r4
      else
        ans = ans + h4 * ( y(n) - h4 * ( r1/2.d0 + h4 * ( r2/6.d0 + (2.d0*h3 + h4) * r3 / 12.d0 ) ) )
        er  = er - h4**3 * r4 * ( h4*(3.d0*h4 + 5.d0*h2) + 10.d0*h3*(h2 + h3 + h4) ) / 60.d0 + s * r4
        ans = ans + er
      end if
    end do
  end subroutine d01gaf

  subroutine e02baf(n, x, y, a, b, c, d, info)
    !! Fit a natural cubic spline passing through (x, y).
    !! The spline on interval [x_i, x_{i+1}] is
    !!   S_i(t) = a_i + b_i * t + c_i * t**2 + d_i * t**3,  t = (x - x_i)
    !! Inputs:
    !!   n   - number of data points (n >= 2)
    !!   x   - strictly increasing abscissae (length n)
    !!   y   - ordinates (length n)
    !! Outputs:
    !!   a,b,c,d - spline coefficients for each of the n-1 intervals
    !!   info    - 0 on success, non-zero otherwise
    integer, intent(in) :: n
    real(8), intent(in) :: x(n), y(n)
    real(8), intent(out) :: a(n-1), b(n-1), c(n-1), d(n-1)
    integer, intent(out) :: info

    real(8), allocatable :: h(:), alpha(:), l(:), mu(:), z(:), c_full(:)
    integer :: i

    info = 0

    if (n < 2) then
      info = 1
      return
    end if

    do i = 2, n
      if (x(i) <= x(i-1)) then
        info = 2
        return
      end if
    end do

    allocate(h(n-1), alpha(max(1,n-2)))
    do i = 1, n-1
      h(i) = x(i+1) - x(i)
      if (h(i) <= 0.d0) then
        info = 2
        deallocate(h, alpha)
        return
      end if
    end do

    if (n > 2) then
      do i = 2, n-1
        alpha(i-1) = (3.d0 / h(i)) * (y(i+1) - y(i)) - (3.d0 / h(i-1)) * (y(i) - y(i-1))
      end do
    end if

    allocate(l(n), mu(n), z(n), c_full(n))
    l(1) = 1.d0
    mu(1) = 0.d0
    z(1) = 0.d0

    if (n > 2) then
      do i = 2, n - 1
        l(i) = 2.d0 * (x(i+1) - x(i-1)) - h(i-1) * mu(i-1)
        if (l(i) == 0.d0) then
          info = 3
          deallocate(h, alpha, l, mu, z, c_full)
          return
        end if
        mu(i) = h(i) / l(i)
        z(i) = (alpha(i-1) - h(i-1) * z(i-1)) / l(i)
      end do
    end if

    l(n) = 1.d0
    z(n) = 0.d0
    c_full(n) = 0.d0

    do i = n - 1, 1, -1
      c_full(i) = z(i) - mu(i) * c_full(i+1)
      b(i) = (y(i+1) - y(i)) / h(i) - h(i) * (c_full(i+1) + 2.d0 * c_full(i)) / 3.d0
      d(i) = (c_full(i+1) - c_full(i)) / (3.d0 * h(i))
      a(i) = y(i)
      c(i) = c_full(i)
    end do

    deallocate(h, alpha, l, mu, z, c_full)
  end subroutine e02baf

  subroutine e02bbf(n, x, a, b, c, d, aint, bint, result, info)
    !! Integrate the cubic spline defined by (x,a,b,c,d) from aint to bint.
    !! Inputs:
    !!   n      - number of data points defining the spline (n >= 2)
    !!   x      - knot positions (length n, strictly increasing)
    !!   a,b,c,d- spline coefficients as returned by e02baf (length n-1 each)
    !!   aint   - lower limit of integration
    !!   bint   - upper limit of integration
    !! Outputs:
    !!   result - definite integral of the spline from aint to bint
    !!   info   - 0 on success, non-zero otherwise
    integer, intent(in) :: n
    real(8), intent(in) :: x(n), a(n-1), b(n-1), c(n-1), d(n-1)
    real(8), intent(in) :: aint, bint
    real(8), intent(out) :: result
    integer, intent(out) :: info

    integer :: i, start_idx, end_idx
    real(8) :: lower, upper, segment_a, segment_b
    real(8) :: dl, du

    info = 0
    result = 0.d0

    if (n < 2) then
      info = 1
      return
    end if

    do i = 2, n
      if (x(i) <= x(i-1)) then
        info = 2
        return
      end if
    end do

    if (aint == bint) return

    lower = min(aint, bint)
    upper = max(aint, bint)

    if (lower < x(1) .or. upper > x(n)) then
      info = 3
      return
    end if

    start_idx = 1
    do while (start_idx < n .and. x(start_idx+1) <= lower)
      start_idx = start_idx + 1
    end do

    end_idx = start_idx
    do while (end_idx < n .and. x(end_idx+1) < upper)
      end_idx = end_idx + 1
    end do

    do i = start_idx, end_idx
      segment_a = merge(lower, x(i), i == start_idx)
      segment_b = merge(upper, x(i+1), i == end_idx)

      dl = segment_a - x(i)
      du = segment_b - x(i)

      result = result + a(i) * (du - dl)                                     &
                      + b(i) * (du**2 - dl**2) / 2.d0                        &
                      + c(i) * (du**3 - dl**3) / 3.d0                        &
                      + d(i) * (du**4 - dl**4) / 4.d0
    end do

    if (aint > bint) result = -result
  end subroutine e02bbf

  subroutine d02pcf(f, neqn, y, yp, t, tout, relerr, abserr, flag, step_count, out, max_step)
  ! rkh45
    implicit none
    procedure(nag_ode_rhs) :: f
    integer, intent(in) :: neqn
    real(8), intent(inout) :: y(neqn)
    real(8), intent(inout) :: yp(neqn)
    real(8), intent(inout) :: t
    real(8), intent(in) :: tout
    real(8), intent(inout) :: relerr
    real(8), intent(inout) :: abserr
    integer, intent(inout) :: flag
    integer, intent(out) :: step_count
    logical, intent(in) :: out
    real(8), intent(in), optional :: max_step

    real(8) :: distance, direction, h, hmin, hmax
    real(8) :: err, fac, scale, tol_small
    integer :: max_steps, max_evals, nfe, i
    logical :: eval_limit
    real(8) :: y4(neqn), y5(neqn)
    real(8) :: k1(neqn), k2(neqn), k3(neqn), k4(neqn), k5(neqn), k6(neqn)


    max_steps = 200000
    max_evals = 6 * max_steps
    step_count = 0
    nfe = 0
    eval_limit = .false.

    if (neqn <= 0) then
      flag = 8
      return
    end if

    relerr = max(relerr, 1.d-12)
    abserr = max(abserr, 1.d-18)

    if (flag == 0 .or. abs(flag) > 2) then
      flag = 8
      return
    end if

    distance = tout - t
    if (distance == 0.d0) then
      call f(t, y, yp)
      flag = 2
      return
    end if

    direction = sign(1.d0, distance)
    call f(t, y, yp)
    nfe = nfe + 1
    if (nfe > max_evals) then
      flag = 4
      return
    end if

    hmin = 1.d-12
    if (present(max_step)) then
      hmax = merge(min(abs(distance), max_step), abs(distance), max_step > 0.d0)
    else
      hmax = abs(distance)
    end if
    h = direction * max(1.d-6, min(abs(distance) / 10.d0, hmax))
    tol_small = 10.d0 * epsilon(t)

    do
      distance = tout - t
      if (direction * distance <= 0.d0) exit
      if (step_count >= max_steps) then
        flag = 4
        return
      end if

      if (abs(h) > abs(distance)) h = direction * abs(distance)
      if (abs(h) > hmax) h = direction * hmax
      if (abs(h) < hmin) then
        flag = 6
        return
      end if
      call rkf45_step(f, neqn, t, y, h, y4, y5, k1, k2, k3, k4, k5, k6, eval_limit, max_evals, nfe)
      if (eval_limit) then
        flag = 4
        return
      end if

      err = 0.d0
      do i = 1, neqn
        scale = abserr + relerr * max(abs(y(i)), abs(y5(i)))
        if (scale > tol_small) then
          err = max(err, abs(y5(i) - y4(i)) / scale)
        end if
      end do
      err = err / sqrt( dble(neqn) )

      if (err <= 1.d0) then
        t = t + h
        y = y5
        step_count = step_count + 1
        call f(t, y, yp)
        nfe = nfe + 1
        if (nfe > max_evals) then
          flag = 4
          return
        end if
        if (out) then
          write(*,"(99es27.17e3)") t, y, yp
        end if
        if (abs(distance) <= 1.d-15) exit
        fac = min(5.d0, 0.9d0 * err**(-0.2d0))
        h = direction * min( abs(h) * fac, hmax )
      else
        fac = max(0.1d0, 0.9d0 * err**(-0.25d0))
        h = direction * min( abs(h) * fac, hmax )
        if (abs(h) < hmin) then
          flag = 6
          return
        end if
      end if
    end do

    flag = 2
  contains
    subroutine rkf45_step(f, neqn, t, y, h, y4loc, y5loc, k1, k2, k3, k4, k5, k6, limit_flag, max_eval, nfe_loc)
      implicit none
      procedure(nag_ode_rhs) :: f
      integer, intent(in) :: neqn, max_eval
      real(8), intent(in) :: t, h
      real(8), intent(in) :: y(neqn)
      real(8), intent(out) :: y4loc(neqn), y5loc(neqn)
      real(8), intent(out) :: k1(neqn), k2(neqn), k3(neqn), k4(neqn), k5(neqn), k6(neqn)
      logical, intent(out) :: limit_flag
      integer, intent(inout) :: nfe_loc

      real(8) :: ywork(neqn)
      real(8), parameter :: a2 = 0.25d0, a3 = 3.d0/8.d0, a4 = 12.d0/13.d0, a5 = 1.d0, a6 = 0.5d0
      real(8), parameter :: b21 = 0.25d0
      real(8), parameter :: b31 = 3.d0/32.d0, b32 = 9.d0/32.d0
      real(8), parameter :: b41 = 1932.d0/2197.d0, b42 = -7200.d0/2197.d0, b43 = 7296.d0/2197.d0
      real(8), parameter :: b51 = 439.d0/216.d0, b52 = -8.d0, b53 = 3680.d0/513.d0, b54 = -845.d0/4104.d0
      real(8), parameter :: b61 = -8.d0/27.d0, b62 = 2.d0, b63 = -3544.d0/2565.d0, b64 = 1859.d0/4104.d0, b65 = -11.d0/40.d0
      real(8), parameter :: c1 = 16.d0/135.d0, c3 = 6656.d0/12825.d0, c4 = 28561.d0/56430.d0, c5 = -9.d0/50.d0, c6 = 2.d0/55.d0
      real(8), parameter :: ch1 = 25.d0/216.d0, ch3 = 1408.d0/2565.d0, ch4 = 2197.d0/4104.d0, ch5 = -1.d0/5.d0

      limit_flag = .false.

      call f(t, y, k1)
      nfe_loc = nfe_loc + 1
      if (nfe_loc > max_eval) then
        limit_flag = .true.
        return
      end if

      ywork = y + h * b21 * k1
      call f(t + a2*h, ywork, k2)
      nfe_loc = nfe_loc + 1
      if (nfe_loc > max_eval) then
        limit_flag = .true.
        return
      end if

      ywork = y + h * (b31 * k1 + b32 * k2)
      call f(t + a3*h, ywork, k3)
      nfe_loc = nfe_loc + 1
      if (nfe_loc > max_eval) then
        limit_flag = .true.
        return
      end if

      ywork = y + h * (b41 * k1 + b42 * k2 + b43 * k3)
      call f(t + a4*h, ywork, k4)
      nfe_loc = nfe_loc + 1
      if (nfe_loc > max_eval) then
        limit_flag = .true.
        return
      end if

      ywork = y + h * (b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4)
      call f(t + a5*h, ywork, k5)
      nfe_loc = nfe_loc + 1
      if (nfe_loc > max_eval) then
        limit_flag = .true.
        return
      end if

      ywork = y + h * (b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5)
      call f(t + a6*h, ywork, k6)
      nfe_loc = nfe_loc + 1
      if (nfe_loc > max_eval) then
        limit_flag = .true.
        return
      end if

      y5loc = y + h * (c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5 + c6 * k6)
      y4loc = y + h * (ch1 * k1 + ch3 * k3 + ch4 * k4 + ch5 * k5)
    end subroutine rkf45_step

  end subroutine d02pcf

end module nag_compat_mod
