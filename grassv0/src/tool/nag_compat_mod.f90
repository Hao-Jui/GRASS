module nag_compat_mod
  implicit none

  ! to include: F06QWF, D02NBF
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

end module nag_compat_mod
