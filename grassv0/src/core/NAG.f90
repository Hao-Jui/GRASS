!----------------------------------------------------------------------------!
!  d01gaf  -  Numerical integration of tabulated data (Gill & Miller method)
!----------------------------------------------------------------------------!
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

  ! Ensure the abscissae are strictly monotonic.
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

  ! Initial interval (points 1-4)
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
