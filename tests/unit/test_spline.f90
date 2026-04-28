program test_spline
  use precision_mod, only: wp
  use test_utils
  use spline_mod
  implicit none

  integer, parameter :: N = 10
  real(wp) :: xp(N), yp(N)
  type(spline_coeff) :: coeff
  integer :: i
  real(wp) :: t, sx, exact
  real(wp), parameter :: tol_cubic  = 1.0e-10_wp
  real(wp), parameter :: tol_linear = 1.0e-12_wp

  ! -----------------------------------------------------------------
  ! Test 1: Cubic reproduction  y = x^3 on 10 equally-spaced points
  ! A not-a-knot spline reproduces cubics exactly.
  ! -----------------------------------------------------------------
  do i = 1, N
    xp(i) = real(i - 1, wp)
    yp(i) = xp(i)**3
  end do
  call build_spline(xp, yp, N, coeff)

  ! Evaluate at interval midpoints: x = i - 0.5 for i = 1 .. N-1
  do i = 1, N - 1
    t  = 0.5_wp
    sx = yp(i) + coeff%b(i)*t + coeff%c(i)*t**2 + coeff%d(i)*t**3
    exact = (xp(i) + 0.5_wp)**3
    call assert_near("cubic x^3 interval " // char(48+i), &
      exact, sx, tol_cubic)
  end do

  ! -----------------------------------------------------------------
  ! Test 2: Linear data  y = 2x  =>  b = 2, c = 0, d = 0
  ! -----------------------------------------------------------------
  do i = 1, N
    xp(i) = real(i, wp)
    yp(i) = 2.0_wp * xp(i)
  end do
  call build_spline(xp, yp, N, coeff)

  call assert_near("linear b(1)",  2.0_wp, coeff%b(1), tol_linear)
  call assert_near("linear c(1)",  0.0_wp, coeff%c(1), tol_linear)
  call assert_near("linear d(1)",  0.0_wp, coeff%d(1), tol_linear)
  call assert_near("linear b(5)",  2.0_wp, coeff%b(5), tol_linear)
  call assert_near("linear c(5)",  0.0_wp, coeff%c(5), tol_linear)
  call assert_near("linear d(5)",  0.0_wp, coeff%d(5), tol_linear)

  ! -----------------------------------------------------------------
  ! Test 3: Knot interpolation — S(xp(i)) = yp(i) for all knots
  ! Use x^3 data again.
  ! -----------------------------------------------------------------
  do i = 1, N
    xp(i) = real(i - 1, wp)
    yp(i) = xp(i)**3
  end do
  call build_spline(xp, yp, N, coeff)

  do i = 1, N - 1
    t  = 0.0_wp
    sx = yp(i) + coeff%b(i)*t + coeff%c(i)*t**2 + coeff%d(i)*t**3
    call assert_near("knot interp " // char(48+i), yp(i), sx, tol_linear)
  end do

  call test_summary("spline_mod")
end program test_spline
