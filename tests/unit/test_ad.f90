program test_ad
  use precision_mod, only: wp
  use test_utils
  use ad_mod
  implicit none

  type(dual) :: x, y, z
  real(wp), parameter :: tol = 1.0e-14_8

  ! --- Addition ---
  x%val = 3.0_wp; x%der = 1.0_wp
  y%val = 4.0_wp; y%der = 2.0_wp
  z = x + y
  call assert_near("add_dd val", 7.0_wp, real(z%val, 8), tol)
  call assert_near("add_dd der", 3.0_wp, real(z%der, 8), tol)

  z = x + 10.0d0
  call assert_near("add_dr val", 13.0_wp, real(z%val, 8), tol)
  call assert_near("add_dr der",  1.0_wp, real(z%der, 8), tol)

  ! --- Subtraction ---
  z = x - y
  call assert_near("sub_dd val", -1.0_wp, real(z%val, 8), tol)
  call assert_near("sub_dd der", -1.0_wp, real(z%der, 8), tol)

  z = -x
  call assert_near("neg_d val", -3.0_wp, real(z%val, 8), tol)
  call assert_near("neg_d der", -1.0_wp, real(z%der, 8), tol)

  ! --- Multiplication ---
  z = x * y
  call assert_near("mul_dd val", 12.0_wp, real(z%val, 8), tol)
  ! d/d? (x*y) = der_x*y + x*der_y = 1*4 + 3*2 = 10
  call assert_near("mul_dd der", 10.0d0, real(z%der, 8), tol)

  z = 5.0_wp * x
  call assert_near("mul_rd val", 15.0_wp, real(z%val, 8), tol)
  call assert_near("mul_rd der",  5.0_wp, real(z%der, 8), tol)

  ! --- Division ---
  ! x=3, der=1; y=4, der=2: val=3/4, der=(1 - (3/4)*2)/4 = (1 - 1.5)/4 = -0.125
  z = x / y
  call assert_near("div_dd val",  0.75d0,   real(z%val, 8), tol)
  call assert_near("div_dd der", -0.125d0,  real(z%der, 8), tol)

  z = x / 2.0_wp
  call assert_near("div_dr val", 1.5d0, real(z%val, 8), tol)
  call assert_near("div_dr der", 0.5_wp, real(z%der, 8), tol)

  ! --- exp ---
  ! exp(x): val=e^3, der=e^3 * 1
  x%val = 1.0_wp; x%der = 1.0_wp
  z = exp(x)
  call assert_near("exp val", exp(1.0_wp), real(z%val, 8), tol)
  call assert_near("exp der", exp(1.0_wp), real(z%der, 8), tol)

  ! --- log ---
  ! log(x): val=log(1)=0, der=1/1=1
  z = log(x)
  call assert_near("log val", 0.0d0, real(z%val, 8), tol)
  call assert_near("log der", 1.0_wp, real(z%der, 8), tol)

  ! --- Chain rule: log(exp(x)) = x, derivative = 1 ---
  x%val = 2.0_wp; x%der = 1.0_wp
  z = log(exp(x))
  call assert_near("chain log(exp(x)) val", 2.0_wp, real(z%val, 8), tol)
  call assert_near("chain log(exp(x)) der", 1.0_wp, real(z%der, 8), tol)

  call test_summary("ad_mod")
end program test_ad
