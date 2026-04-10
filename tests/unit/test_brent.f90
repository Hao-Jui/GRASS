program test_brent
  use precision_mod, only: wp
  use test_utils
  use brent_mod
  implicit none

  real(wp) :: result_val
  real(wp), parameter :: tol_root = 1.0e-8_wp
  real(wp), parameter :: tol_assert = 1.0e-6_wp

  ! -----------------------------------------------------------------
  ! Test 1: f(x) = x - 2, root at x = 2.
  ! Initial guess = 3.0; scale_up = 1.2, scale_down = 1.1.
  ! After iter 1: a = 3.6, b = 3.0/1.1 ≈ 2.727 => both positive.
  ! After several iters b shrinks below 2 and a grows above 2, bracketing.
  ! -----------------------------------------------------------------
  call brent_core(3.0_wp, 1.2_wp, 1.1_wp, tol_root, result_val, linear_residual)
  call assert_near("root of x-2", 2.0_wp, result_val, tol_assert)
  call assert_true("root of x-2 is near 2", abs(result_val - 2.0_wp) < tol_assert)

  ! -----------------------------------------------------------------
  ! Test 2: f(x) = x^2 - 9, positive root at x = 3.
  ! Initial guess = 5.0; as b shrinks below 3 and a stays above 3 the
  ! root is bracketed.
  ! -----------------------------------------------------------------
  call brent_core(5.0_wp, 1.1_wp, 1.1_wp, tol_root, result_val, quadratic_residual)
  call assert_near("root of x^2-9", 3.0_wp, result_val, tol_assert)
  call assert_true("root of x^2-9 is near 3", abs(result_val - 3.0_wp) < tol_assert)

  call test_summary("brent_mod")

contains

  subroutine linear_residual(x, fx)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: fx
    fx = x - 2.0_wp
  end subroutine linear_residual

  subroutine quadratic_residual(x, fx)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: fx
    fx = x**2 - 9.0_wp
  end subroutine quadratic_residual

end program test_brent
