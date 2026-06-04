program test_bh_toroid_green
  use precision_mod, only: wp
  use bh_toroid_green_mod
  use bh_toroid_validation_mod, only: VALID_OK, VALID_BAD_GREEN_ARGS
  use test_utils
  implicit none

  real(wp), parameter :: tol = 1.e-12_wp
  real(wp), parameter :: h0 = 0.25_wp
  real(wp), parameter :: r_left = 0.5_wp, rp_left = 0.75_wp
  real(wp), parameter :: r_equal = 0.6_wp, rp_equal = 0.6_wp
  real(wp), parameter :: r_right = 0.75_wp, rp_right = 0.5_wp

  ! Oracle table from paper eqs. (3.4)-(3.5):
  ! h0=0.25, r=0.5, r'=0.75, r/r'=2/3:
  !   f1_0 = 0, f1_1 = 1/2, f1_2 = 5/12.
  !   f2_0 = 2/3, f2_1 = 7/9, f2_2 = 31/54.
  ! h0=0.25, r=r'=0.6:
  !   f1_1 = 119/144, f1_2 = 20111/20736.
  !   f2_0 = 35/36, f2_1 = 8015/5184, f2_2 = 1228535/746496.
  call assert_status("valid green args", validate_green_args(0, r_left, rp_left, h0), VALID_OK)
  call assert_status("bad green args", validate_green_args(-1, r_left, rp_left, h0), VALID_BAD_GREEN_ARGS)
  call assert_status("bad horizon relation", validate_green_args(0, h0, rp_left, h0), VALID_BAD_GREEN_ARGS)

  call assert_near("f1 n0 left", 0.0_wp, ne_f1_kernel(0, r_left, rp_left, h0), tol)
  call assert_near("f1 n1 left", 0.5_wp, ne_f1_kernel(1, r_left, rp_left, h0), tol)
  call assert_near("f1 n2 left", 5.0_wp / 12.0_wp, ne_f1_kernel(2, r_left, rp_left, h0), tol)
  call assert_near("f2 n0 left", 2.0_wp / 3.0_wp, ne_f2_kernel(0, r_left, rp_left, h0), tol)
  call assert_near("f2 n1 left", 7.0_wp / 9.0_wp, ne_f2_kernel(1, r_left, rp_left, h0), tol)
  call assert_near("f2 n2 left", 31.0_wp / 54.0_wp, ne_f2_kernel(2, r_left, rp_left, h0), tol)

  call assert_near("f1 n1 equality", 119.0_wp / 144.0_wp, ne_f1_kernel(1, r_equal, rp_equal, h0), tol)
  call assert_near("f1 n2 equality", 20111.0_wp / 20736.0_wp, ne_f1_kernel(2, r_equal, rp_equal, h0), tol)
  call assert_near("f2 n0 equality", 35.0_wp / 36.0_wp, ne_f2_kernel(0, r_equal, rp_equal, h0), tol)
  call assert_near("f2 n1 equality", 8015.0_wp / 5184.0_wp, ne_f2_kernel(1, r_equal, rp_equal, h0), tol)
  call assert_near("f2 n2 equality", 1228535.0_wp / 746496.0_wp, ne_f2_kernel(2, r_equal, rp_equal, h0), tol)

  call assert_near("f1 branch symmetry", ne_f1_kernel(1, r_left, rp_left, h0), &
      ne_f1_kernel(1, r_right, rp_right, h0), tol)
  call assert_near("f2 branch symmetry", ne_f2_kernel(1, r_left, rp_left, h0), &
      ne_f2_kernel(1, r_right, rp_right, h0), tol)

  call assert_near("f1 h0=0 limit", 2.0_wp / 3.0_wp, ne_f1_kernel(1, r_left, rp_left, 0.0_wp), tol)
  call assert_near("f2 h0=0 limit", 8.0_wp / 9.0_wp, ne_f2_kernel(1, r_left, rp_left, 0.0_wp), tol)

  call assert_near("f1 horizon image cancellation", 0.0_wp, ne_f1_kernel(2, h0, rp_left, h0), tol)
  call assert_near("f2 horizon image cancellation", 0.0_wp, ne_f2_kernel(2, h0, rp_left, h0), tol)
  call assert_true("finite near horizon", abs(ne_f2_kernel(2, h0 + 1.e-6_wp, rp_left, h0)) < huge(1.0_wp))

  call assert_near("lambda wrapper", ne_f2_kernel(1, r_left, rp_left, h0), &
      lambda_radial_kernel(1, r_left, rp_left, h0), tol)
  call assert_near("B wrapper", ne_f1_kernel(1, r_left, rp_left, h0), &
      b_radial_kernel(1, r_left, rp_left, h0), tol)
  call assert_near("omega wrapper", ne_f2_kernel(1, r_left, rp_left, h0), &
      omega_radial_kernel(1, r_left, rp_left, h0), tol)

  call test_summary("test_bh_toroid_green")

contains

  subroutine assert_status(label, result, expected)
    use bh_toroid_validation_mod, only: validation_result
    character(*), intent(in) :: label
    type(validation_result), intent(in) :: result
    integer, intent(in) :: expected

    call assert_true(label, result%status == expected)
  end subroutine assert_status
end program test_bh_toroid_green
