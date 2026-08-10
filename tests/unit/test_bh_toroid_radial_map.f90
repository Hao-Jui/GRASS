program test_bh_toroid_radial_map
  use precision_mod, only: wp
  use radial_map_mod
  use validation_mod, only: VALID_OK, VALID_BAD_GRID_SIZE, VALID_BAD_HORIZON, VALID_BAD_RADIAL_ORDER
  use test_utils
  implicit none

  integer, parameter :: n = 5
  real(wp), parameter :: h0 = 0.2_wp, rin = 0.55_wp, rout = 10.0_wp
  real(wp) :: rhat(n), weights(n), dr

  call assert_status("valid domain", validate_radial_domain(h0, rin), VALID_OK)
  call assert_status("bad horizon domain", validate_radial_domain(0.0_wp, rin), VALID_BAD_HORIZON)
  call assert_status("bad radial order domain", validate_radial_domain(rin, rin), VALID_BAD_RADIAL_ORDER)
  call assert_status("bad outer domain", validate_radial_domain(h0, 1.0_wp), VALID_BAD_RADIAL_ORDER)

  call assert_status("grid builds", build_rhat_trapezoid_grid(n, h0, rhat, weights), VALID_OK)
  dr = (1.0_wp - h0) / real(n - 1, wp)
  call assert_near("left endpoint", h0, rhat(1), 1.e-14_wp)
  call assert_near("right endpoint", 1.0_wp, rhat(n), 1.e-14_wp)
  call assert_true("monotone grid", all(rhat(2:n) > rhat(1:n-1)))
  call assert_true("no point below horizon", all(rhat >= h0))
  call assert_near("left half weight", 0.5_wp * dr, weights(1), 1.e-14_wp)
  call assert_near("interior weight", dr, weights(2), 1.e-14_wp)
  call assert_near("right half weight", 0.5_wp * dr, weights(n), 1.e-14_wp)
  call assert_near("weight sum", 1.0_wp - h0, sum(weights), 1.e-14_wp)
  call assert_true("positive weights", all(weights > 0.0_wp))

  call assert_near("horizon radius", rout * h0, radius_from_rhat(rout, h0), 1.e-14_wp)
  call assert_near("outer radius", rout, radius_from_rhat(rout, 1.0_wp), 1.e-14_wp)
  call assert_near("jacobian", rout, radial_jacobian(rout), 1.e-14_wp)

  call assert_true("classify horizon", classify_rhat(h0, h0, rin) == ZONE_HORIZON)
  call assert_true("classify vacuum gap", classify_rhat(0.4_wp, h0, rin) == ZONE_VACUUM_GAP)
  call assert_true("classify torus boundary", classify_rhat(rin, h0, rin) == ZONE_TORUS)
  call assert_true("classify outside", classify_rhat(1.1_wp, h0, rin) == ZONE_OUTSIDE)

  call assert_status("bad grid size", build_rhat_trapezoid_grid(1, h0, rhat, weights), VALID_BAD_GRID_SIZE)
  call assert_status("bad grid horizon", build_rhat_trapezoid_grid(n, 1.0_wp, rhat, weights), VALID_BAD_HORIZON)

  call test_summary("test_bh_toroid_radial_map")

contains

  subroutine assert_status(label, result, expected)
    use validation_mod, only: validation_result
    character(*), intent(in) :: label
    type(validation_result), intent(in) :: result
    integer, intent(in) :: expected

    call assert_true(label, result%status == expected)
  end subroutine assert_status
end program test_bh_toroid_radial_map
