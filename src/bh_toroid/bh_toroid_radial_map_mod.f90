module bh_toroid_radial_map_mod
  use precision_mod, only: wp
  use bh_toroid_validation_mod, only: validation_result, validation_ok, validation_error, &
      VALID_BAD_GRID_SIZE, VALID_BAD_HORIZON, VALID_BAD_RADIAL_ORDER
  implicit none
  private

  integer, parameter, public :: ZONE_HORIZON = 1
  integer, parameter, public :: ZONE_VACUUM_GAP = 2
  integer, parameter, public :: ZONE_TORUS = 3
  integer, parameter, public :: ZONE_OUTSIDE = 4

  public :: validate_radial_domain, build_rhat_trapezoid_grid, classify_rhat
  public :: radius_from_rhat, radial_jacobian

contains

  pure function validate_radial_domain(h0_hat, rin_hat) result(res)
    real(wp), intent(in) :: h0_hat, rin_hat
    type(validation_result) :: res

    if (h0_hat <= 0.0_wp .or. h0_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_HORIZON, "horizon coordinate must satisfy 0 < h0_hat < 1")
    else if (rin_hat <= h0_hat .or. rin_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_RADIAL_ORDER, "radial domain must satisfy h0_hat < rin_hat < 1")
    else
      res = validation_ok()
    end if
  end function validate_radial_domain

  function build_rhat_trapezoid_grid(n, h0_hat, rhat, weights) result(res)
    integer, intent(in) :: n
    real(wp), intent(in) :: h0_hat
    real(wp), intent(out) :: rhat(:), weights(:)
    type(validation_result) :: res
    integer :: i
    real(wp) :: dr

    if (n < 2 .or. size(rhat) < n .or. size(weights) < n) then
      res = validation_error(VALID_BAD_GRID_SIZE, "rhat grid requires n >= 2 and matching output arrays")
      return
    end if
    if (h0_hat <= 0.0_wp .or. h0_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_HORIZON, "horizon coordinate must satisfy 0 < h0_hat < 1")
      return
    end if

    dr = (1.0_wp - h0_hat) / real(n - 1, wp)
    do i = 1, n
      rhat(i) = h0_hat + real(i - 1, wp) * dr
      weights(i) = dr
    end do
    weights(1) = 0.5_wp * dr
    weights(n) = 0.5_wp * dr
    res = validation_ok()
  end function build_rhat_trapezoid_grid

  pure elemental function classify_rhat(rhat, h0_hat, rin_hat) result(zone)
    real(wp), intent(in) :: rhat, h0_hat, rin_hat
    integer :: zone
    real(wp) :: tol

    tol = 64.0_wp * epsilon(1.0_wp)
    if (rhat <= h0_hat + tol) then
      zone = ZONE_HORIZON
    else if (rhat < rin_hat) then
      zone = ZONE_VACUUM_GAP
    else if (rhat <= 1.0_wp + tol) then
      zone = ZONE_TORUS
    else
      zone = ZONE_OUTSIDE
    end if
  end function classify_rhat

  pure elemental function radius_from_rhat(rout, rhat) result(radius)
    real(wp), intent(in) :: rout, rhat
    real(wp) :: radius

    radius = rout * rhat
  end function radius_from_rhat

  pure function radial_jacobian(rout) result(jacobian)
    real(wp), intent(in) :: rout
    real(wp) :: jacobian

    if (rout > 0.0_wp) then
      jacobian = rout
    else
      jacobian = 0.0_wp
    end if
  end function radial_jacobian

end module bh_toroid_radial_map_mod
