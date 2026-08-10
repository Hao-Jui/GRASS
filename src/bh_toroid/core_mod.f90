module validation_mod
  use precision_mod, only: wp
  implicit none
  private

  integer, parameter, public :: VALID_OK = 0
  integer, parameter, public :: VALID_BAD_MODEL_FAMILY = 1
  integer, parameter, public :: VALID_BAD_HORIZON = 2
  integer, parameter, public :: VALID_BAD_RADIAL_ORDER = 3
  integer, parameter, public :: VALID_BAD_SCALE = 4
  integer, parameter, public :: VALID_BAD_POLYTROPE = 5
  integer, parameter, public :: VALID_BAD_ROTATION = 6
  integer, parameter, public :: VALID_BAD_GRID_SIZE = 7
  integer, parameter, public :: VALID_BAD_GREEN_ARGS = 8

  type, public :: validation_result
    integer :: status = VALID_OK
    character(len=160) :: message = "ok"
  end type validation_result

  public :: validation_ok, validation_error, is_finite

contains

  pure function validation_ok() result(res)
    type(validation_result) :: res
  end function validation_ok

  pure function validation_error(status, message) result(res)
    integer, intent(in) :: status
    character(*), intent(in) :: message
    type(validation_result) :: res
    res%status = status
    res%message = message
  end function validation_error

  pure elemental function is_finite(value) result(ok)
    real(wp), intent(in) :: value
    logical :: ok
    ok = value == value .and. abs(value) < huge(value)
  end function is_finite

end module validation_mod


module params_mod
  use precision_mod, only: wp
  use validation_mod, only: validation_result, validation_ok, validation_error, &
      VALID_BAD_HORIZON, VALID_BAD_MODEL_FAMILY, VALID_BAD_POLYTROPE, &
      VALID_BAD_RADIAL_ORDER, VALID_BAD_ROTATION, VALID_BAD_SCALE
  implicit none
  private

  ! Mirrors para_panel's model-family values
  integer, parameter, public :: MODEL_NS_VALUE = 1
  integer, parameter, public :: MODEL_BH_TOROID_VALUE = 2

  type, public :: bh_toroid_params
    real(wp) :: h0_hat = 0.1_wp
    real(wp) :: omega_h = 0.0_wp
    real(wp) :: rin_hat = 0.4_wp
    real(wp) :: kappa_ratio = 0.1_wp
    real(wp) :: emax = 1.0_wp
    real(wp) :: poly_n = 3.0_wp
    real(wp) :: rotation_A = 1.0_wp
    real(wp) :: rout_scale = 1.0_wp
  end type bh_toroid_params

  public :: default_bh_toroid_params, validate_bh_toroid_params, validate_model_family

contains

  pure function default_bh_toroid_params() result(params)
    type(bh_toroid_params) :: params
  end function default_bh_toroid_params

  pure function validate_bh_toroid_params(params) result(res)
    type(bh_toroid_params), intent(in) :: params
    type(validation_result) :: res

    if (params%h0_hat <= 0.0_wp .or. params%h0_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_HORIZON, "h0_hat must satisfy 0 < h0_hat < 1")
    else if (params%rin_hat <= params%h0_hat .or. params%rin_hat >= 1.0_wp) then
      res = validation_error(VALID_BAD_RADIAL_ORDER, "radial domain must satisfy h0_hat < rin_hat < 1")
    else if (params%rout_scale <= 0.0_wp .or. params%emax <= 0.0_wp .or. params%kappa_ratio <= 0.0_wp) then
      res = validation_error(VALID_BAD_SCALE, "rout_scale, emax, and kappa_ratio must be positive")
    else if (params%poly_n <= 0.0_wp) then
      res = validation_error(VALID_BAD_POLYTROPE, "poly_n must be positive")
    else if (params%rotation_A <= 0.0_wp) then
      res = validation_error(VALID_BAD_ROTATION, "rotation_A must be positive")
    else
      res = validation_ok()
    end if
  end function validate_bh_toroid_params

  pure function validate_model_family(model_family) result(res)
    integer, intent(in) :: model_family
    type(validation_result) :: res

    select case (model_family)
    case (MODEL_NS_VALUE, MODEL_BH_TOROID_VALUE)
      res = validation_ok()
    case default
      res = validation_error(VALID_BAD_MODEL_FAMILY, "unknown model family")
    end select
  end function validate_model_family

end module params_mod


module radial_map_mod
  use precision_mod, only: wp
  use validation_mod, only: validation_result, validation_ok, validation_error, &
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
    end do
    weights(1:n) = dr
    weights(1) = 0.5_wp * dr
    weights(n) = 0.5_wp * dr
    res = validation_ok()
  end function build_rhat_trapezoid_grid

  pure elemental function classify_rhat(rhat, h0_hat, rin_hat) result(zone)
    real(wp), intent(in) :: rhat, h0_hat, rin_hat
    integer :: zone
    real(wp), parameter :: tol = 64.0_wp * epsilon(1.0_wp)

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
    jacobian = merge(rout, 0.0_wp, rout > 0.0_wp)
  end function radial_jacobian

end module radial_map_mod


module green_mod
  use precision_mod, only: wp
  use validation_mod, only: validation_result, validation_ok, validation_error, VALID_BAD_GREEN_ARGS
  implicit none
  private

  public :: validate_green_args, ne_f1_kernel, ne_f2_kernel

  ! Nishida & Eriguchi (1994):
  !   eq. (3.1) lambda = exp(nu)  -> f^2
  !   eq. (3.2) B = exp(gamma)    -> f^1
  !   eq. (3.3) omega             -> f^2

contains

  pure function validate_green_args(n, r, rp, h0) result(res)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    type(validation_result) :: res

    if (n < 0 .or. h0 < 0.0_wp .or. r <= h0 .or. rp <= h0 .or. r <= 0.0_wp .or. rp <= 0.0_wp) then
      res = validation_error(VALID_BAD_GREEN_ARGS, "Green-kernel arguments must satisfy n >= 0 and r, rp > h0 >= 0")
    else
      res = validation_ok()
    end if
  end function validate_green_args

  ! eq. (3.4): f^1_n(r,r') = (r</r>)^n - (h0^2/(r r'))^n
  pure function ne_f1_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value, ratio

    ratio = merge(rp / r, r / rp, rp <= r)
    value = ratio**n - (h0 * h0 / (r * rp))**n
  end function ne_f1_kernel

  ! eq. (3.5): f^2_n(r,r') = (1/r>) (r</r>)^n - h0^(2n+1)/(r r')^(n+1)
  pure function ne_f2_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value, ratio, inv_outer

    if (rp <= r) then
      ratio = rp / r
      inv_outer = 1.0_wp / r
    else
      ratio = r / rp
      inv_outer = 1.0_wp / rp
    end if
    value = ratio**n * inv_outer - h0 * (h0 * h0 / (r * rp))**n / (r * rp)
  end function ne_f2_kernel

end module green_mod
