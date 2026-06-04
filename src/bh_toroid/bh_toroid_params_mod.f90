module bh_toroid_params_mod
  use precision_mod, only: wp
  use bh_toroid_validation_mod, only: validation_result, validation_ok, validation_error, &
      VALID_BAD_HORIZON, VALID_BAD_MODEL_FAMILY, VALID_BAD_POLYTROPE, &
      VALID_BAD_RADIAL_ORDER, VALID_BAD_ROTATION, VALID_BAD_SCALE
  implicit none
  private

  ! Mirrors para_mod's inert model-family values without coupling this
  ! pure validation module to the global runtime state.
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

    params = bh_toroid_params()
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

end module bh_toroid_params_mod
