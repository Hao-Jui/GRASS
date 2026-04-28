program test_bh_toroid_params
  use precision_mod, only: wp
  use para_mod, only: MODEL_NS, MODEL_BH_TOROID, model_family
  use bh_toroid_params_mod, only: bh_toroid_params, default_bh_toroid_params, &
      validate_bh_toroid_params, validate_model_family
  use bh_toroid_validation_mod, only: validation_result, VALID_OK, &
      VALID_BAD_HORIZON, VALID_BAD_MODEL_FAMILY, VALID_BAD_POLYTROPE, &
      VALID_BAD_RADIAL_ORDER, VALID_BAD_ROTATION, VALID_BAD_SCALE
  use test_utils
  implicit none

  type(bh_toroid_params) :: params
  type(validation_result) :: res

  params = default_bh_toroid_params()
  res = validate_bh_toroid_params(params)
  call assert_true("default parameters valid", res%status == VALID_OK)
  call assert_true("default model family remains neutron star", model_family == MODEL_NS)
  call assert_status("BH model family accepted", validate_model_family(MODEL_BH_TOROID), VALID_OK)

  params = bh_toroid_params(h0_hat=0.15_wp, omega_h=0.02_wp, rin_hat=0.45_wp, &
      kappa_ratio=0.2_wp, emax=2.0_wp, poly_n=2.5_wp, rotation_A=1.3_wp, rout_scale=4.0_wp)
  call assert_status("sample parameters valid", validate_bh_toroid_params(params), VALID_OK)

  params = default_bh_toroid_params(); params%h0_hat = 0.0_wp
  call assert_status("bad horizon status", validate_bh_toroid_params(params), VALID_BAD_HORIZON)

  params = default_bh_toroid_params(); params%h0_hat = params%rin_hat
  call assert_status("bad radial order status", validate_bh_toroid_params(params), VALID_BAD_RADIAL_ORDER)

  params = default_bh_toroid_params(); params%rin_hat = 1.0_wp
  call assert_status("bad outer radial order status", validate_bh_toroid_params(params), VALID_BAD_RADIAL_ORDER)

  params = default_bh_toroid_params(); params%emax = 0.0_wp
  call assert_status("bad emax status", validate_bh_toroid_params(params), VALID_BAD_SCALE)

  params = default_bh_toroid_params(); params%kappa_ratio = 0.0_wp
  call assert_status("bad kappa status", validate_bh_toroid_params(params), VALID_BAD_SCALE)

  params = default_bh_toroid_params(); params%rout_scale = 0.0_wp
  call assert_status("bad rout status", validate_bh_toroid_params(params), VALID_BAD_SCALE)

  params = default_bh_toroid_params(); params%poly_n = 0.0_wp
  call assert_status("bad polytrope status", validate_bh_toroid_params(params), VALID_BAD_POLYTROPE)

  params = default_bh_toroid_params(); params%rotation_A = 0.0_wp
  call assert_status("bad rotation status", validate_bh_toroid_params(params), VALID_BAD_ROTATION)

  call assert_status("bad model family status", validate_model_family(99), VALID_BAD_MODEL_FAMILY)

  call test_summary("test_bh_toroid_params")

contains

  subroutine assert_status(label, result, expected)
    character(*), intent(in) :: label
    type(validation_result), intent(in) :: result
    integer, intent(in) :: expected

    call assert_true(label, result%status == expected)
  end subroutine assert_status
end program test_bh_toroid_params
