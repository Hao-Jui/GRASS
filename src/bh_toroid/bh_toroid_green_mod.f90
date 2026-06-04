module bh_toroid_green_mod
  use precision_mod, only: wp
  use bh_toroid_validation_mod, only: validation_result, validation_ok, validation_error, VALID_BAD_GREEN_ARGS
  implicit none
  private

  public :: validate_green_args
  public :: ne_f1_kernel, ne_f2_kernel
  public :: lambda_radial_kernel, b_radial_kernel, omega_radial_kernel

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

  pure function ne_f1_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value, image

    ! Nishida & Eriguchi (1994), eq. (3.4):
    ! f^1_n(r,r') = (r'/r)^n - h0^(2n)/(r r')^n, for r'/r <= 1,
    !              = (r/r')^n - h0^(2n)/(r r')^n, for r'/r > 1.
    image = (h0 * h0 / (r * rp))**n
    if (rp <= r) then
      value = (rp / r)**n - image
    else
      value = (r / rp)**n - image
    end if
  end function ne_f1_kernel

  pure function ne_f2_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value, image

    ! Nishida & Eriguchi (1994), eq. (3.5):
    ! f^2_n(r,r') = (1/r)  (r'/r)^n - h0^(2n+1)/(r r')^(n+1), for r'/r <= 1,
    !              = (1/r') (r/r')^n - h0^(2n+1)/(r r')^(n+1), for r'/r > 1.
    image = h0 * (h0 * h0 / (r * rp))**n / (r * rp)
    if (rp <= r) then
      value = (rp / r)**n / r - image
    else
      value = (r / rp)**n / rp - image
    end if
  end function ne_f2_kernel

  pure function lambda_radial_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value

    ! Equation (3.1) uses f^2 for the lambda = exp(nu) radial kernel.
    ! The caller supplies the already-selected angular/radial order.
    value = ne_f2_kernel(n, r, rp, h0)
  end function lambda_radial_kernel

  pure function b_radial_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value

    ! Equation (3.2) uses f^1 for the B = exp(gamma) radial kernel.
    ! The caller supplies the already-selected angular/radial order.
    value = ne_f1_kernel(n, r, rp, h0)
  end function b_radial_kernel

  pure function omega_radial_kernel(n, r, rp, h0) result(value)
    integer, intent(in) :: n
    real(wp), intent(in) :: r, rp, h0
    real(wp) :: value

    ! Equation (3.3) uses f^2 for the omega radial kernel.
    ! The caller supplies the already-selected angular/radial order.
    value = ne_f2_kernel(n, r, rp, h0)
  end function omega_radial_kernel

end module bh_toroid_green_mod
