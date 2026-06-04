module bh_toroid_validation_mod
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

  public :: validation_ok, validation_error

contains

  pure function validation_ok() result(res)
    type(validation_result) :: res

    res%status = VALID_OK
    res%message = "ok"
  end function validation_ok

  pure function validation_error(status, message) result(res)
    integer, intent(in) :: status
    character(*), intent(in) :: message
    type(validation_result) :: res

    res%status = status
    res%message = message
  end function validation_error

end module bh_toroid_validation_mod
