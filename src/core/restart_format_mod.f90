module restart_format_mod
  use precision_mod, only: wp
  use iso_fortran_env, only: int32
  implicit none
  private

  character(len=*), parameter, public :: RESTART_MAGIC = "GRASSRST01"
  integer(int32), parameter, public :: RESTART_FORMAT_VERSION = 1_int32
  integer, parameter, public :: RESTART_NFIELDS = 10

  integer, parameter, public :: F_ALPHA = 1, F_GAMA = 2, F_RHO = 3, F_WW = 4
  integer, parameter, public :: F_PRESSURE = 5, F_ENERGY = 6, F_ENTHALPY = 7
  integer, parameter, public :: F_VELOCITY_SQ = 8, F_OMG = 9, F_SPHI = 10

  integer, parameter, public :: RESTART_OK = 0
  integer, parameter, public :: RESTART_BAD_MAGIC = 1
  integer, parameter, public :: RESTART_VERSION_MISMATCH = 2
  integer, parameter, public :: RESTART_BAD_HEADER = 3

  type, public :: restart_meta_t
    integer  :: sdiv = 0, mdiv = 0, spwr = 0
    real(wp) :: r_e = 0._wp, e_center = 0._wp, r_ratio = 0._wp
    real(wp) :: omega_e = 0._wp, omega_c = 0._wp
  end type

  type, public :: restart_status_t
    integer :: code = RESTART_OK
    integer :: expected = 0
    integer :: found = 0
  end type

  public :: pack_header_ints, unpack_header_ints
  public :: pack_header_meta, unpack_header_meta
  public :: validate_magic, validate_version, validate_storage_size, validate_field_count

contains

  pure function pack_header_ints(meta) result(arr)
    type(restart_meta_t), intent(in) :: meta
    integer(int32) :: arr(6)
    arr = [RESTART_FORMAT_VERSION, &
           int(storage_size(1.0_wp), int32), &
           int(RESTART_NFIELDS, int32), &
           int(meta%sdiv, int32), &
           int(meta%mdiv, int32), &
           int(meta%spwr, int32)]
  end function

  pure function pack_header_meta(meta) result(arr)
    type(restart_meta_t), intent(in) :: meta
    real(wp) :: arr(5)
    arr = [meta%r_e, meta%e_center, meta%r_ratio, meta%omega_e, meta%omega_c]
  end function

  pure subroutine unpack_header_ints(arr, meta)
    integer(int32), intent(in) :: arr(6)
    type(restart_meta_t), intent(inout) :: meta
    meta%sdiv = arr(4)
    meta%mdiv = arr(5)
    meta%spwr = arr(6)
  end subroutine

  pure subroutine unpack_header_meta(arr, meta)
    real(wp), intent(in) :: arr(5)
    type(restart_meta_t), intent(inout) :: meta
    meta%r_e      = arr(1)
    meta%e_center = arr(2)
    meta%r_ratio  = arr(3)
    meta%omega_e  = arr(4)
    meta%omega_c  = arr(5)
  end subroutine

  pure function validate_magic(magic) result(status)
    character(len=*), intent(in) :: magic
    type(restart_status_t) :: status
    if (magic == RESTART_MAGIC) then
      status%code = RESTART_OK
    else
      status%code = RESTART_BAD_MAGIC
    end if
  end function

  pure function validate_version(version) result(status)
    integer(int32), intent(in) :: version
    type(restart_status_t) :: status
    status%expected = RESTART_FORMAT_VERSION
    status%found = version
    if (version == RESTART_FORMAT_VERSION) then
      status%code = RESTART_OK
    else
      status%code = RESTART_VERSION_MISMATCH
    end if
  end function

  pure function validate_storage_size(bits) result(status)
    integer(int32), intent(in) :: bits
    type(restart_status_t) :: status
    status%expected = int(storage_size(1.0_wp), int32)
    status%found = bits
    if (bits == status%expected) then
      status%code = RESTART_OK
    else
      status%code = RESTART_BAD_HEADER
    end if
  end function

  pure function validate_field_count(n) result(status)
    integer(int32), intent(in) :: n
    type(restart_status_t) :: status
    status%expected = int(RESTART_NFIELDS, int32)
    status%found = n
    if (n == status%expected) then
      status%code = RESTART_OK
    else
      status%code = RESTART_BAD_HEADER
    end if
  end function

end module restart_format_mod
