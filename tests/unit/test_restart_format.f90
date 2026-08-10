program test_restart_format
  use precision_mod, only: wp
  use iso_fortran_env, only: int32
  use restart_format_mod, only: &
      restart_meta_t, restart_status_t, &
      RESTART_MAGIC, RESTART_FORMAT_VERSION, RESTART_NFIELDS, &
      RESTART_OK, RESTART_BAD_MAGIC, RESTART_VERSION_MISMATCH, RESTART_BAD_HEADER, &
      pack_header_ints, unpack_header_ints, &
      pack_header_meta, unpack_header_meta, &
      validate_magic, validate_version, validate_storage_size, validate_field_count
  use test_utils
  implicit none

  call test_round_trip()
  call test_validate_magic()
  call test_validate_version()
  call test_validate_storage_size()
  call test_validate_field_count()
  call test_golden_binary()

  call test_summary("test_restart_format")

contains

  subroutine test_round_trip()
    type(restart_meta_t) :: original, recovered
    integer(int32) :: hi(6)
    real(wp) :: hm(5)

    original%sdiv = 401
    original%mdiv = 65
    original%spwr = 2
    original%r_e      = 1.234_wp
    original%e_center = 5.6789e14_wp
    original%r_ratio  = 0.7_wp
    original%omega_e  = 0.0123_wp
    original%omega_c  = 0.0456_wp

    hi = pack_header_ints(original)
    hm = pack_header_meta(original)
    call unpack_header_ints(hi, recovered)
    call unpack_header_meta(hm, recovered)

    call assert_true("round-trip sdiv",     recovered%sdiv == original%sdiv)
    call assert_true("round-trip mdiv",     recovered%mdiv == original%mdiv)
    call assert_true("round-trip spwr",     recovered%spwr == original%spwr)
    call assert_near ("round-trip r_e",      original%r_e,      recovered%r_e,      0._wp)
    call assert_near ("round-trip e_center", original%e_center, recovered%e_center, 0._wp)
    call assert_near ("round-trip r_ratio",  original%r_ratio,  recovered%r_ratio,  0._wp)
    call assert_near ("round-trip omega_e",  original%omega_e,  recovered%omega_e,  0._wp)
    call assert_near ("round-trip omega_c",  original%omega_c,  recovered%omega_c,  0._wp)
  end subroutine

  subroutine test_validate_magic()
    type(restart_status_t) :: s
    s = validate_magic(RESTART_MAGIC)
    call assert_true("magic ok",  s%code == RESTART_OK)
    s = validate_magic("WRONGMAGIC")
    call assert_true("magic bad", s%code == RESTART_BAD_MAGIC)
  end subroutine

  subroutine test_validate_version()
    type(restart_status_t) :: s
    s = validate_version(RESTART_FORMAT_VERSION)
    call assert_true("version ok", s%code == RESTART_OK)
    s = validate_version(RESTART_FORMAT_VERSION + 1_int32)
    call assert_true("version mismatch code",     s%code == RESTART_VERSION_MISMATCH)
    call assert_true("version mismatch expected", s%expected == RESTART_FORMAT_VERSION)
    call assert_true("version mismatch found",    s%found == RESTART_FORMAT_VERSION + 1_int32)
  end subroutine

  subroutine test_validate_storage_size()
    type(restart_status_t) :: s
    integer(int32) :: bits_wp
    bits_wp = int(storage_size(1.0_wp), int32)
    s = validate_storage_size(bits_wp)
    call assert_true("storage_size ok",  s%code == RESTART_OK)
    s = validate_storage_size(bits_wp + 1_int32)
    call assert_true("storage_size bad", s%code == RESTART_BAD_HEADER)
  end subroutine

  subroutine test_validate_field_count()
    type(restart_status_t) :: s
    s = validate_field_count(int(RESTART_NFIELDS, int32))
    call assert_true("nfields ok",  s%code == RESTART_OK)
    s = validate_field_count(int(RESTART_NFIELDS + 1, int32))
    call assert_true("nfields bad", s%code == RESTART_BAD_HEADER)
  end subroutine

  subroutine test_golden_binary()
    character(len=*), parameter :: golden_path = "tests/reference/restart_v1_header.bin"
    integer :: unit, ios
    integer(int32) :: hi(6)
    real(wp) :: hm(5)
    character(len=len(RESTART_MAGIC)) :: magic
    type(restart_meta_t) :: meta
    type(restart_status_t) :: s
    logical :: exists

    inquire(file=golden_path, exist=exists)
    if (.not. exists) then
      call write_golden(golden_path)
    end if

    open(newunit=unit, file=golden_path, status="old", action="read", &
         access="stream", form="unformatted", iostat=ios)
    call assert_true("golden open", ios == 0)
    if (ios /= 0) return

    read(unit, iostat=ios) magic
    call assert_true("golden magic read", ios == 0)
    s = validate_magic(magic)
    call assert_true("golden magic validates", s%code == RESTART_OK)

    read(unit, iostat=ios) hi
    call assert_true("golden hi read", ios == 0)
    s = validate_version(hi(1))
    call assert_true("golden version validates", s%code == RESTART_OK)
    s = validate_storage_size(hi(2))
    call assert_true("golden storage_size validates", s%code == RESTART_OK)
    s = validate_field_count(hi(3))
    call assert_true("golden field count validates", s%code == RESTART_OK)

    read(unit, iostat=ios) hm
    call assert_true("golden hm read", ios == 0)

    call unpack_header_ints(hi, meta)
    call unpack_header_meta(hm, meta)

    call assert_true("golden sdiv",     meta%sdiv == 401)
    call assert_true("golden mdiv",     meta%mdiv == 65)
    call assert_true("golden spwr",     meta%spwr == 2)
    call assert_near ("golden r_e",      1.234_wp,       meta%r_e,      0._wp)
    call assert_near ("golden e_center", 5.6789e14_wp,   meta%e_center, 0._wp)
    call assert_near ("golden r_ratio",  0.7_wp,         meta%r_ratio,  0._wp)
    call assert_near ("golden omega_e",  0.0123_wp,      meta%omega_e,  0._wp)
    call assert_near ("golden omega_c",  0.0456_wp,      meta%omega_c,  0._wp)
    close(unit)
  end subroutine

  subroutine write_golden(path)
    character(len=*), intent(in) :: path
    type(restart_meta_t) :: meta
    integer :: unit, ios
    integer(int32) :: hi(6)
    real(wp) :: hm(5)

    meta%sdiv = 401
    meta%mdiv = 65
    meta%spwr = 2
    meta%r_e      = 1.234_wp
    meta%e_center = 5.6789e14_wp
    meta%r_ratio  = 0.7_wp
    meta%omega_e  = 0.0123_wp
    meta%omega_c  = 0.0456_wp

    hi = pack_header_ints(meta)
    hm = pack_header_meta(meta)

    open(newunit=unit, file=path, status="replace", action="write", &
         access="stream", form="unformatted", iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_golden: failed to open ", trim(path)
      return
    end if
    write(unit) RESTART_MAGIC
    write(unit) hi
    write(unit) hm
    close(unit)
    write(*,*) "wrote golden fixture: ", trim(path)
  end subroutine

end program test_restart_format
