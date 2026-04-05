module brent_mod
  use precision_mod, only: wp
  implicit none
  public :: find_omege_e, zbrent_rot

  abstract interface
    subroutine brent_func(x, fx)
      import :: wp
      real(wp), intent(in)  :: x
      real(wp), intent(out) :: fx
    end subroutine brent_func
  end interface

  integer, parameter, private :: MAX_ITER_BRACKET = 200
  integer, parameter, private :: MAX_ITER_BRENT   = 100
  integer, parameter, private :: BRACKET_LOG_SAMPLES = 200
  integer, parameter, private :: ROOT_SUCCESS = 0
  integer, parameter, private :: ROOT_BRACKET_FAIL = 1
  integer, parameter, private :: ROOT_ITER_FAIL = 2
  integer, parameter, private :: ROOT_BAD_INPUT = 3
  real(wp), parameter, private :: ZEPS = 1.0e-8_wp
  real(wp), parameter, private :: MIN_POSITIVE_GUESS = 1.0e-12_wp

contains

subroutine brent_core(x_guess, scale_up, scale_down, tol, return_value, f, ierr, errmsg, small_guess, logfile)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(wp), intent(in)  :: x_guess, scale_up, scale_down, tol
  real(wp), intent(out) :: return_value
  procedure(brent_func) :: f
  integer, intent(out) :: ierr
  character(len=*), intent(out), optional :: errmsg
  real(wp), intent(in), optional :: small_guess
  character(*), intent(in), optional :: logfile

  real(wp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, xm, tol1, x0
  real(wp) :: candidate, f_candidate
  integer  :: iter
  logical  :: bracketed, ok_eval

  ierr = ROOT_SUCCESS
  if (present(errmsg)) errmsg = ""
  return_value = x_guess

  if (tol <= 0.0_wp .or. scale_up <= 1.0_wp .or. scale_down <= 1.0_wp) then
    ierr = ROOT_BAD_INPUT
    if (present(errmsg)) errmsg = "brent_core: invalid solver parameters"
    return
  end if

  x0 = abs(x_guess)
  if (present(small_guess)) x0 = max(x0, abs(small_guess))
  x0 = max(x0, max(tol, MIN_POSITIVE_GUESS))

  call find_bracket(x0, scale_up, scale_down, a, b, fa, fb, bracketed)
  if (.not. bracketed) then
    if (present(logfile)) call write_bracket_log(logfile, x0, a, b, f)
    ierr = ROOT_BRACKET_FAIL
    if (present(errmsg)) errmsg = "brent_core: failed to bracket root"
    return
  end if

  ! ====================================================================
  ! PHASE 3: Brent's method — superlinear convergence with guaranteed
  !          bisection fallback.  Variables follow NR convention:
  !            b  = current best estimate (|fb| <= |fc|)
  !            c  = contrapoint (opposite sign to fb)
  !            a  = previous iterate
  !            d  = proposed step, e = step before last
  ! ====================================================================
  c = b; fc = fb; e = 0.0_wp; d = 0.0_wp

  brent_loop: do iter = 1, MAX_ITER_BRENT
    ! Ensure |fb| <= |fc| for interpolation stability
    if (fb * fc > 0.0_wp) then
      c = a; fc = fa; d = b - a; e = d
    end if
    if (abs(fc) < abs(fb)) then
      a = b; b = c; c = a
      fa = fb; fb = fc; fc = fa
    end if

    tol1 = 2.0_wp * ZEPS * abs(b) + 0.5_wp * tol
    xm = 0.5_wp * (c - b)

    ! Convergence: bracket width or residual below tolerance
    if (abs(xm) <= tol1 .or. fb == 0.0_wp) then
      return_value = b
      return
    end if

    ! Attempt inverse quadratic / linear interpolation
    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s = fb / fa
      if (abs(a - c) < spacing(a)) then
        ! Linear (secant) interpolation
        p = 2.0_wp * xm * s
        q = 1.0_wp - s
      else
        ! Inverse quadratic interpolation
        q = fa / fc; r = fb / fc
        p = s * (2.0_wp * xm * q * (q - r) - (b - a) * (r - 1.0_wp))
        q = (q - 1.0_wp) * (r - 1.0_wp) * (s - 1.0_wp)
      end if
      if (p > 0.0_wp) q = -q
      p = abs(p)
      ! Accept interpolation only if step is small enough
      if (2.0_wp * p < min(3.0_wp * xm * q - abs(tol1 * q), abs(e * q))) then
        e = d; d = p / q
      else
        d = xm; e = d
      end if
    else
      ! Bisection: bounds decreasing too slowly
      d = xm; e = d
    end if

    a = b
    fa = fb
    candidate = b + merge(d, sign(tol1, xm), abs(d) > tol1)
    call evaluate_point(candidate, f_candidate, ok_eval)
    if (.not. ok_eval) then
      candidate = b + merge(xm, sign(tol1, xm), abs(xm) > tol1)
      call evaluate_point(candidate, f_candidate, ok_eval)
      if (.not. ok_eval) then
        ierr = ROOT_ITER_FAIL
        if (present(errmsg)) errmsg = "brent_core: function evaluation failed inside bracket"
        return_value = b
        return
      end if
      d = xm
      e = d
    end if
    b = candidate
    fb = f_candidate
  end do brent_loop

  ierr = ROOT_ITER_FAIL
  if (present(errmsg)) errmsg = "brent_core: exceeded maximum iterations"
  return_value = b

contains
  subroutine evaluate_point(x, fx, ok)
    real(wp), intent(in) :: x
    real(wp), intent(out) :: fx
    logical, intent(out) :: ok

    call f(x, fx)
    ok = ieee_is_finite(fx)
  end subroutine evaluate_point

  subroutine find_bracket(x_seed, grow_up, grow_down, ax, bx, fax, fbx, found)
    real(wp), intent(in) :: x_seed, grow_up, grow_down
    real(wp), intent(out) :: ax, bx, fax, fbx
    logical, intent(out) :: found
    real(wp) :: x_left, x_mid, x_right, f_left, f_mid, f_right
    logical :: ok_left, ok_mid, ok_right
    integer :: k

    found = .false.
    x_mid = x_seed
    call evaluate_point(x_mid, f_mid, ok_mid)
    if (ok_mid .and. f_mid == 0.0_wp) then
      ax = x_mid; bx = x_mid; fax = f_mid; fbx = f_mid
      found = .true.
      return
    end if

    x_left = max(x_seed / grow_down, MIN_POSITIVE_GUESS)
    x_right = max(x_seed * grow_up, x_seed + MIN_POSITIVE_GUESS)

    do k = 1, MAX_ITER_BRACKET
      call evaluate_point(x_left, f_left, ok_left)
      if (ok_left .and. f_left == 0.0_wp) then
        ax = x_left; bx = x_left; fax = f_left; fbx = f_left
        found = .true.
        return
      end if
      if (ok_left .and. ok_mid .and. f_left * f_mid <= 0.0_wp) then
        ax = x_left; bx = x_mid; fax = f_left; fbx = f_mid
        found = .true.
        return
      end if

      call evaluate_point(x_right, f_right, ok_right)
      if (ok_right .and. f_right == 0.0_wp) then
        ax = x_right; bx = x_right; fax = f_right; fbx = f_right
        found = .true.
        return
      end if
      if (ok_mid .and. ok_right .and. f_mid * f_right <= 0.0_wp) then
        ax = x_mid; bx = x_right; fax = f_mid; fbx = f_right
        found = .true.
        return
      end if
      if (ok_left .and. ok_right .and. f_left * f_right <= 0.0_wp) then
        ax = x_left; bx = x_right; fax = f_left; fbx = f_right
        found = .true.
        return
      end if

      x_left = max(x_left / grow_down, MIN_POSITIVE_GUESS)
      x_right = x_right * grow_up
    end do

    ax = x_left
    bx = x_right
    fax = 0.0_wp
    fbx = 0.0_wp
  end subroutine find_bracket

  subroutine write_bracket_log(fname, xg, ax, bx, func)
    character(*), intent(in)  :: fname
    real(wp), intent(in)      :: xg
    real(wp), intent(in)      :: ax, bx
    procedure(brent_func)     :: func
    integer  :: u, i
    real(wp) :: cur_f, cur_x, x_lo, x_hi

    x_lo = max(MIN_POSITIVE_GUESS, min(ax, bx, xg))
    x_hi = max(x_lo * 1.0001_wp, max(ax, bx, xg))
    open(newunit=u, file=fname, status='replace', action='write')
    do i = 1, BRACKET_LOG_SAMPLES
      cur_x = x_lo + (x_hi - x_lo) * real(i - 1, wp) / real(BRACKET_LOG_SAMPLES - 1, wp)
      call func(cur_x, cur_f)
      write(u, "(2ES15.6)") cur_x, cur_f
    end do
    close(u)
  end subroutine write_bracket_log

end subroutine brent_core

subroutine find_omege_e(x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol, return_value, f, ierr, errmsg)
  implicit none
  real(wp), intent(in)  :: x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol
  real(wp), intent(out) :: return_value
  integer, intent(out), optional :: ierr
  character(len=*), intent(out), optional :: errmsg
  interface
    subroutine f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
      import :: wp
      real(wp), intent(in)  :: x, re, rho_h, g_h, w_h, rho_p, g_p
      real(wp), intent(out) :: fx
    end subroutine f
  end interface
  integer :: status
  character(len=256) :: message

  call brent_core(x_guess, 1.2_wp, 1.1_wp, tol, return_value, wrapped, status, message, &
                  small_guess=1.0e-4_wp, logfile="./Cont/diff_rotation.dat")
  call finalize_wrapper_status("find_omege_e", status, message, ierr, errmsg)

contains
  subroutine wrapped(x, fx)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: fx
    call f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
  end subroutine wrapped

end subroutine find_omege_e

subroutine zbrent_rot(x_guess, re, rho_p, ww_p, sgp, mugp, tol, return_value, f, ierr, errmsg)
  implicit none
  real(wp), intent(in)  :: x_guess, re, rho_p, ww_p, sgp, mugp, tol
  real(wp), intent(out) :: return_value
  integer, intent(out), optional :: ierr
  character(len=*), intent(out), optional :: errmsg
  interface
    subroutine f(x, fx, re, rho_p, ww_p, sgp, mugp)
      import :: wp
      real(wp), intent(in)  :: x, re, rho_p, ww_p, sgp, mugp
      real(wp), intent(out) :: fx
    end subroutine f
  end interface
  integer :: status
  character(len=256) :: message

  call brent_core(x_guess, 1.1_wp, 1.1_wp, tol, return_value, wrapped, status, message, &
                  small_guess=1.0e-4_wp, logfile="./Cont/rotation_law.dat")
  call finalize_wrapper_status("zbrent_rot", status, message, ierr, errmsg)

contains
  subroutine wrapped(x, fx)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: fx
    call f(x, fx, re, rho_p, ww_p, sgp, mugp)
  end subroutine wrapped

end subroutine zbrent_rot

subroutine finalize_wrapper_status(routine_name, status, message, ierr, errmsg)
  character(len=*), intent(in) :: routine_name, message
  integer, intent(in) :: status
  integer, intent(out), optional :: ierr
  character(len=*), intent(out), optional :: errmsg

  if (present(ierr)) ierr = status
  if (present(errmsg)) errmsg = trim(message)
  if (status == ROOT_SUCCESS) return
  if (.not. present(ierr)) then
    write(*,'(A,": ",A)') trim(routine_name), trim(message)
    stop trim(routine_name)//": root solve failed"
  end if
end subroutine finalize_wrapper_status

end module brent_mod
