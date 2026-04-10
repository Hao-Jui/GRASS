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
  real(wp), parameter, private :: ZEPS = 1.0e-8_wp

contains

subroutine brent_core(x_guess, scale_up, scale_down, tol, return_value, f, small_guess, logfile)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none

  real(wp), intent(in)  :: x_guess, scale_up, scale_down, tol
  real(wp), intent(out) :: return_value
  procedure(brent_func) :: f
  real(wp), intent(in), optional :: small_guess
  character(*), intent(in), optional :: logfile

  real(wp) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, xm, tol1, x0
  integer  :: iter
  logical  :: bracketed

  x0 = x_guess
  if (present(small_guess)) then
    if (abs(x0) < epsilon(x0)) x0 = small_guess
  end if

  a = x0
  b = x0
  bracketed = .false.

  ! ====================================================================
  ! PHASE 1: Adaptive geometric bracketing — scale outward until
  !          f(a) and f(b) have opposite signs.
  !          scale_up/scale_down are typically 1.1–2.0.
  ! ====================================================================
  bracket_loop: do iter = 1, MAX_ITER_BRACKET
    a = a * scale_up
    b = b / scale_down

    call f(a, fa)
    call f(b, fb)

    if (ieee_is_nan(fa)) then
      a = a / scale_up
      call f(a, fa)
    end if
    if (ieee_is_nan(fb)) then
      b = b * scale_down
      call f(b, fb)
    end if

    if (fa * fb <= 0.0_wp) then
      bracketed = .true.
      exit bracket_loop
    end if
  end do bracket_loop

  ! ====================================================================
  ! PHASE 2: Fallback — sweep a fine grid and log f(x) for diagnostics.
  !          If a sign change is found, use that bracket.
  ! ====================================================================
  if (.not. bracketed) then
    if (present(logfile)) then
      call write_bracket_log(logfile, x_guess, f, a, b, bracketed)
    end if
    if (.not. bracketed) then
      if (present(logfile)) write(*,*) "check  ", logfile
      error stop "brent_core: failed to bracket root"
    end if
    call f(a, fa)
    call f(b, fb)
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
    if (abs(xm) <= tol1 .or. abs(fb) < epsilon(fb)) then
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

    a = b; fa = fb
    b = b + merge(d, sign(tol1, xm), abs(d) > tol1)
    call f(b, fb)
  end do brent_loop

  error stop "brent_core: exceeding maximum iterations"

contains

  subroutine write_bracket_log(fname, xg, func, ax, bx, found)
    character(*), intent(in)  :: fname
    real(wp), intent(in)      :: xg
    procedure(brent_func)     :: func
    real(wp), intent(out)     :: ax, bx
    logical, intent(out)      :: found
    integer  :: u, i
    real(wp) :: cur_f, prev_f, cur_x

    found = .false.
    open(newunit=u, file=fname, status='replace', action='write')
    do i = 1, MAX_ITER_BRACKET
      cur_x = xg * 0.01_wp * i
      call func(cur_x, cur_f)
      if (i > 1 .and. cur_f * prev_f <= 0.0_wp) then
        bx = cur_x
        ax = xg * 0.01_wp * (i - 1)
        found = .true.
      end if
      prev_f = cur_f
      write(u, "(2ES15.6)") cur_x, cur_f
    end do
    close(u)
  end subroutine write_bracket_log

end subroutine brent_core

subroutine find_omege_e(x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol, return_value, f)
  implicit none
  real(wp), intent(in)  :: x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol
  real(wp), intent(out) :: return_value
  interface
    subroutine f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
      import :: wp
      real(wp), intent(in)  :: x, re, rho_h, g_h, w_h, rho_p, g_p
      real(wp), intent(out) :: fx
    end subroutine f
  end interface

  call brent_core(x_guess, 1.2_wp, 1.1_wp, tol, return_value, wrapped, &
                  logfile="./Cont/diff_rotation.dat")

contains
  subroutine wrapped(x, fx)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: fx
    call f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
  end subroutine wrapped

end subroutine find_omege_e

subroutine zbrent_rot(x_guess, re, rho_p, ww_p, sgp, mugp, tol, return_value, f)
  implicit none
  real(wp), intent(in)  :: x_guess, re, rho_p, ww_p, sgp, mugp, tol
  real(wp), intent(out) :: return_value
  interface
    subroutine f(x, fx, re, rho_p, ww_p, sgp, mugp)
      import :: wp
      real(wp), intent(in)  :: x, re, rho_p, ww_p, sgp, mugp
      real(wp), intent(out) :: fx
    end subroutine f
  end interface

  call brent_core(x_guess, 1.1_wp, 1.1_wp, tol, return_value, wrapped, &
                  small_guess=1.0e-4_wp, logfile="./Cont/rotation_law.dat")

contains
  subroutine wrapped(x, fx)
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: fx
    call f(x, fx, re, rho_p, ww_p, sgp, mugp)
  end subroutine wrapped

end subroutine zbrent_rot

end module brent_mod
