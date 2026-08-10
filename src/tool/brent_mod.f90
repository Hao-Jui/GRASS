module brent_mod
  use precision_mod, only: wp
  implicit none
  public :: find_omega_e, zbrent_rot, find_omega_e_admissible, zbrent_rot_admissible
  private :: signs_differ

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

pure logical function signs_differ(a, b)
  real(wp), intent(in) :: a, b
  signs_differ = (a < 0._wp .and. b >= 0._wp) .or. (a >= 0._wp .and. b < 0._wp)
end function signs_differ

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

subroutine find_omega_e(x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol, return_value, f)
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

end subroutine find_omega_e

subroutine find_omega_e_admissible(x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol, return_value, f, success)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
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
  logical, intent(out), optional :: success

  integer, parameter :: N_SCAN = 512
  integer :: i
  real(wp) :: speed_scale, margin, lower_bound, upper_bound, x, fx
  real(wp) :: previous_x, previous_fx, bracket_a, bracket_b, bracket_distance
  real(wp) :: candidate_distance
  logical :: previous_valid, bracket_found, converged

  if (present(success)) success = .false.
  return_value = ieee_value(0._wp, ieee_quiet_nan)
  speed_scale = exp(re**2*rho_h)
  if (.not. ieee_is_finite(speed_scale) .or. speed_scale <= 0._wp .or. &
      .not. ieee_is_finite(tol) .or. tol <= 0._wp) then
    call refuse("find_omega_e_admissible: invalid search interval")
    return
  end if

  margin = max(32._wp*epsilon(speed_scale)*max(1._wp, abs(w_h), speed_scale), 0.01_wp*tol)
  lower_bound = w_h + margin
  upper_bound = w_h + speed_scale - margin
  if (.not. ieee_is_finite(lower_bound) .or. .not. ieee_is_finite(upper_bound) .or. &
      lower_bound >= upper_bound) then
    call refuse("find_omega_e_admissible: empty admissible interval")
    return
  end if

  previous_valid = .false.
  bracket_found = .false.
  bracket_distance = huge(bracket_distance)
  do i = 0, N_SCAN
    x = lower_bound + (upper_bound-lower_bound)*real(i, wp)/real(N_SCAN, wp)
    call f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    if (.not. ieee_is_finite(fx)) then
      previous_valid = .false.
      cycle
    end if
    if (abs(fx) <= epsilon(fx)) then
      return_value = x
      if (present(success)) success = .true.
      return
    end if
    if (previous_valid) then
      if (signs_differ(previous_fx, fx)) then
        candidate_distance = abs(0.5_wp*(previous_x+x)-x_guess)
        if (.not. bracket_found .or. candidate_distance < bracket_distance) then
          bracket_a = previous_x
          bracket_b = x
          bracket_distance = candidate_distance
          bracket_found = .true.
        end if
      end if
    end if
    previous_x = x
    previous_fx = fx
    previous_valid = .true.
  end do

  if (.not. bracket_found) then
    call refuse("find_omega_e_admissible: no root in admissible interval")
    return
  end if

  call bisect_admissible(bracket_a, bracket_b, return_value, converged)
  if (.not. converged) then
    call refuse("find_omega_e_admissible: admissible root did not converge")
    return
  end if
  if (present(success)) success = .true.

contains
  subroutine refuse(message)
    character(*), intent(in) :: message
    return_value = ieee_value(0._wp, ieee_quiet_nan)
    if (.not. present(success)) error stop message
  end subroutine refuse

  subroutine bisect_admissible(left, right, root, root_converged)
    real(wp), intent(in) :: left, right
    real(wp), intent(out) :: root
    logical, intent(out) :: root_converged
    real(wp) :: a, b, midpoint, fa, fb, fmid
    integer :: iteration

    a = left
    b = right
    call f(a, fa, re, rho_h, g_h, w_h, rho_p, g_p)
    call f(b, fb, re, rho_h, g_h, w_h, rho_p, g_p)
    root_converged = .false.
    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb) .or. &
        .not. signs_differ(fa, fb)) return

    do iteration = 1, MAX_ITER_BRENT
      midpoint = a + 0.5_wp*(b-a)
      call f(midpoint, fmid, re, rho_h, g_h, w_h, rho_p, g_p)
      if (.not. ieee_is_finite(fmid)) return
      if (abs(fmid) <= epsilon(fmid) .or. 0.5_wp*abs(b-a) <= tol) then
        root = midpoint
        root_converged = .true.
        return
      end if
      if (signs_differ(fa, fmid)) then
        b = midpoint
        fb = fmid
      else
        a = midpoint
        fa = fmid
      end if
    end do
  end subroutine bisect_admissible

end subroutine find_omega_e_admissible

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

subroutine zbrent_rot_admissible(x_guess, re, rho_p, ww_p, sgp, mugp, tol, return_value, f, success)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
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
  logical, intent(out), optional :: success

  real(wp), parameter :: BRACKET_SCALE = 1.2_wp
  real(wp) :: sin2m, speed_limit, margin, lower_delta, upper_delta, guess_delta
  real(wp) :: left_x, right_x, left_fx, right_fx, new_x, new_fx
  real(wp) :: bracket_a, bracket_b
  integer :: iteration
  logical :: left_open, right_open, bracket_found, converged

  if (present(success)) success = .false.
  return_value = ieee_value(0._wp, ieee_quiet_nan)
  sin2m = 1._wp - mugp**2
  if (.not. ieee_is_finite(sin2m) .or. sin2m <= 0._wp .or. sgp <= 0._wp .or. sgp >= 1._wp .or. &
      .not. ieee_is_finite(tol) .or. tol <= 0._wp) then
    call refuse("zbrent_rot_admissible: invalid rotation-law search interval")
    return
  end if
  speed_limit = exp(re**2*rho_p)*(1._wp-sgp)/(sgp*sqrt(sin2m))
  margin = max(32._wp*epsilon(speed_limit)*max(1._wp, abs(ww_p), speed_limit), 0.01_wp*tol)
  lower_delta = margin
  upper_delta = speed_limit - margin
  if (.not. ieee_is_finite(speed_limit) .or. lower_delta >= upper_delta) then
    call refuse("zbrent_rot_admissible: empty rotation-law search interval")
    return
  end if

  guess_delta = min(max(x_guess-ww_p, max(1.0e-4_wp, lower_delta)), upper_delta)
  left_x = ww_p + guess_delta
  right_x = left_x
  call f(left_x, left_fx, re, rho_p, ww_p, sgp, mugp)
  if (.not. ieee_is_finite(left_fx)) then
    call refuse("zbrent_rot_admissible: initial trial is inadmissible")
    return
  end if
  if (abs(left_fx) <= epsilon(left_fx)) then
    return_value = left_x
    if (present(success)) success = .true.
    return
  end if

  right_fx = left_fx
  left_open = guess_delta > lower_delta
  right_open = guess_delta < upper_delta
  bracket_found = .false.
  do iteration = 1, MAX_ITER_BRACKET
    if (left_open) then
      new_x = ww_p + max(lower_delta, (left_x-ww_p)/BRACKET_SCALE)
      left_open = new_x > ww_p + lower_delta
      call f(new_x, new_fx, re, rho_p, ww_p, sgp, mugp)
      if (.not. ieee_is_finite(new_fx)) then
        call refuse("zbrent_rot_admissible: lower trial left admissible domain")
        return
      end if
      if (signs_differ(new_fx, left_fx)) then
        bracket_a = new_x
        bracket_b = left_x
        bracket_found = .true.
        exit
      end if
      left_x = new_x
      left_fx = new_fx
    end if

    if (right_open) then
      new_x = ww_p + min(upper_delta, (right_x-ww_p)*BRACKET_SCALE)
      right_open = new_x < ww_p + upper_delta
      call f(new_x, new_fx, re, rho_p, ww_p, sgp, mugp)
      if (.not. ieee_is_finite(new_fx)) then
        call refuse("zbrent_rot_admissible: upper trial left admissible domain")
        return
      end if
      if (signs_differ(right_fx, new_fx)) then
        bracket_a = right_x
        bracket_b = new_x
        bracket_found = .true.
        exit
      end if
      right_x = new_x
      right_fx = new_fx
    end if
    if (.not. left_open .and. .not. right_open) exit
  end do

  if (.not. bracket_found) then
    call refuse("zbrent_rot_admissible: no root in admissible interval")
    return
  end if
  call bisect_admissible(bracket_a, bracket_b, return_value, converged)
  if (.not. converged) then
    call refuse("zbrent_rot_admissible: admissible root did not converge")
    return
  end if
  if (present(success)) success = .true.

contains
  subroutine refuse(message)
    character(*), intent(in) :: message
    return_value = ieee_value(0._wp, ieee_quiet_nan)
    if (.not. present(success)) error stop message
  end subroutine refuse

  subroutine bisect_admissible(left, right, root, root_converged)
    real(wp), intent(in) :: left, right
    real(wp), intent(out) :: root
    logical, intent(out) :: root_converged
    real(wp) :: a, b, midpoint, fa, fb, fmid
    integer :: bisect_iteration

    a = left
    b = right
    call f(a, fa, re, rho_p, ww_p, sgp, mugp)
    call f(b, fb, re, rho_p, ww_p, sgp, mugp)
    root_converged = .false.
    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb) .or. &
        .not. signs_differ(fa, fb)) return
    do bisect_iteration = 1, MAX_ITER_BRENT
      midpoint = a + 0.5_wp*(b-a)
      call f(midpoint, fmid, re, rho_p, ww_p, sgp, mugp)
      if (.not. ieee_is_finite(fmid)) return
      if (abs(fmid) <= epsilon(fmid) .or. 0.5_wp*abs(b-a) <= tol) then
        root = midpoint
        root_converged = .true.
        return
      end if
      if (signs_differ(fa, fmid)) then
        b = midpoint
        fb = fmid
      else
        a = midpoint
        fa = fmid
      end if
    end do
  end subroutine bisect_admissible

end subroutine zbrent_rot_admissible

end module brent_mod
