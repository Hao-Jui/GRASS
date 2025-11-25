subroutine zbrent_diff(x_guess, re, rho_h, g_h, w_h, rho_p, g_p, dh, tol, return_value, f)
  ! Brent-style root finder with adaptive bracketing for the differential rotation solver.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  implicit none
  real(8), intent(in) :: x_guess, re, rho_h, g_h, w_h, rho_p, g_p, dh, tol
  real(8), intent(out) :: return_value
  real(8), parameter :: bracket_growth = 1.6d0
  real(8), parameter :: min_step       = 1.d-8
  real(8), parameter :: max_step       = 1.d6
  integer, parameter :: max_bracket_iter = 200
  integer, parameter :: max_iter         = 150
  real(8), parameter :: zero_eps = 10.d0 * epsilon(1.d0)
  real(8) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
  real(8) :: left, right, f_left, f_right, step, fx
  integer :: iter
  logical :: bracketed, ok_left, ok_right, ok_new
  external :: f

  call evaluate(x_guess, fx, ok_new)
  if (abs(fx) <= tol) then
    return_value = x_guess
    return
  end if

  left  = x_guess
  right = x_guess
  f_left  = fx
  f_right = fx
  step = max(min_step, abs(x_guess)*0.1d0)
  bracketed = .false.

  do iter = 1, max_bracket_iter
    if (step > max_step) exit

    left = x_guess - step
    call evaluate(left, f_left, ok_left)
    if (.not. ok_left) then
      step = step * bracket_growth
      cycle
    end if
    if (abs(f_left) <= tol) then
      return_value = left
      return
    end if

    right = x_guess + step
    call evaluate(right, f_right, ok_right)
    if (.not. ok_right) then
      step = step * bracket_growth
      cycle
    end if
    if (abs(f_right) <= tol) then
      return_value = right
      return
    end if

    if (f_left * f_right <= 0.d0) then
      bracketed = .true.
      exit
    end if

    step = step * bracket_growth
  end do

  if (.not. bracketed) then
    stop "zbrent_diff: failed to bracket root"
  end if

  a = left
  b = right
  fa = f_left
  fb = f_right

  if (abs(fa) < abs(fb)) then
    call swap(a, b)
    call swap(fa, fb)
  end if

  c = a
  fc = fa
  d = b - a
  e = d

  do iter = 1, max_iter
    if (abs(fc) < abs(fb)) then
      call cycle_points(a, b, c, fa, fb, fc)
    end if

    tol1 = 2.d0 * zero_eps * abs(b) + 0.5d0 * tol
    xm = 0.5d0 * (c - b)

    if (abs(xm) <= tol1 .or. fb == 0.d0) then
      return_value = b
      return
    end if

    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s = fb / fa
      if (a == c) then
        p = 2.d0 * xm * s
        q = 1.d0 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * (2.d0 * xm * q * (q - r) - (b - a) * (r - 1.d0))
        q = (q - 1.d0) * (r - 1.d0) * (s - 1.d0)
      end if
      if (p > 0.d0) q = -q
      p = abs(p)

      if (2.d0 * p < min(3.d0 * xm * q - abs(tol1 * q), abs(e * q))) then
        e = d
        d = p / q
      else
        d = xm
        e = d
      end if
    else
      d = xm
      e = d
    end if

    a = b
    fa = fb

    if (abs(d) > tol1) then
      b = b + d
    else
      b = b + sign(tol1, xm)
    end if

    call evaluate(b, fb, ok_new)
    if (.not. ok_new) then
      b = b - d * 0.5d0
      call evaluate(b, fb, ok_new)
      if (.not. ok_new) stop "zbrent_diff: invalid function evaluation"
      e = c - b
      d = e
    end if
  end do

  stop "zbrent_diff: exceeded maximum iterations"

contains

  subroutine evaluate(x, fx, ok)
    real(8), intent(in) :: x
    real(8), intent(out) :: fx
    logical, intent(out) :: ok
    call f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p, dh)
    ok = .true.
    if (ieee_is_nan(fx)) ok = .false.
    if (ok) ok = ieee_is_finite(fx)
    if (.not. ok) fx = huge(1.d0)
  end subroutine evaluate

  subroutine swap(x, y)
    real(8), intent(inout) :: x, y
    real(8) :: tmp
    tmp = x
    x = y
    y = tmp
  end subroutine swap

  subroutine cycle_points(a_in, b_in, c_in, fa_in, fb_in, fc_in)
    real(8), intent(inout) :: a_in, b_in, c_in, fa_in, fb_in, fc_in
    call swap(a_in, b_in)
    call swap(fa_in, fb_in)
    call swap(a_in, c_in)
    call swap(fa_in, fc_in)
  end subroutine cycle_points

end subroutine zbrent_diff

subroutine zbrent_rot(x_guess, re, rho_p, ww_p, sgp, mugp, tol, return_value, f)
  ! Brent-style root finder with adaptive bracketing for the local rotation solver.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
  implicit none
  real(8), intent(in) :: x_guess, re, rho_p, ww_p, sgp, mugp, tol
  real(8), intent(out) :: return_value
  real(8), parameter :: bracket_growth = 1.6d0
  real(8), parameter :: min_step       = 1.d-8
  real(8), parameter :: max_step       = 1.d6
  integer, parameter :: max_bracket_iter = 200
  integer, parameter :: max_iter         = 150
  real(8), parameter :: zero_eps = 10.d0 * epsilon(1.d0)
  real(8) :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
  real(8) :: left, right, f_left, f_right, step, fx
  integer :: iter
  logical :: bracketed, ok_left, ok_right, ok_new
  external :: f

  call evaluate(x_guess, fx, ok_new)
  if (abs(fx) <= tol) then
    return_value = x_guess
    return
  end if

  left  = x_guess
  right = x_guess
  f_left  = fx
  f_right = fx
  step = max(min_step, abs(x_guess)*0.1d0)
  bracketed = .false.

  do iter = 1, max_bracket_iter
    if (step > max_step) exit

    left = x_guess - step
    call evaluate(left, f_left, ok_left)
    if (.not. ok_left) then
      step = step * bracket_growth
      cycle
    end if
    if (abs(f_left) <= tol) then
      return_value = left
      return
    end if

    right = x_guess + step
    call evaluate(right, f_right, ok_right)
    if (.not. ok_right) then
      step = step * bracket_growth
      cycle
    end if
    if (abs(f_right) <= tol) then
      return_value = right
      return
    end if

    if (f_left * f_right <= 0.d0) then
      bracketed = .true.
      exit
    end if

    step = step * bracket_growth
  end do

  if (.not. bracketed) then
    stop "zbrent_rot: failed to bracket root"
  end if

  a = left
  b = right
  fa = f_left
  fb = f_right

  if (abs(fa) < abs(fb)) then
    call swap(a, b)
    call swap(fa, fb)
  end if

  c = a
  fc = fa
  d = b - a
  e = d

  do iter = 1, max_iter
    if (abs(fc) < abs(fb)) then
      call cycle_points(a, b, c, fa, fb, fc)
    end if

    tol1 = 2.d0 * zero_eps * abs(b) + 0.5d0 * tol
    xm = 0.5d0 * (c - b)

    if (abs(xm) <= tol1 .or. fb == 0.d0) then
      return_value = b
      return
    end if

    if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
      s = fb / fa
      if (a == c) then
        p = 2.d0 * xm * s
        q = 1.d0 - s
      else
        q = fa / fc
        r = fb / fc
        p = s * (2.d0 * xm * q * (q - r) - (b - a) * (r - 1.d0))
        q = (q - 1.d0) * (r - 1.d0) * (s - 1.d0)
      end if
      if (p > 0.d0) q = -q
      p = abs(p)

      if (2.d0 * p < min(3.d0 * xm * q - abs(tol1 * q), abs(e * q))) then
        e = d
        d = p / q
      else
        d = xm
        e = d
      end if
    else
      d = xm
      e = d
    end if

    a = b
    fa = fb

    if (abs(d) > tol1) then
      b = b + d
    else
      b = b + sign(tol1, xm)
    end if

    call evaluate(b, fb, ok_new)
    if (.not. ok_new) then
      b = b - d * 0.5d0
      call evaluate(b, fb, ok_new)
      if (.not. ok_new) stop "zbrent_rot: invalid function evaluation"
      e = c - b
      d = e
    end if
  end do

  stop "zbrent_rot: exceeded maximum iterations"

contains

  subroutine evaluate(x, fx, ok)
    real(8), intent(in) :: x
    real(8), intent(out) :: fx
    logical, intent(out) :: ok
    call f(x, fx, re, rho_p, ww_p, sgp, mugp)
    ok = .true.
    if (ieee_is_nan(fx)) ok = .false.
    if (ok) ok = ieee_is_finite(fx)
    if (.not. ok) fx = huge(1.d0)
  end subroutine evaluate

  subroutine swap(x, y)
    real(8), intent(inout) :: x, y
    real(8) :: tmp
    tmp = x
    x = y
    y = tmp
  end subroutine swap

  subroutine cycle_points(a_in, b_in, c_in, fa_in, fb_in, fc_in)
    real(8), intent(inout) :: a_in, b_in, c_in, fa_in, fb_in, fc_in
    call swap(a_in, b_in)
    call swap(fa_in, fb_in)
    call swap(a_in, c_in)
    call swap(fa_in, fc_in)
  end subroutine cycle_points

end subroutine zbrent_rot

subroutine brent(ax, bx, cx, f, tol, xmin, ymin, ierr)
  !===============================================================================
  ! Brent's Method for Function Minimization
  !===============================================================================
  ! Finds the minimum of a function f using Brent's method, which combines
  ! golden section search with parabolic interpolation for improved efficiency.
  !
  ! Input:
  !   ax, bx, cx : Bracketing triplet (bx between ax and cx)
  !   f          : External function to minimize
  !   tol        : Fractional precision tolerance (e.g., 1.0d-8)
  !
  ! Output:
  !   xmin       : Abscissa of the minimum
  !   ymin       : Function value at minimum
  !   ierr       : Error flag (0=success, 1=max iterations, 2=invalid bracket)
  !===============================================================================
  implicit none
  
  ! Arguments
  real(8), intent(inout) :: ax, bx, cx
  real(8), intent(in)    :: tol
  real(8), intent(out)   :: xmin, ymin
  integer, intent(out)   :: ierr
  procedure(brent_func_1d) :: f  
  
  ! Parameters
  real(8), parameter :: GOLDEN_RATIO_COMPLEMENT = 0.3819660d0
  real(8), parameter :: ZEPS  = 1.0d-10      ! Small number for protection
  integer, parameter :: ITMAX = 100          ! Maximum iterations

  ! Interface for the function to be minimized
  abstract interface
    subroutine brent_func_1d(x, y)
      real(8), intent(in)  :: x
      real(8), intent(out) :: y
    end subroutine brent_func_1d
  end interface
  
  ! Local variables
  integer :: iter
  real(8) :: a, b, d, e, etemp, fu, fv, fw, fx
  real(8) :: p, q, r, tol1, tol2, u, v, w, x, xm
  real(8) :: fa, fb, fc
  logical :: accept_parabolic
  
  ierr = 0
  
  !-----------------------------------------------------------------------------
  ! Initialize and validate the bracketing triplet
  !-----------------------------------------------------------------------------
  call f(ax, fa)
  call f(bx, fb)
  call f(cx, fc)
  
  ! Initialize points: x is best, w is second best, v is previous w
  x = bx; w = bx; v = bx
  fx = fb; fw = fb; fv = fb

  if (fa < fx) then
    x = ax; fx = fa
  else if (fc < fx) then
    x = cx; fx = fc
  end if

  if (fb >= fa .or. fb >= fc) then
    ierr = 2
    xmin = bx
    ymin = fb
    return
  endif

  a = min(ax, cx)
  b = max(ax, cx)
  e = 0.0d0  ! Distance moved on step before last
  d = 0.0d0  ! Distance moved on last step
  
  !-----------------------------------------------------------------------------
  ! Main iteration loop
  !-----------------------------------------------------------------------------
  do iter = 1, ITMAX
    xm = 0.5d0 * (a + b)
    tol1 = tol * abs(x) + ZEPS
    tol2 = 2.0d0 * tol1
    
    ! Check convergence
    if (abs(x - xm) <= (tol2 - 0.5d0 * (b - a))) then
      xmin = x
      ymin = fx
      return
    endif
    
    !---------------------------------------------------------------------------
    ! Construct parabolic fit if possible
    !---------------------------------------------------------------------------
    accept_parabolic = .false.
    
    if (abs(e) > tol1) then
      ! Compute parabolic interpolation
      r = (x - w) * (fx - fv)
      q = (x - v) * (fx - fw)
      p = (x - v) * q - (x - w) * r
      q = 2.0d0 * (q - r)
      
      if (q > 0.0d0) then
        p = -p
      endif
      q = abs(q)
      etemp = e
      e = d
      
      ! Determine acceptability of parabolic fit
      if (abs(p) < abs(0.5d0 * q * etemp) .and. &
          p > q * (a - x) .and. &
          p < q * (b - x)) then
        ! Accept parabolic interpolation
        accept_parabolic = .true.
        d = p / q
        u = x + d
        
        ! f must not be evaluated too close to a or b
        if ((u - a) < tol2 .or. (b - u) < tol2) then
          d = sign(tol1, xm - x)
        endif
      endif
    endif
    
    !---------------------------------------------------------------------------
    ! Use golden section step if parabolic not acceptable
    !---------------------------------------------------------------------------
    if (.not. accept_parabolic) then
      if (x >= xm) then
        e = a - x
      else
        e = b - x
      endif
      d = GOLDEN_RATIO_COMPLEMENT * e
    endif
    
    !---------------------------------------------------------------------------
    ! Evaluate function at new point
    !---------------------------------------------------------------------------
    if (abs(d) >= tol1) then
      u = x + d
    else
      u = x + sign(tol1, d)
    endif
    
    call f(u, fu)
    
    !---------------------------------------------------------------------------
    ! Update a, b, v, w, x
    !---------------------------------------------------------------------------
    if (fu <= fx) then
      ! New minimum found
      if (u >= x) then
        a = x
      else
        b = x
      endif
      
      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu
    else
      ! Not a new minimum
      if (u < x) then
        a = u
      else
        b = u
      endif
      
      if (fu <= fw .or. w == x) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if (fu <= fv .or. v == x .or. v == w) then
        v = u
        fv = fu
      endif
    endif
  enddo
  
  !-----------------------------------------------------------------------------
  ! Maximum iterations reached
  !-----------------------------------------------------------------------------
  ierr = 1
  xmin = x
  ymin = fx
  
end subroutine brent