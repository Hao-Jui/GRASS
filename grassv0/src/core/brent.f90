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

subroutine brent(ax,bx,cx,f,tol,xmin,ymin)
  ! Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
  ! (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)), 
  ! this routine performs a golden section search for the minimum, isolating it to 
  ! a fractional precision of about tol. The abscissa of the minimum is returned as xmin, 
  ! and the minimum function value is returned as golden, the returned function value.
  
    implicit none
    real(8), intent(inout) :: ax,bx,cx,tol
    real(8), intent(out) :: xmin,ymin
    integer :: iter
    real(8), parameter :: CGOLD = 0.3819660d0
    real(8), parameter :: ZEPS  = 1.d-10
    real(8) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    real(8) :: fa_, fb_, fc_, dx
    external f
  
  
    call f(ax,fa_)
    call f(bx,fb_)
    call f(cx,fc_)
    
    if (fa_ < fb_ .and. fa_ < fc_) then 
      dx = ax
      ax = bx
      bx = dx
    elseif  (fc_ < fa_ .and. fc_ < fa_) then
      dx = cx
      cx = bx
      bx = dx
    endif
    a = min(ax,cx)
    b = max(ax,cx)
    v = bx
    w = v 
    x = v
    e = 0 ! This will be the distance moved on the step before last.
    d = 0.d0
  
    call f(x,fx)
    fv=fx
    fw=fx
  
    do iter = 1, 100
      xm = 0.5d0 * (a+b)
      tol1 = tol * abs(x)+ ZEPS
      tol2 = 2.d0*tol1 
      if( abs(x-xm) .le. (tol2-0.5d0*(b-a)) ) exit
      if( abs(e) > tol1 ) then
        r=(x-w)*(fx-fv) 
        q=(x-v)*(fx-fw) 
        p=(x-v)*q-(x-w)*r 
        q=2.*(q-r) 
        if( q > 0.d0 ) p=-p 
        q    = abs(q)
        etemp= e 
        e    = d
        if( abs(p).ge.abs(0.5d0*q*etemp).or.p.le.q*(a-x).or. p.ge.q*(b-x)) then 
          if(x.ge.xm) then
            e=a-x
          else 
            e=b-x
          endif
        else
          d    = p/q
          u    = x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d = sign(tol1,xm-x)
        endif
      endif
  
      if(abs(d).ge.tol1) then
        u=x+d
      else 
        u=x+sign(tol1,d)
      endif
  
      call  f(u, fu)
      if(fu.le.fx) then
        if(u.ge.x) then 
          a=x
        else 
          b=x
        endif
        v=w
        fv=fw
        w=x
        fw=fx
        x=u
        fx=fu
      else
        if(u.lt.x) then
          a=u 
        else
          b=u 
        endif
        if(fu.le.fw .or. w.eq.x) then 
          v=w
          fv=fw
          w=u
          fw=fu
        elseif(fu.le.fv .or. v.eq.x .or. v.eq.w) then 
          v=u
          fv=fu 
        endif
      endif
    enddo
  
    xmin=x
    ymin=fx
  
end subroutine brent  
