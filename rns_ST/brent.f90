subroutine zbrent_diff(x_guess, re, rho_e, g_e, w_e, rho_p, g_p, tol, return_value,f)
! Solve for Omega_e
  implicit none
  real(8), intent(in) :: x_guess, re, rho_e, g_e, w_e, rho_p, g_p, tol
  real(8), intent(out):: return_value
  integer :: iter, n
  real(8), parameter :: ZEPS  = 3.d-8
  real(8) :: a,b,c,d,e, fa, fb, fc, p,q,r,s, xm, tol1
  real(8) ::  ax, bx
  
  external f

  ax   = x_guess
  bx   = x_guess

  do iter = 1, 50
    ax = ax * 1.2d0
    bx = bx / 1.05d0
    call f(ax, fa, re, rho_e, g_e, w_e, rho_p, g_p)
    call f(bx, fb, re, rho_e, g_e, w_e, rho_p, g_p)
    if ( fa.ne.fa ) then 
      ax = ax / 1.2d0
      call f(ax, fa, re, rho_e, g_e, w_e, rho_p, g_p)
    endif
    if ( fb.ne.fb ) then 
      bx = bx * 1.05d0
      call f(bx, fb, re, rho_e, g_e, w_e, rho_p, g_p)
    endif
    if ( fa*fb <= 0.d0 ) then
      exit
    endif
    if ( iter == 50 ) then 
      open(43,file="./Cont/diff_rotation.dat")
      do n = 1, 100
        ax = x_guess * 8.d-3 * n
        call f(ax, fa, re, rho_e, g_e, w_e, rho_p, g_p)
        write(43,"(2es15.6)") 8.d-3 * n, fa
      enddo
      close(43)
      stop "root of Omega_e must be bracketed for zbrent, L32"
    endif
  enddo

  a = ax
  b = bx

  c = b
  fc= fb

  do iter = 1, 100
    if ( fb*fc > 0.d0) then ! Rename a, b, c and adjust bounding interval d
      c = a
      fc= fa
      d = b-a
      e = d
    endif
    if ( abs(fc) < abs(fb) ) then
      a = b
      b = c
      c = a
      fa= fb
      fb= fc
      fc= fa
    endif

    tol1 = 2.d0 * zeps *abs(b) + 0.5d0 * tol ! convergence check
    xm = (c-b)/2.d0
    if ( abs(xm) <= tol1 .or. fb==0.d0 ) then
      return_value = b
      return
    endif
    if ( abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then ! Attempt inverse quadratic interpolation
      s = fb/fa
      if (a==c) then
        p = 2.d0 * xm * s
        q = 1.d0 - s
      else
        q = fa/fc
        r = fb/fc
        p = s * (2.d0 * xm * q * (q-r) - (b-a) * (r-1.d0))
        q = (q-1.d0) * (r-1.d0) * (s-1.d0)
      endif
      if (p>0.d0) q = -q ! Check whether in bounds
      p = abs(p)
      if( 2.d0*p < min( 3.d0*xm*q - abs(tol1*q), abs(e*q) ) ) then
        e=d ! Accept interpolation
        d=p/q
      else
        d=xm ! Interpolation failed, use bisection
        e=d
      endif
    else ! Bounds decreasing too slowly, use bisection
      d = xm
      e = d
    endif
    a = b ! Move last best guess to a
    fa= fb
    if( abs(d) > tol1 ) then ! Evaluate new trial root
      b = b+d
    else
      b = b + sign(tol1,xm)
    endif

    call f(b, fb, re, rho_e, g_e, w_e, rho_p, g_p)
  enddo

  stop "zbrent exceeding maximum iterations, L191"

end subroutine zbrent_diff

subroutine zbrent_rot(x_guess, re, rho_p, ww_p, sgp, mugp, tol, return_value, f)
! Point-wisely solve for Omega profile
  implicit none
  real(8), intent(in) :: x_guess, re, rho_p, ww_p, sgp, mugp, tol
  real(8), intent(out):: return_value
  integer :: iter, n
  real(8), parameter :: ZEPS  = 3.d-8
  real(8) :: a,b,c,d,e, fa, fb, fc, p,q,r,s, xm, tol1
  real(8) ::  ax, bx
  
  external f

  ax   = x_guess
  bx   = x_guess

  do iter = 1, 40
    ax = ax * 1.1d0
    bx = bx / 1.1d0
    call f(ax, fa, re, rho_p, ww_p, sgp, mugp)
    call f(bx, fb, re, rho_p, ww_p, sgp, mugp)
    if ( fa.ne.fa ) then 
      ax = ax / 1.1d0
      call f(ax, fa, re, rho_p, ww_p, sgp, mugp)
    endif
    if ( fb.ne.fb ) then 
      bx = bx * 1.1d0
      call f(bx, fb, re, rho_p, ww_p, sgp, mugp)
    endif
    if ( fa*fb <= 0.d0 ) then
      exit
    endif
    if ( iter == 40 ) then 
      open(43,file="./Cont/rotation_law.dat")
      do n = 1, 100
        ax = x_guess * 8.d-3 * n
        call f(ax, fa, re, rho_p, ww_p, sgp, mugp)
        write(43,"(2es15.6)") 8.d-3 * n, fa
      enddo
      close(43)
      stop "root of omg(s,m) must be bracketed for zbrent, L122"
    endif
  enddo
  
  a = ax
  b = bx

  c = b
  fc= fb

  do iter = 1, 100
    if ( fb*fc > 0.d0) then ! Rename a, b, c and adjust bounding interval d
      c = a
      fc= fa
      d = b-a
      e = d
    endif
    if ( abs(fc) < abs(fb) ) then
      a = b
      b = c
      c = a
      fa= fb
      fb= fc
      fc= fa
    endif

    tol1 = 2.d0 * zeps *abs(b) + 0.5d0 * tol ! convergence check
    xm = (c-b)/2.d0
    if ( abs(xm) <= tol1 .or. fb==0.d0 ) then
      return_value = b; return
    endif
    if ( abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then ! Attempt inverse quadratic interpolation
      s = fb/fa
      if (a==c) then
        p = 2.d0 * xm * s
        q = 1.d0 - s
      else
        q = fa/fc
        r = fb/fc
        p = s * (2.d0 * xm * q * (q-r) - (b-a) * (r-1.d0))
        q = (q-1.d0) * (r-1.d0) * (s-1.d0)
      endif
      if (p>0.d0) q = -q ! Check whether in bounds
      p = abs(p)
      if( 2.d0*p < min( 3.d0*xm*q - abs(tol1*q), abs(e*q) ) ) then
        e=d ! Accept interpolation
        d=p/q
      else
        d=xm ! Interpolation failed, use bisection
        e=d
      endif
    else ! Bounds decreasing too slowly, use bisection
      d = xm
      e = d
    endif
    a = b ! Move last best guess to a
    fa= fb
    if( abs(d) > tol1 ) then ! Evaluate new trial root
      b = b+d
    else
      b = b + sign(tol1,xm)
    endif
    call f(b, fb, re, rho_p, ww_p, sgp, mugp)
    !write(*,"(2es18.9)") b,fb
  enddo

  stop "zbrent exceeding maximum iterations, L188"

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