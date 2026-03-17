module brent_mod
  implicit none
  public :: find_omege_e, zbrent_rot

  abstract interface
    subroutine brent_func(x, fx)
      real(8), intent(in)  :: x
      real(8), intent(out) :: fx
    end subroutine brent_func
  end interface

contains

subroutine brent_core(x_guess, scale_up, scale_down, tol, return_value, f, small_guess, logfile)
! Shared Brent-style solver with adaptive bracketing.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  real(8), intent(in) :: x_guess, scale_up, scale_down, tol
  real(8), intent(out) :: return_value
  procedure(brent_func) :: f
  real(8), intent(in), optional :: small_guess
  character(*), intent(in), optional :: logfile

  integer :: iter, n, i_max, unit, ios
  logical :: bracket_found
  real(8), parameter :: ZEPS  = 1.d-8
  real(8) :: a,b,c,d,e, fa, fb, fc, p,q,r,s, xm, tol1
  real(8) :: ax, bx, cx, cx1, x0

  i_max = 200
  x0 = x_guess
  if (present(small_guess)) then
    if (abs(x0) < epsilon(x0)) x0 = small_guess
  end if

  ax = x0
  bx = x0
  bracket_found = .false.

  do iter = 1, i_max
    ax = ax * scale_up
    bx = bx / scale_down
    call f(ax, fa)
    call f(bx, fb)
    if (ieee_is_nan(fa)) then 
      ax = ax / scale_up
      call f(ax, fa)
    endif
    if (ieee_is_nan(fb)) then 
      bx = bx * scale_down
      call f(bx, fb)
    endif
    if ( fa*fb <= 0.d0 ) then
      exit
    endif
    if ( iter == i_max .and. present(logfile) ) then 
      open(newunit=unit,file=logfile, status='replace', action='write', iostat=ios)
      if (ios /= 0) then
        write(*,*) "write_eq_profile: failed to open file ", trim(logfile)
        return
      end if
      do n = 1, i_max
        ax = x_guess * 1.d-2 * n
        call f(ax, fa)
        cx = fa
        if ( n > 1 .and. cx * cx1 <= 0.d0 ) then
            bx = x_guess * 1.d-2 * n
            ax = x_guess * 1.d-2 * (n-1)
            bracket_found = .true.
            exit
        endif
        cx1 = cx
        write(unit,"(2es15.6)") ax, fa
      enddo
      close(unit)
      if (bracket_found) then
        call f(ax, fa)
        call f(bx, fb)
        exit  ! leave outer bracketing loop with new [ax, bx]
      else
        write(*,*) "check  ", logfile
        stop "brent_core: failed to bracket root"
      end if
    endif
  enddo

  a = ax
  b = bx

  c = b
  fc= fb
  e = 0.d0 
  d = 0.d0

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
    xm = (c-b) / 2.d0
    if ( abs(xm) <= tol1 .or. abs(fb) < epsilon(fb) ) then
      return_value = b
      return
    endif
    if ( abs(e) >= tol1 .and. abs(fa) > abs(fb) ) then ! Attempt inverse quadratic interpolation
      s = fb/fa
      if (abs(a - c) < epsilon(a)) then
        p = 2.d0 * xm * s
        q = 1.d0 - s
      else
        q = fa/fc
        r = fb/fc
        p = s * (2.d0 * xm * q * (q-r) - (b-a) * (r-1.d0))
        q = (q-1.d0) * (r-1.d0) * (s-1.d0)
      endif
      if ( p > 0.d0 ) q = -q ! Check whether in bounds
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

    call f(b, fb)
  enddo

  stop "brent_core: exceeding maximum iterations"

end subroutine brent_core

subroutine find_omege_e(x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol, return_value,f)
! Solve for Omega_e
  implicit none
  real(8), intent(in) :: x_guess, re, rho_h, g_h, w_h, rho_p, g_p, tol
  real(8), intent(out):: return_value
  interface
    subroutine f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
      real(8), intent(in)  :: x, re, rho_h, g_h, w_h, rho_p, g_p
      real(8), intent(out) :: fx
    end subroutine f
  end interface

  call brent_core(x_guess, 1.2d0, 1.1d0, tol, return_value, wrapped, logfile="./Cont/diff_rotation.dat")

  contains
    subroutine wrapped(x, fx)
      real(8), intent(in) :: x
      real(8), intent(out) :: fx
      call f(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    end subroutine wrapped

end subroutine find_omege_e

subroutine zbrent_rot(x_guess, re, rho_p, ww_p, sgp, mugp, tol, return_value, f)
! Point-wisely solve for Omega profile
  implicit none
  real(8), intent(in) :: x_guess, re, rho_p, ww_p, sgp, mugp, tol
  real(8), intent(out):: return_value
  interface
    subroutine f(x, fx, re, rho_p, ww_p, sgp, mugp)
      real(8), intent(in)  :: x, re, rho_p, ww_p, sgp, mugp
      real(8), intent(out) :: fx
    end subroutine f
  end interface

  call brent_core(x_guess, 1.1d0, 1.1d0, tol, return_value, wrapped, small_guess=1.d-4, logfile="./Cont/rotation_law.dat")

  contains
    subroutine wrapped(x, fx)
      real(8), intent(in) :: x
      real(8), intent(out):: fx
      call f(x, fx, re, rho_p, ww_p, sgp, mugp)
    end subroutine wrapped

end subroutine zbrent_rot

end module brent_mod
