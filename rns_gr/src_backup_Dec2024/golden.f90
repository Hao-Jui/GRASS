subroutine golden(ax,bx,cx,f,tol,xmin,ymin)
! Given a function f, and given a bracketing triplet of abscissas ax, bx, cx 
! (such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)), 
! this routine performs a golden section search for the minimum, isolating it to 
! a fractional precision of about tol. The abscissa of the minimum is returned as xmin, 
! and the minimum function value is returned as golden, the returned function value.

  implicit none
  real(8), intent(inout) :: ax,bx,cx,tol
  real(8), intent(out) :: xmin,ymin
  real(8), parameter :: Rat = 0.61803399d0 
  real(8), parameter :: CR = 1.d0 - Rat
  real(8) :: f1, f2, f3, x0, x1, x2, x3
  real(8) :: dx, fa_, fb_, fc_

  external f


  call f(ax,fa_)
  call f(bx,fb_)
  call f(cx,fc_)
  
  ! need to ensure that ax > bx > cx, and f(bx) is the smallest
  if (fa_ < fb_ .and. fa_ < fc_) then 
    dx = ax
    ax = bx
    bx = dx
  elseif  (fc_ < fa_ .and. fc_ < fa_) then
    dx = cx
    cx = bx
    bx = dx
  endif
  if (ax < cx) then
    dx = ax
    ax = cx 
    cx = dx
  endif
  if (ax < bx) then 
    write(*,"(3es15.6)") ax, bx, cx
    write(*,"(3es15.6)") fa_, fb_, fc_
    stop "condition fails; L41 golden"
  endif
  x0 = ax ! At any given time we will keep track of four points, x0,x1,x2,x3
  x3 = cx 

  if ( abs(cx-bx) > abs(bx-ax) ) then ! Make x0 to x1 the smaller segment
    x1 = bx
    x2 = bx + CR*(cx-bx) ! and fill in the new point to be tried.
  else
    x2 = bx
    x1 = bx - CR*(bx-ax)
  endif

  ! The initial function evaluations. Note that we never need to
  ! evaluate the function at the original endpoints
  call f(x1,f1)
  call f(x2,f2)

  do while ( abs(x3-x0) > tol*( abs(x1)+abs(x2) ) )
    write(*,*) x2, f2
    if ( f2 < f1 ) then
      x0 = x1 
      x1 = x2
      x2 = Rat*x1 + CR*x3
      f1 = f2
      call f(x2,f2)
    else
      x3 = x2
      x2 = x1
      x1 = Rat*x2 + CR*x0
      f2 = f1
      call f(x1,f1)
    endif
  enddo

  if ( f1 < f2 ) then
    ymin = f1 
    xmin = x1
  else 
    ymin = f2
    xmin = x2
  endif

end subroutine golden

subroutine fun_hc(hc, fx)
  use para_mod
  implicit none
  real(8), intent(in) :: hc
  real(8), intent(out):: fx
  
  h_center = hc
  call spin
  call mass_radius
  
  fx = abs(mass_0/MSUN-Mb_goal)
  !write(*,*) mass_0/MSUN,fx
end subroutine fun_hc
  
subroutine fun_rerp(rerp, fx)
  use para_mod, only: J_goal,ang_mom,r_ratio
  implicit none
  real(8), intent(in) :: rerp
  real(8), intent(out):: fx
  
  r_ratio = rerp
  call spin
  call mass_radius
  
  fx = abs( ang_mom - J_goal )
  !write(*,*) ang_mom
end subroutine fun_rerp
  
