subroutine interp0(xp,yp,np, xb,yb)

  implicit none
  integer,intent(in) :: np
  integer :: n_nearest_pt,kk
  integer :: n_order = 4
  real(8),intent(in)  :: xp(np),yp(np)
  real(8),intent(in)  :: xb
  real(8),intent(out) :: yb
  real(8) :: fr
  
  
  n_nearest_pt = minloc(abs(xb-xp),np)
  
  kk = min( max( n_nearest_pt-(n_order-1)/2,1 ),np+1-n_order)
!  if( xb==xp(kk) .or.  xb==xp(kk+1) .or. xb==xp(kk+2) .or. xb==xp(kk+3)) &
!    xb = xb + DBL_EPSILON
  
  yb = (xb-xp(kk+1))*(xb-xp(kk+2))*(xb-xp(kk+3))*yp(kk) / &
        ((xp(kk)-xp(kk+1))*(xp(kk)-xp(kk+2))*(xp(kk)-xp(kk+3))) & 
    +(xb-xp(kk))*(xb-xp(kk+2))*(xb-xp(kk+3))*yp(kk+1)/ &
       ((xp(kk+1)-xp(kk))*(xp(kk+1)-xp(kk+2))*(xp(kk+1)-xp(kk+3))) & 
    +(xb-xp(kk))*(xb-xp(kk+1))*(xb-xp(kk+3))*yp(kk+2)/ &
       ((xp(kk+2)-xp(kk))*(xp(kk+2)-xp(kk+1))*(xp(kk+2)-xp(kk+3))) & 
    +(xb-xp(kk))*(xb-xp(kk+1))*(xb-xp(kk+2))*yp(kk+3)/ &
       ((xp(kk+3)-xp(kk))*(xp(kk+3)-xp(kk+1))*(xp(kk+3)-xp(kk+2)))

end subroutine interp0

subroutine interp(xp,yp,np, xb,yb)

  implicit none
  integer,intent(in) :: np
  integer :: n_nearest_pt,ii,kk,ir
  integer :: n_order = 4
  real(8),intent(in)  :: xp(np),yp(np)
  real(8),intent(in)  :: xb
  real(8),intent(out) :: yb
  real(8) :: fr 
  
  n_nearest_pt = minloc(abs(xb-xp),1)
  
  ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))

  yb = 0.d0
  do ii = -n_order,n_order
    fr = 1.d0
    do kk = -n_order,n_order
      if (ii==kk) cycle
      fr = fr * ( xb - xp(ir+kk) ) / ( xp(ir+ii) - xp(ir+kk) )
    enddo
    yb  = yb + fr * yp(ir+ii)
  enddo

end subroutine interp


real(8) function deriv_s(f,s,m)

  use para_mod, only : ds,SDIV,MDIV
  implicit none
  real(8) :: f(SDIV,MDIV)
  integer,intent(in):: s,m

  ! ------- v1
  if (s==1) then
    deriv_s = (f(s+1,m)-f(s,m))/ds
  elseif (s==SDIV) then
    deriv_s = (f(s,m)-f(s-1,m))/ds
  else
    deriv_s = (f(s+1,m)-f(s-1,m))/(dble(2)*ds)
  endif
 
  ! ------- v2
  !if (s < 3) then
    !	deriv_s = (-25.d0*f(s,m)+48.d0*f(s+1,m)-36.d0*f(s+2,m)+16.d0*f(s+3,m)-3.d0*f(s+4,m)) / (12.d0*ds)
  !elseif (s > SDIV-2) then
  !	deriv_s = (-25.d0*f(s,m)+48.d0*f(s-1,m)-36.d0*f(s-2,m)+16.d0*f(s-3,m)-3.d0*f(s-4,m)) / (12.d0*ds)
  !else
  !	deriv_s = (f(s+2,m)-8.d0*f(s+1,m)+8.d0*f(s-1,m)-f(s-2,m))/(12.d0*ds)
  !endif

end function deriv_s

real(8) function deriv_ss(f,s,m)

  use para_mod, only : SDIV,MDIV,ds
  implicit none
  real(8) :: f(SDIV,MDIV)
  integer :: s,m

  ! ------- v1
  if (s < 4) s = 4
  if (s > (SDIV-2)) s = SDIV-2
  deriv_ss = (f(s+2,m)-2.d0*f(s,m)+f(s-2,m))/(dble(4)*ds**2)
  
  ! ------- v2
  !if (s < 3) then 
  !	deriv_ss = (2.d0*f(s,m)-5.d0*f(s+1,m)+4.d0*f(s+2,m)-f(s+3,m)) / ds**2
  !elseif (s > SDIV-2) then
  !	deriv_ss = (2.d0*f(s,m)-5.d0*f(s-1,m)+4.d0*f(s-2,m)-f(s-3,m)) / ds**2
  !else
  !	deriv_ss = (-f(s+2,m)+16.d0*f(s+1,m)-30.d0*f(s,m)+16.d0*f(s-1,m)-f(s-2,m)) &
  !		/ (12.d0*ds**2)
  !endif

end function deriv_ss

real(8) function deriv_m(f,s,m)

  use para_mod, only : dm,SDIV,MDIV
  implicit none
  real(8) f(SDIV,MDIV)
  integer,intent(in) :: s,m

  ! ------- v1
  if (m==1) then
    deriv_m = (f(s,m+1)-f(s,m))/dm
  elseif (m==MDIV) then
    deriv_m = (f(s,m)-f(s,m-1))/dm
  else
    deriv_m = (f(s,m+1)-f(s,m-1))/(dble(2)*dm)
  endif

  ! ------- v2
  !if (m < 3) then
    !	deriv_m = (-25.d0*f(s,m)+48.d0*f(s,m+1)-36.d0*f(s,m+2)+16.d0*f(s,m+3)-3.d0*f(s,m+4)) / (12.d0*dm)
  !elseif (m > MDIV-2) then
  !	deriv_m = (-25.d0*f(s,m)+48.d0*f(s,m-1)-36.d0*f(s,m-2)+16.d0*f(s,m-3)-3.d0*f(s,m-4)) / (12.d0*dm)
  !else
  !	deriv_m = (-f(s,m+2)+8.d0*f(s,m+1)-8.d0*f(s,m-1)+f(s,m-2))/(12.d0*dm)
  !endif

end function deriv_m

real(8) function deriv_mm(f,s,m)

  use para_mod, only : SDIV,MDIV,dm
  implicit none
  real(8) :: f(SDIV,MDIV)
  integer :: s,m

  ! ------- v1
  if (m == 1) m = 2
  if (m == MDIV) m = MDIV-1
  deriv_mm = (f(s,m+1)-2.d0*f(s,m)+f(s,m-1))/dm**2
  
  ! ------- v2
  !if (m < 3) then 
  !	deriv_mm = (2.d0*f(s,m)-5.d0*f(s,m+1)+4.d0*f(s,m+2)-f(s,m+3)) / dm**2
  !elseif (m > MDIV-2) then
  !	deriv_mm = (2.d0*f(s,m)-5.d0*f(s,m-1)+4.d0*f(s,m-2)-f(s,m-3)) / dm**2
  !else
  !	deriv_mm = (-f(s,m+2)+16.d0*f(s,m+1)-30.d0*f(s,m)+16.d0*f(s,m-1)-f(s,m-2)) &
  !		/ (12.d0*dm**2)
  !endif

end function deriv_mm

real(8) function deriv_sm(f,s,m)

  use para_mod, only : MDIV,SDIV,dm,ds
  implicit none
  real(8) :: f(SDIV,MDIV)
  integer,intent(in):: s,m


  if (s==1) then
      if(m==1) then
        deriv_sm = (f(s+1,m+1)-f(s,m+1)-f(s+1,m)+f(s,m))/(dm*ds)
      else
        if (m==MDIV) then
          deriv_sm = (f(s+1,m)-f(s,m)-f(s+1,m-1)+f(s,m-1))/(dm*ds)
        else
          deriv_sm = (f(s+1,m+1)-f(s+1,m-1)-f(s,m+1)+f(s,m-1))/(dble(2)*dm*ds)
        endif
      endif

  else if (s==SDIV) then
    if(m.eq.1) then
      deriv_sm = (f(s,m+1)-f(s,m)-f(s-1,m+1)+f(s-1,m))/(dm*ds)
    else
      if (m.eq.MDIV) then
        deriv_sm = (f(s,m)-f(s-1,m)-f(s,m-1)+f(s-1,m-1))/(dm*ds)
      else
        deriv_sm = (f(s,m+1)-f(1,m-1)-f(s-1,m+1)+f(s-1,m-1))/(dble(2)*dm*ds)
      endif
    endif

  else
    if(m.eq.1) then
      deriv_sm = (f(s+1,m+1)-f(s-1,m+1)-f(s+1,m)+f(s-1,m))/(dble(2)*dm*ds)
    else
      if (m.eq.MDIV) then
        deriv_sm = (f(s+1,m)-f(s-1,m)-f(s+1,m-1)+f(s-1,m-1))/(dble(2)*dm*ds)
      else
        deriv_sm = (f(s+1,m+1)-f(s-1,m+1)-f(s+1,m-1)+f(s-1,m-1))/(dble(4)*dm*ds)
      endif
    endif
  endif

end function deriv_sm

real(8) function legendre(n,x)

  implicit none
  integer,intent(in) :: n
  real(8),intent(in) :: x
  integer :: i
  real(8) :: p,p_1,p_2

  p_2 = 1.d0
  p_1 = x

  if (n >= 2) then
    do i=2,n
      p = (x*(2.d0*dble(i)-1.d0)*p_1 - (dble(i)-1.d0)*p_2)/dble(i)
      p_2 = p_1
      p_1 = p
    enddo
    legendre = p
  else
    if (n == 1) then
      legendre = p_1
    else
      legendre = p_2
    endif
  endif

end function legendre

real(8) function plgndr(l,m,x)

  implicit none
  integer,intent(in) :: l,m
  real(8),intent(in) :: x
  integer :: ll
  real(8) :: fact,pmm,pmmp1,somx2,pll

  if(m<0 .or. m>l .or. abs(x)>1.d0) then
    write(*,*) m,l,x
    stop "Bad arguments in routine PLGNDR"
  endif

  pmm = 1.d0
  if ( m > 0 ) then
    somx2 = dsqrt((1.d0-x)*(1.d0+x))
    fact = 1.d0
    do ll=1,m
      pmm = pmm*(-fact*somx2)
      fact = fact + 2.d0
    enddo
  endif

  if (l == m) then
    plgndr = pmm
  else
    pmmp1 = x * dble(2*m+1) * pmm
    if(l==(m+1)) then
      plgndr = pmmp1
    else
      do ll=(m+2),l
        pll = (x * dble(2*ll-1) * pmmp1 - dble(ll+m-1)*pmm ) / dble(ll-m)
        pmm = pmmp1
        pmmp1 = pll
      enddo
      plgndr = pll
    endif
  endif

end function plgndr
