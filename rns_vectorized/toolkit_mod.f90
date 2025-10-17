module toolkit_mod
  implicit none

contains
  elemental function interp_log_h_to_p(x) result(y)
    use para_mod
    implicit none
    real(8), intent(in) :: x
    real(8) :: y
    integer :: i

    if (x <= log_h(1)) then
      y = log_p(1)
    elseif (x >= log_h(num_tab)) then
      y = log_p(num_tab)
    else
      do i = 1, num_tab - 1
          if (x < log_h(i+1)) then
            y = log_p(i) + (log_p(i+1)-log_p(i)) * &
                  (x - log_h(i)) / (log_h(i+1)-log_h(i))
            return
          end if
      end do
    end if
  end function interp_log_h_to_p

  elemental function interp_log_p_to_e(x) result(y)
    use para_mod
    implicit none
    real(8), intent(in) :: x
    real(8) :: y
    integer :: i

    if (x <= log_p(1)) then
      y = log_e(1)
    elseif (x >= log_p(num_tab)) then
      y = log_e(num_tab)
    else
      do i = 1, num_tab - 1
          if (x < log_p(i+1)) then
            y = log_e(i) + (log_e(i+1)-log_e(i)) * &
                  (x - log_p(i)) / (log_p(i+1)-log_p(i))
            return
          end if
      end do
    end if
  end function interp_log_p_to_e

  subroutine interp(xp,yp,np, xb,yb)

    implicit none
    integer,intent(in) :: np
    integer :: n_nearest_pt, ii, kk, ir
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

  subroutine interp_pt(xp,yp,np, xb,yb)
    use para_mod, only: p_at_PT, C, KSCALE
    implicit none
    integer,intent(in) :: np
    integer :: n_nearest_pt, ii, kk, ir
    integer :: n_order = 4
    real(8),intent(in)  :: xp(np), yp(np)
    real(8),intent(in)  :: xb
    real(8),intent(out) :: yb
    real(8) :: fr 
    
    n_nearest_pt = minloc( abs(xb-xp), 1 )

    if ( abs(n_nearest_pt-p_at_PT) > n_order+1 ) then
      ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))
      yb = 0.d0
      do ii = -n_order, n_order
        fr = 1.d0
        do kk = -n_order, n_order
          if ( ii == kk ) cycle
          fr = fr * ( xb - xp(ir+kk) ) / ( xp(ir+ii) - xp(ir+kk) )
        enddo
        yb  = yb + fr * yp(ir+ii)
      enddo
    elseif ( n_nearest_pt == p_at_PT .or. n_nearest_pt == p_at_PT-1 ) then
      if (xp(p_at_PT) < xb ) then 
        yb = yp(p_at_PT+1) + ( xb - xp(p_at_PT+1) ) * ( yp(p_at_PT) - yp(p_at_PT+1) ) &
            / ( xp(p_at_PT) - xp(p_at_PT+1) )
      elseif ( xp(p_at_PT-1) < xb .and. xp(p_at_PT) > xb  ) then
        yb = yp(p_at_PT-1)
      else
        yb = yp(p_at_PT-2) + ( xb - xp(p_at_PT-2) ) * ( yp(p_at_PT-1) - yp(p_at_PT-2) ) &
          / ( xp(p_at_PT-1) - xp(p_at_PT-2) )
      endif
    else
      yb = merge( yp(n_nearest_pt-1) + ( xb - xp(n_nearest_pt-1) ) * ( yp(n_nearest_pt) - yp(n_nearest_pt-1) ) &
                          / ( xp(n_nearest_pt) - xp(n_nearest_pt-1) ), &
          yp(n_nearest_pt+1) + ( xb - xp(n_nearest_pt+1) ) * ( yp(n_nearest_pt) - yp(n_nearest_pt+1) ) &
                                      / ( xp(n_nearest_pt) - xp(n_nearest_pt+1) ), &
          xp(n_nearest_pt) > xb )
    endif

  end subroutine interp_pt

  !===========================================================
  ! First and second order derivatives in s and m directions
  ! using 2nd order finite differences
  !===========================================================
  real(8) function deriv_s(f, s, m)
    use para_mod, only : ds, SDIV, MDIV
    implicit none
    real(8), intent(in) :: f(SDIV, MDIV)
    integer, intent(in) :: s, m

    if (s == 1) then
      deriv_s = (f(2, m) - f(1, m)) / ds
    elseif (s == SDIV) then
      deriv_s = (f(SDIV, m) - f(SDIV-1, m)) / ds
    else
      deriv_s = (f(s+1, m) - f(s-1, m)) / (2.d0 * ds)
    end if
  end function deriv_s

  real(8) function deriv_ss(f, s, m)
    use para_mod, only : ds, SDIV, MDIV
    implicit none
    real(8), intent(in) :: f(SDIV, MDIV)
    integer, intent(in) :: s, m
    integer :: si

    si = max(4, min(s, SDIV-2))
    deriv_ss = (f(si+2, m) - 2.d0*f(si, m) + f(si-2, m)) / (4.d0 * ds**2)
  end function deriv_ss

  real(8) function deriv_m(f, s, m)
    use para_mod, only : dm, SDIV, MDIV
    implicit none
    real(8), intent(in) :: f(SDIV, MDIV)
    integer, intent(in) :: s, m

    if (m == 1) then
      deriv_m = (f(s, 2) - f(s, 1)) / dm
    elseif (m == MDIV) then
      deriv_m = (f(s, MDIV) - f(s, MDIV-1)) / dm
    else
      deriv_m = (f(s, m+1) - f(s, m-1)) / (2.d0 * dm)
    end if
  end function deriv_m

  real(8) function deriv_mm(f, s, m)
    use para_mod, only : dm, SDIV, MDIV
    implicit none
    real(8), intent(in) :: f(SDIV, MDIV)
    integer, intent(in) :: s, m
    integer :: mi

    mi = merge(2, m, m == 1)
    mi = merge(MDIV-1, mi, m == MDIV)
    deriv_mm = (f(s, mi+1) - 2.d0*f(s, mi) + f(s, mi-1)) / (dm**2)
  end function deriv_mm

  real(8) function deriv_sm(f, s, m)
    use para_mod, only : dm, ds, SDIV, MDIV
    implicit none
    real(8), intent(in) :: f(SDIV, MDIV)
    integer, intent(in) :: s, m

    select case (s)
    case (1)
      select case (m)
      case (1)
        deriv_sm = (f(2,2)-f(1,2)-f(2,1)+f(1,1))/(dm*ds)
      case (MDIV)
        deriv_sm = (f(2,MDIV)-f(1,MDIV)-f(2,MDIV-1)+f(1,MDIV-1))/(dm*ds)
      case default
        deriv_sm = (f(2,m+1)-f(2,m-1)-f(1,m+1)+f(1,m-1))/(2.d0*dm*ds)
      end select

    case (SDIV)
      select case (m)
      case (1)
        deriv_sm = (f(SDIV,2)-f(SDIV,1)-f(SDIV-1,2)+f(SDIV-1,1))/(dm*ds)
      case (MDIV)
        deriv_sm = (f(SDIV,MDIV)-f(SDIV-1,MDIV)-f(SDIV,MDIV-1)+f(SDIV-1,MDIV-1))/(dm*ds)
      case default
        deriv_sm = (f(SDIV,m+1)-f(SDIV,m-1)-f(SDIV-1,m+1)+f(SDIV-1,m-1))/(2.d0*dm*ds)
      end select

    case default
      select case (m)
      case (1)
        deriv_sm = (f(s+1,2)-f(s-1,2)-f(s+1,1)+f(s-1,1))/(2.d0*dm*ds)
      case (MDIV)
        deriv_sm = (f(s+1,MDIV)-f(s-1,MDIV)-f(s+1,MDIV-1)+f(s-1,MDIV-1))/(2.d0*dm*ds)
      case default
        deriv_sm = (f(s+1,m+1)-f(s-1,m+1)-f(s+1,m-1)+f(s-1,m-1))/(4.d0*dm*ds)
      end select
    end select
  end function deriv_sm
  
  real(8) function deriv_s_1d(f, s)
    use para_mod, only : ds, SDIV
    implicit none
    real(8), intent(in) :: f(SDIV)
    integer, intent(in) :: s

    if (s == 1) then
      deriv_s_1d = (f(2) - f(1)) / ds
    elseif (s == SDIV) then
      deriv_s_1d = (f(SDIV) - f(SDIV-1)) / ds
    else
      deriv_s_1d = (f(s+1) - f(s-1)) / (2.d0 * ds)
    end if
  end function deriv_s_1d

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
    real(8) :: fact, pmm, pmmp1, somx2, pll

    if(m<0 .or. m>l .or. abs(x)>1.d0) then
      write(*,*) m,l,x
      stop "Bad arguments in routine PLGNDR"
    endif

    pmm = 1.d0
    pll = 0.d0
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
        do ll = (m+2), l
          pll = (x * dble(2*ll-1) * pmmp1 - dble(ll+m-1)*pmm ) / dble(ll-m)
          pmm = pmmp1
          pmmp1 = pll
        enddo
        plgndr = pll
      endif
    endif

  end function plgndr

end module toolkit_mod