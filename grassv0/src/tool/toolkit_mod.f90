module toolkit_mod
  implicit none

contains
  pure elemental function interp_log_h_to_p(x) result(y)
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

  pure elemental function interp_log_p_to_e(x) result(y)
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

  pure subroutine interp(xp,yp,np, xb,yb)

    implicit none
    integer,intent(in) :: np
    integer :: n_nearest_pt, ii, kk, ir
    integer, parameter :: n_order = 4
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

  pure subroutine interp_pt(xp,yp,np, xb,yb)
    use para_mod, only: p_at_PT, C, KSCALE
    implicit none
    integer,intent(in) :: np
    integer :: n_nearest_pt, ii, kk, ir
    integer, parameter :: n_order = 4
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

  pure subroutine interp_dual(xp, yp, np, xb, yb)
    use ad_mod
    implicit none
    integer, intent(in) :: np
    real(8), intent(in) :: xp(np), yp(np)
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb
    integer :: n_nearest_pt, ii, kk, ir
    integer, parameter :: n_order = 4
    type(dual) :: fr

    n_nearest_pt = minloc(abs(xb%val - xp),1)
    ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))

    yb = dual_const(0.d0)
    do ii = -n_order, n_order
      fr = dual_const(1.d0)
      do kk = -n_order, n_order
        if (ii == kk) cycle
        fr = fr * (xb - dual_const(xp(ir+kk))) / dual_const(xp(ir+ii) - xp(ir+kk))
      end do
      yb = yb + fr * dual_const(yp(ir+ii))
    end do
  end subroutine interp_dual

  ! **********************************************************************
  ! Integration helpers
  ! **********************************************************************
  subroutine integrate_profiles(x, values, results, err_estimates, ifails)
    use nag_compat_mod, only: d01gaf
    implicit none
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: values(:, :)
    real(8), intent(out) :: results(:)
    real(8), intent(out), optional :: err_estimates(:)
    integer, intent(out), optional :: ifails(:)
    integer :: n_points, n_cols, i
    real(8) :: err_local
    integer :: ifail_local

    n_points = size(x)
    if (size(values, 1) /= n_points) then
      write(*,*) "integrate_profiles: mismatched x and values dimensions"
      return
    end if

    n_cols = size(values, 2)
    if (size(results) < n_cols) then
      write(*,*) "integrate_profiles: results array too small"
      return
    end if
    if (present(err_estimates)) then
      if (size(err_estimates) < n_cols) then
        write(*,*) "integrate_profiles: err_estimates array too small"
        return
      end if
    end if
    if (present(ifails)) then
      if (size(ifails) < n_cols) then
        write(*,*) "integrate_profiles: ifails array too small"
        return
      end if
    end if

    do i = 1, n_cols
      err_local = 0.d0
      ifail_local = 0
      call d01gaf(x, values(:, i), n_points, results(i), err_local, ifail_local)
      if (present(err_estimates)) err_estimates(i) = err_local
      if (present(ifails)) ifails(i) = ifail_local
    end do
  end subroutine integrate_profiles

  subroutine integrate_column_spline(values, knots, result, status)
    use nag_compat_mod, only: e02baf, e02bbf
    use simpson_mod, only: simpson_1d
    real(8), intent(in)  :: values(:)
    real(8), intent(in)  :: knots(:)
    real(8), intent(out) :: result
    integer, intent(out) :: status
    integer :: n_points
    integer :: info_fit, info_int
    real(8), allocatable :: coeff_a(:), coeff_b(:), coeff_c(:), coeff_d(:)
    real(8) :: temp_mat(size(values),1), temp_vec(1)

    result = 0.d0
    status = 0
    n_points = size(values)
    if (n_points <= 1 .or. size(knots) /= n_points) then
      status = 1
      return
    end if

    allocate(coeff_a(n_points-1), coeff_b(n_points-1), coeff_c(n_points-1), coeff_d(n_points-1))
    info_fit = 0
    info_int = 0

    call e02baf(n_points, knots, values, coeff_a, coeff_b, coeff_c, coeff_d, info_fit)
    if (info_fit == 0) then
      call e02bbf(n_points, knots, coeff_a, coeff_b, coeff_c, coeff_d, knots(1), knots(n_points), result, info_int)
    end if

    if (info_fit /= 0 .or. info_int /= 0) then
      temp_mat(:,1) = values
      temp_vec = simpson_1d(temp_mat, knots(1), knots(n_points))
      result = temp_vec(1)
      status = 2
    end if
    deallocate(coeff_a, coeff_b, coeff_c, coeff_d)
  end subroutine integrate_column_spline

  pure function cumsum(x) result(y)
    real(8), intent(in) :: x(:)
    real(8) :: y(size(x))
    integer :: i
    y(1) = x(1)
    do i = 2, size(x)
      y(i) = y(i-1) + x(i)
    end do
  end function

  ! **********************************************************************
  ! First and second order derivatives in s and m directions
  ! using 2nd order finite differences
  ! **********************************************************************
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

    if (s == 1) then
      if (m == 1) then
        deriv_sm = (f(2,2)-f(1,2)-f(2,1)+f(1,1))/(dm*ds)
      elseif (m == MDIV) then
        deriv_sm = (f(2,MDIV)-f(1,MDIV)-f(2,MDIV-1)+f(1,MDIV-1))/(dm*ds)
      else
        deriv_sm = (f(2,m+1)-f(2,m-1)-f(1,m+1)+f(1,m-1))/(2.d0*dm*ds)
      end if
    elseif (s == SDIV) then
      if (m == 1) then
        deriv_sm = (f(SDIV,2)-f(SDIV,1)-f(SDIV-1,2)+f(SDIV-1,1))/(dm*ds)
      elseif (m == MDIV) then
        deriv_sm = (f(SDIV,MDIV)-f(SDIV-1,MDIV)-f(SDIV,MDIV-1)+f(SDIV-1,MDIV-1))/(dm*ds)
      else
        deriv_sm = (f(SDIV,m+1)-f(SDIV,m-1)-f(SDIV-1,m+1)+f(SDIV-1,m-1))/(2.d0*dm*ds)
      end if
    else
      if (m == 1) then
        deriv_sm = (f(s+1,2)-f(s-1,2)-f(s+1,1)+f(s-1,1))/(2.d0*dm*ds)
      elseif (m == MDIV) then
        deriv_sm = (f(s+1,MDIV)-f(s-1,MDIV)-f(s+1,MDIV-1)+f(s-1,MDIV-1))/(2.d0*dm*ds)
      else
        deriv_sm = (f(s+1,m+1)-f(s-1,m+1)-f(s+1,m-1)+f(s-1,m-1))/(4.d0*dm*ds)
      end if
    end if
  end function deriv_sm
  
  real(8) function deriv_s_1d(f, s)
    use para_mod, only : ds, SDIV
    implicit none
    real(8), intent(in) :: f(SDIV)
    integer, intent(in) :: s

    if (SDIV < 5) then
      ! Not enough points for the 4th order stencil; fall back to 2nd order.
      if (s == 1) then
        deriv_s_1d = (f(2) - f(1)) / ds
      elseif (s == SDIV) then
        deriv_s_1d = (f(SDIV) - f(SDIV-1)) / ds
      else
        deriv_s_1d = (f(s+1) - f(s-1)) / (2.d0 * ds)
      end if
      return
    end if

    select case (s)
    case (1)
      deriv_s_1d = (-25.d0*f(1) + 48.d0*f(2) - 36.d0*f(3) + 16.d0*f(4) - 3.d0*f(5)) / (12.d0 * ds)
    case (2)
      deriv_s_1d = (-3.d0*f(1) - 10.d0*f(2) + 18.d0*f(3) - 6.d0*f(4) + f(5)) / (12.d0 * ds)
    case default
      if (s == SDIV-1) then
        deriv_s_1d = (3.d0*f(SDIV) + 10.d0*f(SDIV-1) - 18.d0*f(SDIV-2) + 6.d0*f(SDIV-3) - f(SDIV-4)) / (12.d0 * ds)
      elseif (s == SDIV) then
        deriv_s_1d = (25.d0*f(SDIV) - 48.d0*f(SDIV-1) + 36.d0*f(SDIV-2) - 16.d0*f(SDIV-3) + 3.d0*f(SDIV-4)) / (12.d0 * ds)
      else
        deriv_s_1d = (-f(s+2) + 8.d0*f(s+1) - 8.d0*f(s-1) + f(s-2)) / (12.d0 * ds)
      end if
    end select
  end function deriv_s_1d

  ! **********************************************************************
  ! Special functions: Legendre polynomials, Associated ones, 
  ! Modified Spherical Bessel functions of two kinds
  ! **********************************************************************

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

  real(8) function besseli(n,x) 
    implicit none
    integer,intent(in) :: n
    real(8),intent(in) :: x
    integer :: ell_idx, Lrec, l_idx
    real(8) :: xx, i_prev, i_curr, i_next, two_ell_plus_one
    real(8) :: sinh_x, cosh_x, scale, besseli_pos
    real(8), parameter :: eps_small = 1.d-3
    real(8), parameter :: switch_downward = 20.d0
    real(8) :: i0_exact
    real(8), allocatable :: down_vals(:)
    real(8) :: sign_factor

    if (n < 0) then
      stop "besseli expects n >= 0"
    end if

    xx = abs(x)
    if (x < 0.d0 .and. mod(n,2) == 1) then
      sign_factor = -1.d0
    else
      sign_factor = 1.d0
    end if

    if (xx < eps_small) then
      besseli = clip_bessel(sign_factor * abs(spherical_i_series(n, xx)))
      return
    end if

    if (xx <= switch_downward) then
      Lrec = max(n + 40, 60)
      allocate(down_vals(0:Lrec+1))
      down_vals(Lrec+1) = 0.d0
      down_vals(Lrec)   = 1.d0
      do l_idx = Lrec, 1, -1
        down_vals(l_idx-1) = down_vals(l_idx+1) + ((2.d0*dble(l_idx)+1.d0)/xx) * down_vals(l_idx)
      end do
      sinh_x = sinh(xx)
      i0_exact = sinh_x / xx
      scale = i0_exact / down_vals(0)
      besseli_pos = abs(scale * down_vals(n))
      deallocate(down_vals)
      besseli = clip_bessel( sign_factor * besseli_pos )
      return
    end if

    if (xx > 700.d0) then
      sinh_x = 0.5d0 * exp(xx/2.d0) * exp(xx/2.d0)
      cosh_x = sinh_x
    elseif (xx < -700.d0) then
      sinh_x = -0.5d0 * exp(-xx/2.d0) * exp(-xx/2.d0)
      cosh_x = -sinh_x
    else
      sinh_x = sinh(xx)
      cosh_x = cosh(xx)
    end if

    i_prev = sinh_x / xx
    if (n == 0) then
      besseli = clip_bessel(sign_factor * abs(i_prev))
      return
    end if

    i_curr = (xx*cosh_x - sinh_x) / (xx*xx)
    if (n == 1) then
      besseli = clip_bessel(sign_factor * abs(i_curr))
      return
    end if

    do ell_idx = 1, n-1
      two_ell_plus_one = dble(2*ell_idx + 1)
      i_next = i_prev - (two_ell_plus_one/xx) * i_curr
      i_prev = clip_bessel(i_curr)
      i_curr = clip_bessel(i_next)
    end do

    besseli = clip_bessel(sign_factor * abs(i_curr))
  end function besseli

  real(8) function besselk(n,x)
    implicit none
    integer,intent(in) :: n
    real(8),intent(in) :: x
    integer :: ell
    real(8) :: xx, k_prev, k_curr, k_next, two_ell_plus_one
    real(8), parameter :: eps_small = 1.d-3
    real(8) :: exp_neg, em1

    if (n < 0) then
      stop "besselk expects n >= 0"
    end if

    xx = abs(x)
    if (xx == 0.d0) then
      besselk = 0.d0
      return
    end if

    if (xx < eps_small) then
      em1 = expm1_safe(-xx)
      exp_neg = 1.d0 + em1
    else
      exp_neg = exp(-xx)
    end if

    k_prev = clip_besselk(exp_neg / xx)
    if (n == 0) then
      besselk = k_prev
      return
    end if

    k_curr = clip_besselk(exp_neg * (1.d0 + 1.d0/xx) / xx)
    if (n == 1) then
      besselk = k_curr
      return
    end if

    do ell = 1, n-1
      two_ell_plus_one = dble(2*ell + 1)
      k_next = k_prev + (two_ell_plus_one/xx) * k_curr
      k_prev = k_curr
      k_curr = clip_besselk(k_next)
    end do

    besselk = k_curr
  end function besselk

  pure real(8) function odd_double_factorial(m)
    implicit none
    integer,intent(in) :: m
    integer :: k
    real(8) :: acc

    if (m <= 0) then
      odd_double_factorial = 1.d0
      return
    end if

    acc = 1.d0
    do k = m, 1, -2
      acc = acc * dble(k)
    end do

    odd_double_factorial = acc
  end function odd_double_factorial

  pure real(8) function spherical_i_series(n,x)
    implicit none
    integer,intent(in) :: n
    real(8),intent(in) :: x
    integer :: k
    real(8) :: term_factor, sum_series, x2, base

    if (n < 0) then
      spherical_i_series = 0.d0
      return
    end if

    x2 = x*x

    sum_series = 1.d0
    term_factor = 1.d0

    do k = 1, 64
      term_factor = term_factor * x2 / (2.d0*dble(k)*(2.d0*dble(n+k)+1.d0))
      sum_series = sum_series + term_factor
      if (abs(term_factor) < max(1.d-18, abs(sum_series)*1.d-16)) exit
    end do

    if (n == 0) then
      base = 1.d0
    else
      base = pow_int_real(x, n) / odd_double_factorial(2*n+1)
    end if

    spherical_i_series = clip_bessel(base * sum_series)
  end function spherical_i_series

  pure real(8) function pow_int_real(x,n)
    implicit none
    real(8),intent(in) :: x
    integer,intent(in) :: n
    integer :: k
    real(8) :: result

    result = 1.d0
    if (n <= 0) then
      pow_int_real = result
      return
    end if

    do k = 1, n
      result = result * x
    end do

    pow_int_real = result
  end function pow_int_real

  pure real(8) function expm1_safe(x)
    implicit none
    real(8),intent(in) :: x
    real(8) :: absx, term, sum
    integer :: k

    absx = abs(x)
    if (absx > 1.d-4) then
      expm1_safe = exp(x) - 1.d0
      return
    end if

    term = x
    sum = term
    do k = 2, 20
      term = term * x / dble(k)
      sum = sum + term
      if (abs(term) < 1.d-20) exit
    end do
    expm1_safe = sum
  end function expm1_safe

  pure real(8) function clip_bessel(val)
    implicit none
    real(8),intent(in) :: val
    real(8), parameter :: limit = 1.d98
    clip_bessel = min(limit, max(-limit, val))
  end function clip_bessel

  pure real(8) function clip_besselk(val)
    implicit none
    real(8),intent(in) :: val
    real(8) :: tmp
    tmp = abs(clip_bessel(val))
    if (tmp < tiny(1.d0)) tmp = tiny(1.d0)
    clip_besselk = tmp
  end function clip_besselk

  subroutine debug_mod_bessel
    use para_mod, only : LMAX, SDIV, s_gp
    implicit none
    integer :: n, s, unit_id
    real(8) :: xx, denom, val_i, val_k
    character(len=64) :: fname
    logical :: warned_large_x

    do n = 0, LMAX
      write(fname, "(A,I0,A)") "trig/spher_bessel_", 2*n, ".dat"
      open(newunit=unit_id, file=fname, status="replace", action="write")
      warned_large_x = .false.
      do s = 1, SDIV
        denom = 1.d0 - s_gp(s)
        if (abs(denom) < 1.d-12) denom = sign(1.d-12, denom)
        xx = s_gp(s) / denom
        if (xx > 1.d3) then
          if (.not. warned_large_x) then
            !write(*,"(A,1x,I0,1x,es15.6)") "debug_mod_bessel: skipping points with x >", n, xx
            warned_large_x = .true.
          end if
          val_k = 0.d0
          val_i = 0.d0
        else
          val_k = besselk(2*n, xx)
          val_i = besseli(2*n, xx)
        end if
        write(unit_id, "(4es25.16)") s_gp(s), xx, val_k, val_i
      end do
      close(unit_id)
    end do
  end subroutine debug_mod_bessel

end module toolkit_mod
