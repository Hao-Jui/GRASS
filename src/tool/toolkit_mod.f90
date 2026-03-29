! toolkit_mod — numerical utilities for GRASS
!
! interp / interp_dual use barycentric Lagrange interpolation.
! For a uniform grid the stencil denominator products reduce to
!   prod_{k/=j}(x_j - x_k) = h^(2N) * (-1)^(j+N) * (N+j)! * (N-j)!
! where N = n_order and h is the grid spacing. 
! The h^8 cancels in the ratio num/den, leaving the signed binomial coefficients:
!   bary_w(j) = (-1)^(j+4) / ((4+j)! * (4-j)!)  =  C(8, j+4) * (-1)^(j+4) / 8! 
! whichi are independent of h and stencil position.
! Cost: O(2N+1) per evaluation vs O((2N+1)^2) for naive Lagrange.
module toolkit_mod
  use precision_mod, only: wp
  use spectral_hub, only: legendre_sequence
  implicit none

  ! Stencil order shared by interp, interp_pt, interp_dual
  integer, parameter :: n_order = 4 ! only this order is implemented
  integer, parameter :: pt_window = n_order + 1
  real(wp), parameter :: bary_w(-n_order:n_order) = &
      [1.0_wp, -8.0_wp, 28.0_wp, -56.0_wp, 70.0_wp, -56.0_wp, 28.0_wp, -8.0_wp, 1.0_wp] / 40320.0_wp

contains
  pure elemental logical function same_abscissa(xa, xb) result(is_same)
    implicit none
    real(wp), intent(in) :: xa, xb
    real(wp) :: scale

    scale = max(1.0_wp, abs(xa), abs(xb))
    is_same = abs(xa - xb) <= 16.0_wp * epsilon(scale) * scale
  end function same_abscissa

  pure logical function same_grid(x, grid) result(is_same)
    real(wp), intent(in) :: x(:), grid(:)
    integer :: i

    is_same = size(x) == size(grid)
    if (.not. is_same) return

    do i = 1, size(x)
      if (.not. same_abscissa(x(i), grid(i))) then
        is_same = .false.
        return
      end if
    end do
  end function same_grid

  integer function nearest_transition_point(n_nearest_pt) result(pt)
    use para_mod, only: p_at_PT, n_PT
    implicit none
    integer, intent(in) :: n_nearest_pt
    integer :: j_closest

    pt = 0
    if (n_PT <= 0) return

    j_closest = minloc(abs(p_at_PT - n_nearest_pt), 1)
    if (abs(p_at_PT(j_closest) - n_nearest_pt) <= pt_window) pt = p_at_PT(j_closest)
  end function nearest_transition_point

  ! PT decision logic shared by scalar and dual interp_pt variants.
  ! action=0: full barycentric (no PT nearby or invalid)
  ! action=1: linear segment between il and ir
  ! action=2: constant value at il (inside the PT gap)
  subroutine pt_interp_action(xp, np, xb_val, n_nearest_pt, action, il, ir)
    implicit none
    integer, intent(in) :: np, n_nearest_pt
    real(wp), intent(in) :: xp(np), xb_val
    integer, intent(out) :: action, il, ir
    integer :: pt

    pt = nearest_transition_point(n_nearest_pt)
    if (pt <= 0 .or. pt <= 1 .or. pt > np) then
      action = 0; return
    end if

    if (n_nearest_pt == pt .or. n_nearest_pt == pt - 1) then
      if (xb_val < xp(pt - 1)) then
        action = 1; il = max(1, pt - 2); ir = pt - 1
      elseif (xb_val > xp(pt)) then
        action = 1; il = pt; ir = min(np, pt + 1)
      else
        action = 2; il = pt - 1; ir = il
      end if
    else
      if (n_nearest_pt <= 1) then
        action = 2; il = 1; ir = 1
      elseif (n_nearest_pt >= np) then
        action = 2; il = np; ir = np
      elseif (xb_val < xp(n_nearest_pt)) then
        action = 1; il = n_nearest_pt - 1; ir = n_nearest_pt
      else
        action = 1; il = n_nearest_pt; ir = n_nearest_pt + 1
      end if
    end if
  end subroutine pt_interp_action

  pure subroutine interp_linear_segment(xp, yp, i_left, i_right, xb, yb)
    implicit none
    integer, intent(in) :: i_left, i_right
    real(wp), intent(in) :: xp(:), yp(:), xb
    real(wp), intent(out) :: yb

    if (i_left >= i_right) then
      yb = yp(i_left)
    elseif (same_abscissa(xb, xp(i_left))) then
      yb = yp(i_left)
    elseif (same_abscissa(xb, xp(i_right))) then
      yb = yp(i_right)
    elseif (same_abscissa(xp(i_left), xp(i_right))) then
      yb = yp(i_left)
    else
      yb = yp(i_left) + (xb - xp(i_left)) * (yp(i_right) - yp(i_left)) / (xp(i_right) - xp(i_left))
    end if
  end subroutine interp_linear_segment

  pure subroutine interp(xp, yp, np, xb, yb)
    implicit none
    integer, intent(in)  :: np
    real(wp), intent(in)  :: xp(np), yp(np), xb
    real(wp), intent(out) :: yb
    integer :: n_nearest_pt, ir, ii
    real(wp) :: dx, wi, num, den, ds_uniform

    ! O(1) nearest-index lookup for uniform grids: n_nearest_pt = nint(xb / ds) + 1
    ! Fallback to O(SDIV) minloc for non-uniform grids
    if (np >= 3) then
      ds_uniform = xp(2) - xp(1)
      if (abs((xp(3) - xp(2)) - ds_uniform) < epsilon(ds_uniform) * abs(ds_uniform)) then
        ! Uniform grid detected: use O(1) index calculation
        n_nearest_pt = min(np, max(1, nint(xb / ds_uniform) + 1))
      else
        ! Non-uniform grid: use O(np) minloc
        n_nearest_pt = minloc(abs(xb - xp), 1)
      end if
    else
      n_nearest_pt = minloc(abs(xb - xp), 1)
    end if

    ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))

    num = 0.0_wp
    den = 0.0_wp
    do ii = -n_order, n_order
      dx = xb - xp(ir + ii)
      if (abs(dx) < epsilon(dx)) then   ! xb lands exactly on a node
        yb = yp(ir + ii)
        return
      end if
      wi  = bary_w(ii) / dx
      num = num + wi * yp(ir + ii)
      den = den + wi
    end do
    yb = num / den
  end subroutine interp

  subroutine interp_pt(xp, yp, np, xb, yb)
    implicit none
    integer, intent(in) :: np
    real(wp), intent(in) :: xp(np), yp(np), xb
    real(wp), intent(out) :: yb
    integer :: n_nearest_pt, action, il, ir

    n_nearest_pt = minloc(abs(xb - xp), 1)
    if (same_abscissa(xb, xp(n_nearest_pt))) then
      yb = yp(n_nearest_pt); return
    end if

    call pt_interp_action(xp, np, xb, n_nearest_pt, action, il, ir)
    select case (action)
    case (1); call interp_linear_segment(xp, yp, il, ir, xb, yb)
    case (2); yb = yp(il)
    case default; call interp(xp, yp, np, xb, yb)
    end select
  end subroutine interp_pt

  pure subroutine interp_dual(xp, yp, np, xb, yb)
    use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
    implicit none
    integer, intent(in)     :: np
    real(wp), intent(in)     :: xp(np), yp(np)
    type(dual), intent(in)  :: xb
    type(dual), intent(out) :: yb
    integer    :: n_nearest_pt, ir, ii
    type(dual) :: dx, wi, num, den

    n_nearest_pt = minloc(abs(xb%val - xp), 1)
    ir = min(np - n_order, max(1 + n_order, n_nearest_pt - 1))

    num = dual_const(0.0_wp)
    den = dual_const(0.0_wp)
    do ii = -n_order, n_order
      dx = xb - dual_const(xp(ir + ii))
      if (abs(dx%val) < epsilon(dx%val)) then
        yb = dual_const(yp(ir + ii))
        return
      end if
      wi  = dual_const(bary_w(ii)) / dx
      num = num + wi * dual_const(yp(ir + ii))
      den = den + wi
    end do
    yb = num / den
  end subroutine interp_dual

  subroutine interp_pt_dual(xp, yp, np, xb, yb)
    use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
    implicit none
    integer, intent(in)    :: np
    real(wp), intent(in)   :: xp(np), yp(np)
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb
    integer :: n_nearest_pt, action, il, ir

    n_nearest_pt = minloc(abs(xb%val - xp), 1)
    if (same_abscissa(xb%val, xp(n_nearest_pt))) then
      yb = dual_const(yp(n_nearest_pt)); return
    end if

    call pt_interp_action(xp, np, xb%val, n_nearest_pt, action, il, ir)
    select case (action)
    case (1); call interp_linear_segment_dual(xp, yp, il, ir, xb, yb)
    case (2); yb = dual_const(yp(il))
    case default; call interp_dual(xp, yp, np, xb, yb)
    end select
  end subroutine interp_pt_dual

  pure subroutine interp_linear_segment_dual(xp, yp, i_left, i_right, xb, yb)
    use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
    implicit none
    integer, intent(in) :: i_left, i_right
    real(wp), intent(in) :: xp(:), yp(:)
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb

    if (i_left >= i_right) then
      yb = dual_const(yp(i_left))
    elseif (same_abscissa(xb%val, xp(i_left))) then
      yb = dual_const(yp(i_left))
    elseif (same_abscissa(xb%val, xp(i_right))) then
      yb = dual_const(yp(i_right))
    elseif (same_abscissa(xp(i_left), xp(i_right))) then
      yb = dual_const(yp(i_left))
    else
      yb = dual_const(yp(i_left)) + (xb - dual_const(xp(i_left))) &
         * dual_const((yp(i_right) - yp(i_left)) / (xp(i_right) - xp(i_left)))
    end if
  end subroutine interp_linear_segment_dual

  ! **********************************************************************
  ! Integration helpers
  ! **********************************************************************
  subroutine integrate_profiles(x, values, results, err_estimates, ifails)
    use nag_compat_mod, only: d01gaf
    use para_mod, only: mu, w_mu, angular_collocation, COLLOCATION_UNI
    implicit none
    real(wp), intent(in) :: x(:)
    real(wp), intent(in) :: values(:, :)
    real(wp), intent(out) :: results(:)
    real(wp), intent(out), optional :: err_estimates(:)
    integer, intent(out), optional :: ifails(:)
    integer :: n_points, n_cols, i
    real(wp) :: err_local
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

    if (angular_collocation /= COLLOCATION_UNI .and. same_grid(x, mu)) then
      results(1:n_cols) = matmul(w_mu, values(:,1:n_cols))
      if (present(err_estimates)) err_estimates(1:n_cols) = 0.0_wp
      if (present(ifails)) ifails(1:n_cols) = 0
      return
    end if

    do i = 1, n_cols
      err_local = 0.0_wp
      ifail_local = 0
      call d01gaf(x, values(:, i), n_points, results(i), err_local, ifail_local)
      if (present(err_estimates)) err_estimates(i) = err_local
      if (present(ifails)) ifails(i) = ifail_local
    end do
  end subroutine integrate_profiles

  ! **********************************************************************
  ! 1D derivative
  ! **********************************************************************
  real(wp) function deriv_s_1d(f, s)
    use para_mod, only : ds, SDIV
    implicit none
    real(wp), intent(in) :: f(SDIV)
    integer, intent(in) :: s
    real(wp) :: inv_ds

    if (SDIV < 5) then
      if (s == 1) then
        deriv_s_1d = (f(2) - f(1)) / ds
      elseif (s == SDIV) then
        deriv_s_1d = (f(SDIV) - f(SDIV-1)) / ds
      else
        deriv_s_1d = (f(s+1) - f(s-1)) / (2.0_wp * ds)
      end if
      return
    end if

    inv_ds = 1.0_wp / (12.0_wp * ds)

    if (s == 1) then
      deriv_s_1d = (-25.0_wp*f(1) + 48.0_wp*f(2) - 36.0_wp*f(3) + 16.0_wp*f(4) - 3.0_wp*f(5)) * inv_ds
    elseif (s == 2) then
      deriv_s_1d = (-3.0_wp*f(1) - 10.0_wp*f(2) + 18.0_wp*f(3) - 6.0_wp*f(4) + f(5)) * inv_ds
    elseif (s == SDIV-1) then
      deriv_s_1d = (3.0_wp*f(SDIV) + 10.0_wp*f(SDIV-1) - 18.0_wp*f(SDIV-2) + 6.0_wp*f(SDIV-3) - f(SDIV-4)) * inv_ds
    elseif (s == SDIV) then
      deriv_s_1d = (25.0_wp*f(SDIV) - 48.0_wp*f(SDIV-1) + 36.0_wp*f(SDIV-2) - 16.0_wp*f(SDIV-3) + 3.0_wp*f(SDIV-4)) * inv_ds
    else
      deriv_s_1d = (-f(s+2) + 8.0_wp*f(s+1) - 8.0_wp*f(s-1) + f(s-2)) * inv_ds
    end if
  end function deriv_s_1d

  ! **********************************************************************
  ! Special functions: Legendre polynomials, Associated ones, 
  ! Modified Spherical Bessel functions of two kinds
  ! **********************************************************************

  real(wp) function legendre(n,x)

    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: x
    real(wp) :: pnm1

    call legendre_sequence(n, x, legendre, pnm1)

  end function legendre

  ! associated Legendre polynomials P_l^m(x)
  real(wp) function plgndr(l,m,x)
    implicit none
    integer,intent(in) :: l,m
    real(wp),intent(in) :: x
    integer :: ll
    real(wp) :: fact, pmm, pmmp1, somx2, pll

    if(m<0 .or. m>l .or. abs(x)>1.0_wp) then
      write(*,*) m,l,x
      stop "Bad arguments in routine PLGNDR"
    endif

    pmm = 1.0_wp
    pll = 0.0_wp
    if ( m > 0 ) then
      somx2 = dsqrt((1.0_wp-x)*(1.0_wp+x))
      fact = 1.0_wp
      do ll=1,m
        pmm = pmm*(-fact*somx2)
        fact = fact + 2.0_wp
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

  real(wp) function besseli(n,x) 
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: x
    integer :: ell_idx, Lrec, l_idx
    real(wp) :: xx, i_prev, i_curr, i_next, two_ell_plus_one
    real(wp) :: sinh_x, cosh_x, scale, besseli_pos
    real(wp), parameter :: eps_small = 1.e-3_wp
    real(wp), parameter :: switch_downward = 20.0_wp
    real(wp) :: i0_exact
    real(wp), allocatable :: down_vals(:)
    real(wp) :: sign_factor

    if (n < 0) then
      stop "besseli expects n >= 0"
    end if

    xx = abs(x)
    if (x < 0.0_wp .and. mod(n,2) == 1) then
      sign_factor = -1.0_wp
    else
      sign_factor = 1.0_wp
    end if

    if (xx < eps_small) then
      besseli = clip_bessel(sign_factor * abs(spherical_i_series(n, xx)))
      return
    end if

    if (xx <= switch_downward) then
      Lrec = max(n + 40, 60)
      allocate(down_vals(0:Lrec+1))
      down_vals(Lrec+1) = 0.0_wp
      down_vals(Lrec)   = 1.0_wp
      do l_idx = Lrec, 1, -1
        down_vals(l_idx-1) = down_vals(l_idx+1) + ((2.0_wp*dble(l_idx)+1.0_wp)/xx) * down_vals(l_idx)
      end do
      sinh_x = sinh(xx)
      i0_exact = sinh_x / xx
      scale = i0_exact / down_vals(0)
      besseli_pos = abs(scale * down_vals(n))
      deallocate(down_vals)
      besseli = clip_bessel( sign_factor * besseli_pos )
      return
    end if

    if (xx > 700.0_wp) then
      sinh_x = 0.5e0_wp * exp(xx/2.0_wp) * exp(xx/2.0_wp)
      cosh_x = sinh_x
    elseif (xx < -700.0_wp) then
      sinh_x = -0.5e0_wp * exp(-xx/2.0_wp) * exp(-xx/2.0_wp)
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

  real(wp) function besselk(n,x)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: x
    integer :: ell
    real(wp) :: xx, k_prev, k_curr, k_next, two_ell_plus_one
    real(wp), parameter :: eps_small = 1.e-3_wp
    real(wp) :: exp_neg, em1

    if (n < 0) then
      stop "besselk expects n >= 0"
    end if

    xx = abs(x)
    if (abs(xx) < epsilon(xx)) then
      besselk = 0.0_wp
      return
    end if

    if (xx < eps_small) then
      em1 = expm1_safe(-xx)
      exp_neg = 1.0_wp + em1
    else
      exp_neg = exp(-xx)
    end if

    k_prev = clip_besselk(exp_neg / xx)
    if (n == 0) then
      besselk = k_prev
      return
    end if

    k_curr = clip_besselk(exp_neg * (1.0_wp + 1.0_wp/xx) / xx)
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

  subroutine bessel_even_tables(x, lmax, ivec, kvec)
    real(wp), intent(in) :: x
    integer, intent(in) :: lmax
    real(wp), intent(out) :: ivec(0:lmax), kvec(0:lmax)
    integer :: nord, ell, n, lrec
    real(wp) :: xx, sh, ch, exp_neg, prev, curr, nxt, scale
    real(wp), parameter :: eps_small = 1.e-3_wp, switch_downward = 20.0_wp
    real(wp), allocatable :: work(:)

    if (lmax < 0) return
    xx = abs(x)
    nord = 2 * lmax

    ! --- i_ell (modified spherical Bessel, 1st kind) ---
    if (xx < eps_small) then
      do n = 0, lmax
        ivec(n) = clip_bessel(abs(spherical_i_series(2*n, xx)))
      end do
    else if (xx <= switch_downward) then
      lrec = max(nord + 40, 60)
      allocate(work(0:lrec+1))
      work(lrec+1) = 0.0_wp; work(lrec) = 1.0_wp
      do ell = lrec, 1, -1
        work(ell-1) = work(ell+1) + (2.0_wp*ell + 1.0_wp)/xx * work(ell)
      end do
      scale = sinh(xx) / xx / work(0)
      do n = 0, lmax
        ivec(n) = clip_bessel(abs(scale * work(2*n)))
      end do
      deallocate(work)
    else
      if (xx > 700.0_wp) then
        sh = 0.5e0_wp * exp(xx/2.0_wp) * exp(xx/2.0_wp); ch = sh
      else
        sh = sinh(xx); ch = cosh(xx)
      end if
      prev = sh / xx
      curr = (xx*ch - sh) / (xx*xx)
      ivec(0) = clip_bessel(abs(prev))
      if (lmax >= 1) ivec(1) = clip_bessel(abs(curr))
      do ell = 1, nord - 1
        nxt = prev - (2.0_wp*ell + 1.0_wp)/xx * curr
        prev = clip_bessel(curr); curr = clip_bessel(nxt)
        if (mod(ell+1, 2) == 0) ivec((ell+1)/2) = clip_bessel(abs(curr))
      end do
    end if

    ! --- k_ell (modified spherical Bessel, 2nd kind) ---
    if (xx < epsilon(xx)) then
      kvec = 0.0_wp; return
    end if
    if (xx < eps_small) then
      exp_neg = 1.0_wp + expm1_safe(-xx)
    else
      exp_neg = exp(-xx)
    end if
    prev = clip_besselk(exp_neg / xx)
    curr = clip_besselk(exp_neg * (1.0_wp + 1.0_wp/xx) / xx)
    kvec(0) = prev
    if (lmax >= 1) kvec(1) = curr
    do ell = 1, nord - 1
      nxt = prev + (2.0_wp*ell + 1.0_wp)/xx * curr
      prev = curr; curr = clip_besselk(nxt)
      if (mod(ell+1, 2) == 0) kvec((ell+1)/2) = curr
    end do
  end subroutine bessel_even_tables

  pure real(wp) function odd_double_factorial(m)
    implicit none
    integer,intent(in) :: m
    integer :: k
    real(wp) :: acc

    if (m <= 0) then
      odd_double_factorial = 1.0_wp
      return
    end if

    acc = 1.0_wp
    do k = m, 1, -2
      acc = acc * dble(k)
    end do

    odd_double_factorial = acc
  end function odd_double_factorial

  pure real(wp) function spherical_i_series(n,x)
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: x
    integer :: k
    real(wp) :: term_factor, sum_series, x2, base

    if (n < 0) then
      spherical_i_series = 0.0_wp
      return
    end if

    x2 = x*x

    sum_series = 1.0_wp
    term_factor = 1.0_wp

    do k = 1, 64
      term_factor = term_factor * x2 / (2.0_wp*dble(k)*(2.0_wp*dble(n+k)+1.0_wp))
      sum_series = sum_series + term_factor
      if (abs(term_factor) < max(1.e-18_wp, abs(sum_series)*1.e-16_wp)) exit
    end do

    if (n == 0) then
      base = 1.0_wp
    else
      base = pow_int_real(x, n) / odd_double_factorial(2*n+1)
    end if

    spherical_i_series = clip_bessel(base * sum_series)
  end function spherical_i_series

  pure real(wp) function pow_int_real(x,n)
    implicit none
    real(wp),intent(in) :: x
    integer,intent(in) :: n
    integer :: k
    real(wp) :: result

    result = 1.0_wp
    if (n <= 0) then
      pow_int_real = result
      return
    end if

    do k = 1, n
      result = result * x
    end do

    pow_int_real = result
  end function pow_int_real

  pure real(wp) function expm1_safe(x)
    implicit none
    real(wp),intent(in) :: x
    real(wp) :: absx, term, sum
    integer :: k

    absx = abs(x)
    if (absx > 1.e-4_wp) then
      expm1_safe = exp(x) - 1.0_wp
      return
    end if

    term = x
    sum = term
    do k = 2, 20
      term = term * x / dble(k)
      sum = sum + term
      if (abs(term) < 1.e-20_wp) exit
    end do
    expm1_safe = sum
  end function expm1_safe

  pure real(wp) function clip_bessel(val)
    implicit none
    real(wp),intent(in) :: val
    real(wp), parameter :: limit = 1.e98_wp
    clip_bessel = min(limit, max(-limit, val))
  end function clip_bessel

  pure real(wp) function clip_besselk(val)
    implicit none
    real(wp),intent(in) :: val
    real(wp) :: tmp
    tmp = abs(clip_bessel(val))
    if (tmp < tiny(1.0_wp)) tmp = tiny(1.0_wp)
    clip_besselk = tmp
  end function clip_besselk

  pure function binary_search_index(arr, n, x) result(idx)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: arr(n), x
    integer :: idx, left, right, mid

    left = 1
    right = n - 1
    do while (left < right)
      mid = (left + right) / 2
      if (arr(mid) < x) then
        left = mid + 1
      else
        right = mid
      end if
    end do
    idx = left
  end function binary_search_index

  subroutine debug_mod_bessel
    use para_mod, only : LMAX, SDIV, s_gp
    implicit none
    integer :: n, s, unit_id
    real(wp) :: xx, denom, val_i, val_k
    character(len=64) :: fname
    logical :: warned_large_x

    do n = 0, LMAX
      write(fname, "(A,I0,A)") "trig/spher_bessel_", 2*n, ".dat"
      open(newunit=unit_id, file=fname, status="replace", action="write")
      warned_large_x = .false.
      do s = 1, SDIV
        denom = 1.0_wp - s_gp(s)
        if (abs(denom) < 1.e-12_wp) denom = sign(1.e-12_wp, denom)
        xx = s_gp(s) / denom
        if (xx > 1.e3_wp) then
          if (.not. warned_large_x) then
            !write(*,"(A,1x,I0,1x,es15.6)") "debug_mod_bessel: skipping points with x >", n, xx
            warned_large_x = .true.
          end if
          val_k = 0.0_wp
          val_i = 0.0_wp
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
