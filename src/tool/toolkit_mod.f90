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
  use precision_mod, only: wp
  use spectral_hub_mod, only: legendre_sequence
  implicit none

  integer, parameter :: n_order = 4
  integer, parameter :: pt_window = n_order + 1
  real(wp), parameter :: BARY_W(-n_order:n_order) = &
    [1.0_wp, -8.0_wp, 28.0_wp, -56.0_wp, 70.0_wp, -56.0_wp, 28.0_wp, -8.0_wp, 1.0_wp] / 40320.0_wp

  real(wp), allocatable, save :: bessel_down_workspace(:)
  integer, save :: bessel_down_size = 0

contains

  ! ---------------------------------------------------------------------------
  ! Pure utilities
  ! ---------------------------------------------------------------------------

  pure elemental logical function same_abscissa(xa, xb) result(is_same)
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
        is_same = .false.; return
      end if
    end do
  end function same_grid

  pure function binary_search_index(arr, n, x) result(idx)
    integer, intent(in) :: n
    real(wp), intent(in) :: arr(n), x
    integer :: idx, left, right, mid
    left = 1; right = n - 1
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

  pure real(wp) function pow_int_real(x, n)
    real(wp), intent(in) :: x
    integer, intent(in) :: n
    integer :: k
    pow_int_real = 1.0_wp
    do k = 1, n
      pow_int_real = pow_int_real * x
    end do
  end function pow_int_real

  pure real(wp) function expm1_safe(x)
    real(wp), intent(in) :: x
    real(wp) :: term, sum
    integer :: k
    if (abs(x) > 1.e-4_wp) then
      expm1_safe = exp(x) - 1.0_wp; return
    end if
    term = x; sum = term
    do k = 2, 20
      term = term * x / real(k, wp)
      sum = sum + term
      if (abs(term) < 1.e-20_wp) exit
    end do
    expm1_safe = sum
  end function expm1_safe

  pure real(wp) function clip_bessel(val)
    real(wp), intent(in) :: val
    real(wp), parameter :: LIMIT = 1.e98_wp
    clip_bessel = min(limit, max(-limit, val))
  end function clip_bessel

  pure real(wp) function clip_besselk(val)
    real(wp), intent(in) :: val
    real(wp) :: tmp
    tmp = abs(clip_bessel(val))
    if (tmp < tiny(1.0_wp)) tmp = tiny(1.0_wp)
    clip_besselk = tmp
  end function clip_besselk

  ! ---------------------------------------------------------------------------
  ! Interpolation — barycentric Lagrange with phase-transition awareness
  ! ---------------------------------------------------------------------------

  pure subroutine interp_linear_segment(xp, yp, i_left, i_right, xb, yb)
    integer, intent(in) :: i_left, i_right
    real(wp), intent(in) :: xp(:), yp(:), xb
    real(wp), intent(out) :: yb
    if (i_left >= i_right .or. same_abscissa(xp(i_left), xp(i_right))) then
      yb = yp(i_left)
    elseif (same_abscissa(xb, xp(i_left))) then
      yb = yp(i_left)
    elseif (same_abscissa(xb, xp(i_right))) then
      yb = yp(i_right)
    else
      yb = yp(i_left) + (xb - xp(i_left)) * (yp(i_right) - yp(i_left)) / (xp(i_right) - xp(i_left))
    end if
  end subroutine interp_linear_segment

  pure subroutine interp_linear_segment_dual(xp, yp, i_left, i_right, xb, yb)
    use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
    integer, intent(in) :: i_left, i_right
    real(wp), intent(in) :: xp(:), yp(:)
    type(dual), intent(in) :: xb
    type(dual), intent(out) :: yb
    if (i_left >= i_right .or. same_abscissa(xp(i_left), xp(i_right))) then
      yb = dual_const(yp(i_left))
    elseif (same_abscissa(xb%val, xp(i_left))) then
      yb = dual_const(yp(i_left))
    elseif (same_abscissa(xb%val, xp(i_right))) then
      yb = dual_const(yp(i_right))
    else
      yb = dual_const(yp(i_left)) + (xb - dual_const(xp(i_left))) &
         * dual_const((yp(i_right) - yp(i_left)) / (xp(i_right) - xp(i_left)))
    end if
  end subroutine interp_linear_segment_dual

  integer function nearest_transition_point(n_nearest_pt) result(pt)
    use para_mod, only: p_at_PT, n_PT
    integer, intent(in) :: n_nearest_pt
    integer :: j_closest
    pt = 0
    if (n_PT <= 0) return
    j_closest = minloc(abs(p_at_PT - n_nearest_pt), 1)
    if (abs(p_at_PT(j_closest) - n_nearest_pt) <= pt_window) pt = p_at_PT(j_closest)
  end function nearest_transition_point

  ! action=0: full barycentric  action=1: linear segment  action=2: constant
  subroutine pt_interp_action(xp, np, xb_val, n_nearest_pt, action, il, ir)
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

  pure subroutine interp(xp, yp, np, xb, yb)
    integer, intent(in)  :: np
    real(wp), intent(in)  :: xp(np), yp(np), xb
    real(wp), intent(out) :: yb
    integer :: n_nearest_pt, ir, ii
    real(wp) :: dx, wi, num, den, ds_uniform

    if (np >= 3) then
      ds_uniform = xp(2) - xp(1)
      if (abs((xp(3) - xp(2)) - ds_uniform) < epsilon(ds_uniform) * abs(ds_uniform)) then
        n_nearest_pt = min(np, max(1, nint(xb / ds_uniform) + 1))
      else
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
      if (abs(dx) < epsilon(dx)) then
        yb = yp(ir + ii); return
      end if
      wi  = bary_w(ii) / dx
      num = num + wi * yp(ir + ii)
      den = den + wi
    end do
    yb = num / den
  end subroutine interp

  pure subroutine interp_dual(xp, yp, np, xb, yb)
    use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
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
        yb = dual_const(yp(ir + ii)); return
      end if
      wi  = dual_const(bary_w(ii)) / dx
      num = num + wi * dual_const(yp(ir + ii))
      den = den + wi
    end do
    yb = num / den
  end subroutine interp_dual

  subroutine interp_pt(xp, yp, np, xb, yb)
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

  subroutine interp_pt_dual(xp, yp, np, xb, yb)
    use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
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

  ! ---------------------------------------------------------------------------
  ! Integration
  ! ---------------------------------------------------------------------------

  subroutine integrate_profiles(x, values, results, err_estimates, ifails)
    use nag_compat_mod, only: d01gaf
    use para_mod, only: mu, w_mu, angular_collocation, COLLOCATION_UNI
    real(wp), intent(in) :: x(:), values(:, :)
    real(wp), intent(out) :: results(:)
    real(wp), intent(out), optional :: err_estimates(:)
    integer, intent(out), optional :: ifails(:)
    integer :: n_points, n_cols, i
    real(wp) :: err_local
    integer :: ifail_local

    n_points = size(x)
    n_cols = size(values, 2)

    if (angular_collocation /= COLLOCATION_UNI .and. same_grid(x, mu)) then
      results(1:n_cols) = matmul(w_mu, values(:,1:n_cols))
      if (present(err_estimates)) err_estimates(1:n_cols) = 0.0_wp
      if (present(ifails)) ifails(1:n_cols) = 0
      return
    end if

    do i = 1, n_cols
      err_local = 0.0_wp; ifail_local = 0
      call d01gaf(x, values(:, i), n_points, results(i), err_local, ifail_local)
      if (present(err_estimates)) err_estimates(i) = err_local
      if (present(ifails)) ifails(i) = ifail_local
    end do
  end subroutine integrate_profiles

  ! ---------------------------------------------------------------------------
  ! 1D derivative (4th-order finite difference)
  ! ---------------------------------------------------------------------------

  real(wp) function deriv_s_1d(f, s)
    use para_mod, only : ds, SDIV
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

  ! ---------------------------------------------------------------------------
  ! Legendre polynomials
  ! ---------------------------------------------------------------------------

  real(wp) function legendre(n, x)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    real(wp) :: pnm1
    call legendre_sequence(n, x, legendre, pnm1)
  end function legendre

  real(wp) function plgndr(l, m, x)
    integer, intent(in) :: l, m
    real(wp), intent(in) :: x
    integer :: ll
    real(wp) :: fact, pmm, pmmp1, somx2, pll

    if (m < 0 .or. m > l .or. abs(x) > 1.0_wp) then
      write(*,*) m, l, x
      stop "Bad arguments in routine PLGNDR"
    endif

    pmm = 1.0_wp
    pll = 0.0_wp
    if (m > 0) then
      somx2 = dsqrt((1.0_wp - x) * (1.0_wp + x))
      fact = 1.0_wp
      do ll = 1, m
        pmm = pmm * (-fact * somx2)
        fact = fact + 2.0_wp
      enddo
    endif

    if (l == m) then
      plgndr = pmm
    else
      pmmp1 = x * real(2*m + 1, wp) * pmm
      if (l == (m + 1)) then
        plgndr = pmmp1
      else
        do ll = (m + 2), l
          pll = (x * real(2*ll - 1, wp) * pmmp1 - real(ll + m - 1, wp) * pmm) / real(ll - m, wp)
          pmm = pmmp1
          pmmp1 = pll
        enddo
        plgndr = pll
      endif
    endif
  end function plgndr

  ! ---------------------------------------------------------------------------
  ! Modified spherical Bessel functions
  ! ---------------------------------------------------------------------------

  subroutine ensure_bessel_workspace(needed_size)
    integer, intent(in) :: needed_size
    real(wp), allocatable :: tmp(:)
    if (.not. allocated(bessel_down_workspace)) then
      allocate(bessel_down_workspace(needed_size))
      bessel_down_size = needed_size
    else if (bessel_down_size < needed_size) then
      allocate(tmp(needed_size))
      deallocate(bessel_down_workspace)
      call move_alloc(tmp, bessel_down_workspace)
      bessel_down_size = needed_size
    end if
  end subroutine ensure_bessel_workspace

  pure real(wp) function odd_double_factorial(m)
    integer, intent(in) :: m
    integer :: k
    odd_double_factorial = 1.0_wp
    do k = m, 1, -2
      odd_double_factorial = odd_double_factorial * real(k, wp)
    end do
  end function odd_double_factorial

  pure real(wp) function spherical_i_series(n, x)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    integer :: k
    real(wp) :: term_factor, sum_series, x2, base

    if (n < 0) then
      spherical_i_series = 0.0_wp; return
    end if

    x2 = x*x
    sum_series = 1.0_wp
    term_factor = 1.0_wp
    do k = 1, 64
      term_factor = term_factor * x2 / (2.0_wp*real(k, wp)*(2.0_wp*real(n+k, wp)+1.0_wp))
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

  real(wp) function spherical_i_value(n, x)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    integer :: ell_idx, Lrec, l_idx
    real(wp) :: i_prev, i_curr, i_next, scale
    real(wp) :: sinh_x, cosh_x
    real(wp), parameter :: EPS_SMALL = 1.e-3_wp
    real(wp), parameter :: SWITCH_DOWNWARD = 20.0_wp

    if (x < eps_small) then
      spherical_i_value = spherical_i_series(n, x); return
    end if

    if (x <= switch_downward) then
      Lrec = max(n + 40, 60)
      call ensure_bessel_workspace(Lrec + 2)
      bessel_down_workspace(Lrec+1) = 0.0_wp
      bessel_down_workspace(Lrec)   = 1.0_wp
      do l_idx = Lrec, 1, -1
        bessel_down_workspace(l_idx-1) = bessel_down_workspace(l_idx+1) + &
          ((2.0_wp*real(l_idx, wp)+1.0_wp)/x) * bessel_down_workspace(l_idx)
      end do
      scale = (sinh(x) / x) / bessel_down_workspace(0)
      spherical_i_value = scale * bessel_down_workspace(n)
      return
    end if

    if (x > 700.0_wp) then
      sinh_x = 0.5_wp * exp(x); cosh_x = sinh_x
    else
      sinh_x = sinh(x); cosh_x = cosh(x)
    end if

    i_prev = sinh_x / x
    if (n == 0) then
      spherical_i_value = i_prev; return
    end if
    i_curr = (x*cosh_x - sinh_x) / (x*x)
    if (n == 1) then
      spherical_i_value = i_curr; return
    end if
    do ell_idx = 1, n-1
      i_next = i_prev - (real(2*ell_idx + 1, wp)/x) * i_curr
      i_prev = i_curr; i_curr = i_next
    end do
    spherical_i_value = i_curr
  end function spherical_i_value

  real(wp) function spherical_k_value(n, x)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    integer :: ell
    real(wp) :: k_prev, k_curr, k_next, exp_neg
    real(wp), parameter :: EPS_SMALL = 1.e-3_wp

    if (x < epsilon(x)) then
      spherical_k_value = 0.0_wp; return
    end if

    if (x < eps_small) then
      exp_neg = 1.0_wp + expm1_safe(-x)
    else
      exp_neg = exp(-x)
    end if

    k_prev = clip_besselk(exp_neg / x)
    if (n == 0) then
      spherical_k_value = k_prev; return
    end if
    k_curr = clip_besselk(exp_neg * (1.0_wp + 1.0_wp/x) / x)
    if (n == 1) then
      spherical_k_value = k_curr; return
    end if
    do ell = 1, n-1
      k_next = k_prev + (real(2*ell + 1, wp)/x) * k_curr
      k_prev = k_curr; k_curr = clip_besselk(k_next)
    end do
    spherical_k_value = k_curr
  end function spherical_k_value

  real(wp) function besseli(n, x)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    real(wp) :: xx, sign_factor
    if (n < 0) stop "besseli expects n >= 0"
    xx = abs(x)
    sign_factor = merge(-1.0_wp, 1.0_wp, x < 0.0_wp .and. mod(n, 2) == 1)
    besseli = clip_bessel(sign_factor * spherical_i_value(n, xx))
  end function besseli

  real(wp) function besselk(n, x)
    integer, intent(in) :: n
    real(wp), intent(in) :: x
    if (n < 0) stop "besselk expects n >= 0"
    besselk = spherical_k_value(n, abs(x))
  end function besselk

  subroutine bessel_even_tables(x, lmax, ivec, kvec)
    real(wp), intent(in) :: x
    integer, intent(in) :: lmax
    real(wp), intent(out) :: ivec(0:lmax), kvec(0:lmax)
    integer :: nord, ell, n, lrec, l_idx
    real(wp) :: xx, prev, curr, nxt, scale
    real(wp) :: sinh_x, cosh_x, exp_neg
    real(wp), parameter :: EPS_SMALL = 1.e-3_wp
    real(wp), parameter :: SWITCH_DOWNWARD = 20.0_wp

    if (lmax < 0) return
    xx = abs(x)
    nord = 2 * lmax

    ! --- i_ell ---
    if (xx < eps_small) then
      do n = 0, lmax
        ivec(n) = clip_bessel(spherical_i_series(2*n, xx))
      end do
    else if (xx <= switch_downward) then
      lrec = max(nord + 40, 60)
      call ensure_bessel_workspace(lrec + 2)
      bessel_down_workspace(lrec+1) = 0.0_wp
      bessel_down_workspace(lrec) = 1.0_wp
      do l_idx = lrec, 1, -1
        bessel_down_workspace(l_idx-1) = bessel_down_workspace(l_idx+1) + &
          ((2.0_wp * real(l_idx, wp) + 1.0_wp) / xx) * bessel_down_workspace(l_idx)
      end do
      scale = (sinh(xx) / xx) / bessel_down_workspace(0)
      do n = 0, lmax
        ivec(n) = clip_bessel(scale * bessel_down_workspace(2*n))
      end do
    else
      if (xx > 700.0_wp) then
        sinh_x = 0.5_wp * exp(xx); cosh_x = sinh_x
      else
        sinh_x = sinh(xx); cosh_x = cosh(xx)
      end if
      prev = sinh_x / xx
      ivec(0) = clip_bessel(prev)
      if (nord >= 1) curr = (xx * cosh_x - sinh_x) / (xx * xx)
      if (lmax >= 1) ivec(1) = clip_bessel(curr)
      do ell = 1, nord - 1
        nxt = prev - (2.0_wp * ell + 1.0_wp) / xx * curr
        prev = curr; curr = nxt
        if (mod(ell + 1, 2) == 0) ivec((ell + 1) / 2) = clip_bessel(curr)
      end do
    end if

    ! --- k_ell ---
    if (xx < epsilon(xx)) then
      kvec = 0.0_wp; return
    end if
    if (xx < eps_small) then
      exp_neg = 1.0_wp + expm1_safe(-xx)
    else
      exp_neg = exp(-xx)
    end if
    prev = clip_besselk(exp_neg / xx)
    kvec(0) = prev
    if (nord >= 1) curr = clip_besselk(exp_neg * (1.0_wp + 1.0_wp / xx) / xx)
    if (lmax >= 1) kvec(1) = curr
    do ell = 1, nord - 1
      nxt = prev + (2.0_wp * ell + 1.0_wp) / xx * curr
      prev = curr; curr = clip_besselk(nxt)
      if (mod(ell + 1, 2) == 0) kvec((ell + 1) / 2) = curr
    end do
  end subroutine bessel_even_tables

  ! ---------------------------------------------------------------------------
  ! Debug / diagnostics
  ! ---------------------------------------------------------------------------

  subroutine debug_mod_bessel
    use para_mod, only : LMAX, SDIV, s_gp
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
          if (.not. warned_large_x) warned_large_x = .true.
          val_k = 0.0_wp; val_i = 0.0_wp
        else
          val_k = besselk(2*n, xx); val_i = besseli(2*n, xx)
        end if
        write(unit_id, "(4es25.16)") s_gp(s), xx, val_k, val_i
      end do
      close(unit_id)
    end do
  end subroutine debug_mod_bessel

end module toolkit_mod
