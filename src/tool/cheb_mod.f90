module cheb_mod
  use precision_mod, only: wp
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr
  implicit none
  public

  type :: cheb_fit_stats
    integer :: n_deg = 0
    integer :: n_pts = 0
    integer :: tail_modes = 0
    integer :: solve_info = 0
    logical :: used_dct = .false.
    real(wp) :: rms_abs = 0.0_wp
    real(wp) :: rms_rel = 0.0_wp
    real(wp) :: linf_abs = 0.0_wp
    real(wp) :: linf_rel = 0.0_wp
    real(wp) :: coeff_max_abs = 0.0_wp
    real(wp) :: coeff_tail_max_abs = 0.0_wp
    real(wp) :: coeff_tail_l2_ratio = 0.0_wp
    real(wp) :: coeff_tail_max_ratio = 0.0_wp
  end type cheb_fit_stats
  
  real(wp), parameter, private :: PI = acos(-1.e0_wp)
  integer(c_int), parameter, private :: fftw_redft00 = 3_c_int
  interface
    function fftw_plan_r2r_1d(n, in, out, kind, flags) bind(C, name="fftw_plan_r2r_1d") result(plan)
      import :: c_int, c_double, c_ptr
      integer(c_int), value :: n
      real(c_double) :: in(*), out(*)
      integer(c_int), value :: kind, flags
      type(c_ptr) :: plan
    end function fftw_plan_r2r_1d
    subroutine fftw_execute_r2r(plan, in, out) bind(C, name="fftw_execute_r2r")
      import :: c_ptr, c_double
      type(c_ptr), value :: plan
      real(c_double) :: in(*), out(*)
    end subroutine fftw_execute_r2r
    subroutine fftw_destroy_plan(plan) bind(C, name="fftw_destroy_plan")
      import :: c_ptr
      type(c_ptr), value :: plan
    end subroutine fftw_destroy_plan
  end interface
contains

  ! Fit Chebyshev expansion coefficients (degree n_deg) from n_pts arbitrary
  ! nodes in [a,b].  If the nodes are Chebyshev-Lobatto, use FFTW DCT-I.
  ! Otherwise fall back to a least-squares fit.
  subroutine cheb_std_base(n_deg, n_pts, x_nodes, h_nodes, a, b, coeffs, stats)
    integer, intent(in)  :: n_deg, n_pts
    real(wp), intent(in)  :: x_nodes(n_pts), h_nodes(n_pts), a, b
    real(wp), intent(out) :: coeffs(0:n_deg)
    type(cheb_fit_stats), intent(out), optional :: stats
    real(wp), allocatable :: xi(:), T(:,:), rhs(:,:), work(:)
    integer :: k, info, lwork
    interface
      subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
        import :: wp
        character(len=1), intent(in)    :: trans
        integer, intent(in)             :: m, n, nrhs, lda, ldb, lwork
        integer, intent(out)            :: info
        real(wp), intent(inout)         :: a(lda,*), b(ldb,*)
        real(wp), intent(inout)         :: work(*)
      end subroutine dgels
    end interface

    if (nodes_are_cheb_lobatto(n_deg, n_pts, x_nodes, a, b)) then
      call cheb_coeffs_dct1(n_deg, h_nodes, coeffs)
      if (present(stats)) call fill_cheb_fit_stats(n_deg, n_pts, x_nodes, h_nodes, a, b, coeffs, 0, .true., stats)
      return
    end if

    allocate(xi(n_pts), T(n_pts, 0:n_deg), rhs(n_pts, 1))

    xi = 2.e0_wp*(x_nodes - a)/(b - a) - 1.e0_wp

    T(:,0) = 1.e0_wp
    if (n_deg >= 1) T(:,1) = xi
    do k = 2, n_deg
      T(:,k) = 2.e0_wp*xi*T(:,k-1) - T(:,k-2)
    end do

    ! Workspace query
    lwork = -1
    allocate(work(1))
    rhs(:,1) = h_nodes
    call dgels('N', n_pts, n_deg+1, 1, T, n_pts, rhs, n_pts, work, lwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    rhs(:,1) = h_nodes
    call dgels('N', n_pts, n_deg+1, 1, T, n_pts, rhs, n_pts, work, lwork, info)
    if (info /= 0) then
      write(*,*) 'cheb_std_base: dgels failed, info =', info
      coeffs = 0.e0_wp
    else
      coeffs = rhs(1:n_deg+1, 1)
    end if
    if (present(stats)) call fill_cheb_fit_stats(n_deg, n_pts, x_nodes, h_nodes, a, b, coeffs, info, .false., stats)
    deallocate(xi, T, rhs, work)
  end subroutine cheb_std_base

  logical pure function nodes_are_cheb_lobatto(n_deg, n_pts, x_nodes, a, b)
    integer, intent(in) :: n_deg, n_pts
    real(wp), intent(in) :: x_nodes(n_pts), a, b
    real(wp) :: x_expected, half_width, center, tol
    integer :: j

    nodes_are_cheb_lobatto = .false.
    if (n_pts /= n_deg + 1) return
    if (n_deg < 1) then
      nodes_are_cheb_lobatto = (n_pts == 1)
      return
    end if

    half_width = 0.5e0_wp * (b - a)
    center = 0.5e0_wp * (a + b)
    tol = 100.e0_wp * epsilon(1.e0_wp) * max(1.e0_wp, abs(a), abs(b))

    do j = 0, n_deg
      x_expected = center + half_width * cos(pi * real(j, wp) / real(n_deg, wp))
      if (abs(x_nodes(j+1) - x_expected) > tol) return
    end do
    nodes_are_cheb_lobatto = .true.
  end function nodes_are_cheb_lobatto

  subroutine cheb_coeffs_dct1(n_deg, h_nodes, coeffs)
    integer, intent(in) :: n_deg
    real(wp), intent(in) :: h_nodes(n_deg+1)
    real(wp), intent(out) :: coeffs(0:n_deg)
    real(wp), allocatable :: work_in(:), work_out(:)
    type(c_ptr) :: plan
    integer(c_int), parameter :: fftw_estimate = 64_c_int

    if (n_deg == 0) then
      coeffs(0) = h_nodes(1)
      return
    end if

    allocate(work_in(n_deg+1), work_out(n_deg+1))
    work_in = h_nodes
    plan = fftw_plan_r2r_1d(int(n_deg + 1, c_int), work_in, work_out, fftw_redft00, fftw_estimate)
    call fftw_execute_r2r(plan, work_in, work_out)
    call fftw_destroy_plan(plan)

    coeffs = work_out / real(n_deg, wp)
    coeffs(0) = 0.5e0_wp * coeffs(0)
    coeffs(n_deg) = 0.5e0_wp * coeffs(n_deg)

    deallocate(work_in, work_out)
  end subroutine cheb_coeffs_dct1

  subroutine cheb_diff_matrix(N, D)
  ! Trefethen’s
    integer, intent(in) :: N
    real(wp), intent(out) :: D(0:N,0:N)
    real(wp) :: x(0:N), c(0:N), factor
    integer :: i, j

    do i = 0, N
      x(i) = cos( pi * dble(i) / dble(N) )
      if (i == 0 .or. i == N) then
        c(i) = 2.e0_wp
      else
        c(i) = 1.e0_wp
      end if
    end do

    do i = 0, N
      do j = 0, N
        if (i == j) then
          if (i == 0) then
            D(i,j) = ( 2.e0_wp * N**2 + 1.e0_wp ) / 6.e0_wp
          else if (i == N) then
            D(i,j) = -( 2.e0_wp * N**2 + 1.e0_wp ) / 6.e0_wp
          else
            D(i,j) = -x(i) / ( 2.e0_wp * ( 1.e0_wp - x(i)**2 ) )
          end if
        else
          factor = (c(i)/c(j)) * (-1.e0_wp)**(i+j)
          D(i,j) = factor / (x(i) - x(j))
        end if
      end do
    end do
  end subroutine cheb_diff_matrix

  ! Evaluate Chebyshev expansion at a single point x in [a,b].
  ! Uses Clenshaw's algorithm: O(n), numerically stable.
  pure function cheb_get_val_point(n, coeffs, a, b, x) result(val)
    integer,  intent(in) :: n
    real(wp), intent(in) :: coeffs(0:n), a, b, x
    real(wp) :: val
    real(wp) :: xi, bk, bk1, bk2
    integer  :: k

    xi  = 2.e0_wp*(x - a)/(b - a) - 1.e0_wp
    bk1 = 0.e0_wp
    bk2 = 0.e0_wp
    do k = n, 1, -1
      bk  = coeffs(k) + 2.e0_wp*xi*bk1 - bk2
      bk2 = bk1
      bk1 = bk
    end do
    val = coeffs(0) + xi*bk1 - bk2
  end function cheb_get_val_point

  function cheb_get_deriv_point(n, coeffs, a, b, x, order) result(val)
    integer, intent(in) :: n, order
    real(wp), intent(in) :: coeffs(0:n), a, b, x
    real(wp) :: val
    real(wp), allocatable :: cur(:), next(:)
    real(wp) :: scale
    integer :: deg, m

    if (order < 0) error stop "cheb_get_deriv_point: derivative order must be non-negative"
    if (order == 0) then
      val = cheb_get_val_point(n, coeffs, a, b, x)
      return
    end if
    if (order > n) then
      val = 0.0_wp
      return
    end if

    allocate(cur(0:n))
    cur = coeffs
    deg = n
    scale = 1.0_wp
    do m = 1, order
      allocate(next(0:deg-1))
      call cheb_diff_coeffs(deg, cur, next)
      scale = scale * (2.0_wp / (b - a))
      deallocate(cur)
      allocate(cur(0:deg-1))
      cur = next
      deallocate(next)
      deg = deg - 1
    end do

    val = scale * cheb_get_val_point(deg, cur, a, b, x)
    deallocate(cur)
  end function cheb_get_deriv_point

  pure subroutine cheb_diff_coeffs(n, coeffs, deriv_coeffs)
    integer, intent(in) :: n
    real(wp), intent(in) :: coeffs(0:n)
    real(wp), intent(out) :: deriv_coeffs(0:n-1)
    integer :: k

    deriv_coeffs = 0.0_wp
    deriv_coeffs(n-1) = 2.0_wp * real(n, wp) * coeffs(n)
    if (n > 1) deriv_coeffs(n-2) = 2.0_wp * real(n - 1, wp) * coeffs(n-1)
    do k = n - 3, 0, -1
      deriv_coeffs(k) = deriv_coeffs(k+2) + 2.0_wp * real(k + 1, wp) * coeffs(k+1)
    end do
    deriv_coeffs(0) = 0.5_wp * deriv_coeffs(0)
  end subroutine cheb_diff_coeffs

  subroutine fill_cheb_fit_stats(n_deg, n_pts, x_nodes, h_nodes, a, b, coeffs, solve_info, used_dct, stats)
    integer, intent(in) :: n_deg, n_pts, solve_info
    real(wp), intent(in) :: x_nodes(n_pts), h_nodes(n_pts), a, b, coeffs(0:n_deg)
    logical, intent(in) :: used_dct
    type(cheb_fit_stats), intent(out) :: stats
    real(wp), allocatable :: fitted(:), residual(:)
    real(wp) :: h_rms, h_max, coeff_l2, tail_l2, tiny_scale
    integer :: i, tail_lo

    allocate(fitted(n_pts), residual(n_pts))
    do i = 1, n_pts
      fitted(i) = cheb_get_val_point(n_deg, coeffs, a, b, x_nodes(i))
    end do
    residual = fitted - h_nodes

    tiny_scale = epsilon(1.0_wp)
    h_rms = sqrt(sum(h_nodes**2) / max(1, n_pts))
    h_max = maxval(abs(h_nodes))
    coeff_l2 = sqrt(sum(coeffs**2))

    stats%n_deg = n_deg
    stats%n_pts = n_pts
    stats%tail_modes = min(3, n_deg + 1)
    stats%solve_info = solve_info
    stats%used_dct = used_dct
    stats%rms_abs = sqrt(sum(residual**2) / max(1, n_pts))
    stats%rms_rel = stats%rms_abs / max(h_rms, tiny_scale)
    stats%linf_abs = maxval(abs(residual))
    stats%linf_rel = stats%linf_abs / max(h_max, tiny_scale)
    stats%coeff_max_abs = maxval(abs(coeffs))

    tail_lo = max(0, n_deg - stats%tail_modes + 1)
    stats%coeff_tail_max_abs = maxval(abs(coeffs(tail_lo:n_deg)))
    tail_l2 = sqrt(sum(coeffs(tail_lo:n_deg)**2))
    stats%coeff_tail_l2_ratio = tail_l2 / max(coeff_l2, tiny_scale)
    stats%coeff_tail_max_ratio = stats%coeff_tail_max_abs / max(stats%coeff_max_abs, tiny_scale)

    deallocate(fitted, residual)
  end subroutine fill_cheb_fit_stats

end module cheb_mod
