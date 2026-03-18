module cheb_mod
  use precision_mod, only: wp
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_ptr
  implicit none
  public
  
  real(wp), parameter, private :: pi = acos(-1.e0_wp)
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
  subroutine cheb_std_base(n_deg, n_pts, x_nodes, h_nodes, a, b, coeffs)
    integer, intent(in)  :: n_deg, n_pts
    real(wp), intent(in)  :: x_nodes(n_pts), h_nodes(n_pts), a, b
    real(wp), intent(out) :: coeffs(0:n_deg)
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

end module cheb_mod
