! Spectral quadrature and polynomial helpers on [0, 1]
module spectral_hub_mod
  use precision_mod, only: wp
  implicit none
  private
  public :: legendre_sequence, gauss_legendre, gauss_lobatto, &
            chebyshev_lobatto_points, clenshaw_curtis_weights, barycentric_diff_matrices, &
            integrate_tabulated, integrate_at_gauss_nodes

  integer, parameter :: max_newton_iter = 50
  real(wp), parameter :: NEWTON_TOL = 64.0_wp * epsilon(1.0_wp)
contains
  pure subroutine legendre_sequence(n, z, pn, pnm1)
    integer, intent(in) :: n
    real(wp), intent(in) :: z
    real(wp), intent(out) :: pn, pnm1 ! P_{n}(z) and P_{n-1}(z) 
    integer :: j
    real(wp) :: p_km2, p_km1, pk

    if (n == 0) then
      pn = 1.0_wp;  pnm1 = 0.0_wp;  return
    end if

    p_km2 = 1.0_wp;  p_km1 = z
    if (n == 1) then
      pn = p_km1;  pnm1 = p_km2;  return
    end if

    do j = 2, n
      pk = ((2*j - 1) * z * p_km1 - (j - 1) * p_km2) / real(j, wp)
      p_km2 = p_km1;  p_km1 = pk
    end do
    pn = p_km1;  pnm1 = p_km2
  end subroutine legendre_sequence

  pure subroutine legendre_derivatives(n, z, pn, pnm1, dpn, ddpn)
    integer, intent(in) :: n
    real(wp), intent(in) :: z, pn, pnm1
    real(wp), intent(out) :: dpn, ddpn
    real(wp) :: one_minus_z2

    one_minus_z2 = 1.0_wp - z*z
    dpn = real(n, wp) * (pnm1 - z * pn) / one_minus_z2
    ddpn = (2.0_wp * z * dpn - real(n * (n + 1), wp) * pn) / one_minus_z2
  end subroutine legendre_derivatives

  ! ---- Shared transform: [-1, 1] -> [0, 1] ----
  pure subroutine to_unit_interval(n, x, w)
    integer, intent(in) :: n
    real(wp), intent(inout) :: x(n), w(n)
    x = (x + 1.0_wp) / 2.0_wp
    w = w / 2.0_wp
  end subroutine to_unit_interval

  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(wp), intent(out) :: x(n), w(n)
    integer :: i, iter, m
    real(wp) :: z, z1, pn, pnm1, dp, ddp, pi, scale

    if (n < 1) error stop "gauss_legendre: n must be at least 1"

    pi = 4.0_wp * atan(1.0_wp)
    m = (n + 1) / 2

    do i = 1, m
      z = cos(pi * (i - 0.25_wp) / (n + 0.5_wp))
      do iter = 1, max_newton_iter
        call legendre_sequence(n, z, pn, pnm1)
        call legendre_derivatives(n, z, pn, pnm1, dp, ddp)
        if (abs(dp) <= epsilon(1.0_wp) * max(1.0_wp, abs(pn))) &
          error stop "gauss_legendre: Newton step became singular"
        z1 = z;  z = z1 - pn / dp
        scale = max(1.0_wp, abs(z), abs(z1))
        if (abs(z - z1) <= newton_tol * scale) exit
      end do
      if (iter > max_newton_iter) error stop "gauss_legendre: Newton iteration did not converge"
      x(i) = -z;  x(n+1-i) = z
      w(i) = 2.0_wp / ((1.0_wp - z*z) * dp*dp);  w(n+1-i) = w(i)
    end do
    call to_unit_interval(n, x, w)
  end subroutine gauss_legendre

  subroutine gauss_lobatto(n, x, w)
    integer, intent(in)  :: n
    real(wp), intent(out) :: x(n), w(n)
    integer  :: i, iter, m, order
    real(wp) :: z, z1, dp, ddp, pn1, pnm2, pi, scale

    if (n < 2) error stop "gauss_lobatto: n must be at least 2"

    pi = 4.0_wp * atan(1.0_wp)
    order = n - 1;  m = (n + 1) / 2

    do i = 1, m
      if (i == 1) then
        z = -1.0_wp
      elseif (i == m .and. mod(n, 2) == 1) then
        z = 0.0_wp
      else
        z = -cos(pi * real(i - 1, wp) / real(order, wp))
        do iter = 1, max_newton_iter
          call legendre_sequence(order, z, pn1, pnm2)
          call legendre_derivatives(order, z, pn1, pnm2, dp, ddp)
          if (abs(ddp) <= epsilon(1.0_wp) * max(1.0_wp, abs(dp))) &
            error stop "gauss_lobatto: Newton step became singular"
          z1 = z;  z = z1 - dp / ddp
          scale = max(1.0_wp, abs(z), abs(z1))
          if (abs(z - z1) <= newton_tol * scale) exit
        end do
        if (iter > max_newton_iter) error stop "gauss_lobatto: Newton iteration did not converge"
      end if

      call legendre_sequence(order, z, pn1, pnm2)
      x(i) = z;  x(n+1-i) = -z
      w(i) = 2.0_wp / (real(order * n, wp) * pn1**2);  w(n+1-i) = w(i)
    end do
    call to_unit_interval(n, x, w)
  end subroutine gauss_lobatto

  pure subroutine chebyshev_lobatto_points(n, x)
    integer, intent(in) :: n
    real(wp), intent(out) :: x(n)
    integer :: i
    real(wp) :: pi

    if (n < 2) error stop "chebyshev_lobatto_points: n must be at least 2"

    pi = 4.0_wp * atan(1.0_wp)
    do i = 1, n
      x(i) = 0.5_wp * (1.0_wp - cos(pi * real(i - 1, wp) / real(n - 1, wp)))
    end do
  end subroutine chebyshev_lobatto_points

  ! Clenshaw-Curtis weights on [0, 1] for n Chebyshev-Lobatto points.
  ! Follows Trefethen, Spectral Methods in MATLAB (2000), clencurt.m.
  pure subroutine clenshaw_curtis_weights(n, w)
    integer, intent(in) :: n
    real(wp), intent(out) :: w(n)
    integer :: j, k, N_int, half
    real(wp) :: theta, pi_val

    if (n < 2) error stop "clenshaw_curtis_weights: n must be at least 2"

    pi_val = 4.0_wp * atan(1.0_wp)
    N_int = n - 1

    if (mod(N_int, 2) == 0) then
      w(1) = 1.0_wp / real(N_int*N_int - 1, wp)
    else
      w(1) = 1.0_wp / real(N_int*N_int, wp)
    end if
    w(n) = w(1)

    half = N_int / 2
    do j = 2, N_int
      theta = real(j - 1, wp) * pi_val / real(N_int, wp)
      w(j) = 1.0_wp
      if (mod(N_int, 2) == 0) then
        do k = 1, half - 1
          w(j) = w(j) - 2.0_wp * cos(2*k*theta) / real(4*k*k - 1, wp)
        end do
        w(j) = w(j) - cos(N_int * theta) / real(N_int*N_int - 1, wp)
      else
        do k = 1, half
          w(j) = w(j) - 2.0_wp * cos(2*k*theta) / real(4*k*k - 1, wp)
        end do
      end if
      w(j) = 2.0_wp * w(j) / real(N_int, wp)
    end do

    w = w * 0.5_wp  ! [-1,1] -> [0,1]
  end subroutine clenshaw_curtis_weights

  subroutine barycentric_diff_matrices(x, d1, d2)
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: d1(:,:), d2(:,:)
    real(wp) :: bary_w(size(x))
    integer :: i, j, n

    n = size(x)
    if (size(d1, 1) /= n .or. size(d1, 2) /= n) error stop "barycentric_diff_matrices: D1 size mismatch"
    if (size(d2, 1) /= n .or. size(d2, 2) /= n) error stop "barycentric_diff_matrices: D2 size mismatch"

    if (n <= 1) then
      d1 = 0.0_wp;  d2 = 0.0_wp;  return
    end if

    ! General barycentric weights: w_j = 1 / prod_{k/=j} (x_j - x_k)
    do i = 1, n
      bary_w(i) = 1.0_wp
      do j = 1, n
        if (j /= i) bary_w(i) = bary_w(i) / (x(i) - x(j))
      end do
    end do

    d1 = 0.0_wp
    do i = 1, n
      do j = 1, n
        if (j /= i) d1(i, j) = bary_w(j) / (bary_w(i) * (x(i) - x(j)))
      end do
      d1(i, i) = -sum(d1(i, :))
    end do

    d2 = matmul(d1, d1)
  end subroutine barycentric_diff_matrices

  function integrate_tabulated(f_values, a, b, n) result(integral)
    integer, intent(in) :: n
    real(wp), intent(in) :: f_values(n), a, b
    real(wp) :: integral
    real(wp) :: x(n), w(n)

    call gauss_legendre(n, x, w)
    integral = dot_product(w, f_values) * (b - a)
  end function integrate_tabulated

  subroutine integrate_at_gauss_nodes(f, a, b, n, cumulative)
    real(wp), intent(in) :: a, b
    integer, intent(in) :: n
    real(wp), intent(out) :: cumulative(n)
    interface
      function f(x)
        import :: wp
        real(wp), intent(in) :: x
        real(wp) :: f
      end function f
    end interface
    real(wp) :: x(n), w(n), x_phys(n)
    real(wp) :: scale, sub_scale
    integer :: i, j

    if (n < 1) error stop "integrate_at_gauss_nodes: n must be at least 1"

    call gauss_legendre(n, x, w)
    scale = b - a
    x_phys = a + scale * x

    do i = 1, n
      sub_scale = x_phys(i) - a
      cumulative(i) = 0.0_wp
      do j = 1, n
        cumulative(i) = cumulative(i) + w(j) * f(a + sub_scale * x(j))
      end do
      cumulative(i) = cumulative(i) * sub_scale
    end do
  end subroutine integrate_at_gauss_nodes

end module spectral_hub_mod
