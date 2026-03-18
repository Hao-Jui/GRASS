module grid_mod
contains

subroutine make_grid
  use para_mod, only: SDIV, MDIV, s_gp, mu, DS, DM, pi
  implicit none
  integer :: m, s, i, n
  real(8) :: x, x_prev, p_n, dp_n

  do s = 1, SDIV
    s_gp(s) = (dble(s)-1.d0) * DS
  enddo
  do m = 1, MDIV
    mu(m) = (dble(m)-1.d0) * DM
  enddo
  mu(MDIV) = 1.d0

  ! Angular grid: Gauss-Legendre nodes mapped to [0,1]
  if (.false.) then
    n = MDIV
    do i = 1, (n+1)/2
      x = cos(pi * (dble(i) - 0.25d0) / (dble(n) + 0.5d0)) ! Initial guess
      do
        call legendre_and_deriv(n, x, p_n, dp_n)
        x_prev = x
        x = x_prev - p_n / dp_n
        if (abs(x - x_prev) < 1.d-14) exit
      end do
      mu(i)     = ( -x + 1.d0 ) / 2.d0  ! map from [-1,1] to [0,1], lower half
      mu(n+1-i) = (  x + 1.d0 ) / 2.d0  ! upper half
    end do
    if (mod(n,2) == 1) mu((n+1)/2) = 0.5d0
  endif
end subroutine make_grid

subroutine legendre_and_deriv(n, x, p, dp)
  integer, intent(in) :: n
  real(8), intent(in) :: x
  real(8), intent(out) :: p, dp
  integer :: k
  real(8) :: pkm1, pk, pkm2

  pkm2 = 1.d0
  pkm1 = x
  if (n == 0) then
    p = 1.d0
    dp = 0.d0
    return
  else if (n == 1) then
    p = x
    dp = 1.d0
    return
  end if

  do k = 2, n
    pk = ((2.d0*dble(k)-1.d0)*x*pkm1 - (dble(k)-1.d0)*pkm2) / dble(k)
    pkm2 = pkm1
    pkm1 = pk
  end do
  p = pk
  dp = dble(n) * (x * pk - pkm2) / (x*x - 1.d0)
end subroutine legendre_and_deriv

subroutine GridTrig
  use toolkit_mod, only: legendre, plgndr
  use para_mod, only: MDIV, LMAX, mu, sin_theta, P_2n, P1_2n_1, sin_2n_1_theta
  implicit none
  real(8) :: theta(MDIV)
  integer :: i, n

  P_2n(:,:)           = 0.0d0
  P1_2n_1(:,:)        = 0.0d0
  sin_2n_1_theta(:,:) = 0.0d0

  do i = 1, MDIV
    sin_theta(i) = sqrt(1.0d0 - mu(i)**2)
    theta(i) = asin(sin_theta(i))
  end do

  n = 0
  do i = 1, MDIV
    P_2n(i, n+1) = legendre(2*n, mu(i))
  end do

  do n = 1, LMAX
    do i = 1, MDIV
      P_2n(i, n+1)        = legendre(2*n, mu(i))
      P1_2n_1(i, n+1)     = plgndr(2*n - 1, 1, mu(i))
      sin_2n_1_theta(i, n) = sin((2.0d0 * n - 1.0d0) * theta(i))
    end do
  end do

end subroutine GridTrig

end module grid_mod
