subroutine make_grid
  use para_mod
  implicit none
  integer :: m, s, i, n
  real(8) :: x, x_prev, p_n, dp_n

  interface
    subroutine legendre_and_deriv(n, x, p, dp)
      integer, intent(in) :: n
      real(8), intent(in) :: x
      real(8), intent(out) :: p, dp
    end subroutine legendre_and_deriv
  end interface

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
  !open(72,file="./trig/s_axis.dat")
  !write(72,"(es18.9)") s_gp
  !open(73,file="./trig/mu_axis.dat")
  !write(73,"(es18.9)") mu
  !close(72); close(73)
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

!===============================================================================
! Subroutine: GridTrig
!
! Purpose:
!   Initializes trigonometric and polynomial basis functions on the 
!   computational grid for solving Einstein field equations in rotating
!   neutron star configurations.
!
! Description:
!   This subroutine computes:
!   1. Legendre polynomials P_{2n}(mu) and associated Legendre P^1_{2n-1}(mu)
!   2. Sine functions sin((2n-1)*theta) for multipole expansions
!
! Grid Structure:
!   - s_gp(i): Radial grid points (i=1..SDIV)
!   - mu(m): Angular grid points, mu = cos(theta) (m=1..MDIV)
!   - n: Multipole order (n=0..LMAX)
!
!===============================================================================
subroutine GridTrig
  use toolkit_mod, only: legendre, plgndr
  use para_mod
  implicit none
  
  ! Local variables
  real(8) :: theta(MDIV)            ! Angular grid in theta coordinates
  
  !-----------------------------------------------------------------------------
  ! Initialize arrays to zero
  !-----------------------------------------------------------------------------
  call initialize_grid_arrays()
  
  !-----------------------------------------------------------------------------
  ! Compute angular basis functions
  !-----------------------------------------------------------------------------
  call compute_angular_basis_functions(theta)
  
end subroutine GridTrig


!===============================================================================
! Subroutine: initialize_grid_arrays
!
! Purpose:
!   Initializes all grid arrays to zero before computation.
!===============================================================================
subroutine initialize_grid_arrays()
  use para_mod
  implicit none
  
  P_2n(:,:)           = 0.0d0
  P1_2n_1(:,:)        = 0.0d0
  sin_2n_1_theta(:,:) = 0.0d0
  
end subroutine initialize_grid_arrays


!===============================================================================
! Subroutine: compute_angular_basis_functions
!
! Purpose:
!   Computes Legendre polynomials and trigonometric basis functions
!   on the angular grid.
!
! Basis functions computed:
!   - P_{2n}(mu): Even Legendre polynomials
!   - P^1_{2n-1}(mu): Associated Legendre polynomials
!   - sin((2n-1)*theta): Sine functions for odd multipoles
!===============================================================================
subroutine compute_angular_basis_functions(theta)
  use toolkit_mod, only: legendre, plgndr
  use para_mod
  implicit none
  
  real(8), intent(out) :: theta(MDIV)
  integer :: i, n, m
  
  !-----------------------------------------------------------------------------
  ! Compute angular coordinate theta from mu = cos(theta)
  !-----------------------------------------------------------------------------
  do m = 1, MDIV
    sin_theta(m) = sqrt(1.0d0 - mu(m)**2)
    theta(m) = asin(sin_theta(m))
  end do
  
  !-----------------------------------------------------------------------------
  ! Monopole term (n=0)
  !-----------------------------------------------------------------------------
  n = 0
  do i = 1, MDIV
    P_2n(i, n+1) = legendre(2*n, mu(i))
  end do
  
  !-----------------------------------------------------------------------------
  ! Higher multipoles (n >= 1)
  !-----------------------------------------------------------------------------
  do n = 1, LMAX
    do i = 1, MDIV
      ! Even Legendre polynomial P_{2n}
      P_2n(i, n+1) = legendre(2*n, mu(i))
      
      ! Associated Legendre polynomial P^1_{2n-1}
      P1_2n_1(i, n+1) = plgndr(2*n - 1, 1, mu(i))
      
      ! Sine function for odd multipole expansion
      sin_2n_1_theta(i, n) = sin((2.0d0 * n - 1.0d0) * theta(i))
    end do
  end do
  
end subroutine compute_angular_basis_functions
