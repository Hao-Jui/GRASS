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
!   3. Radial integration kernels f_rho and f_gamma for metric potentials
!   4. Grid transformation functions f2n based on coordinate mapping
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
  real(8) :: f2n(LMAX+1, SDIV)     ! Radial transformation function
  real(8) :: theta(MDIV)            ! Angular grid in theta coordinates
  integer :: i, j, k, n, m          ! Loop indices
  real(8) :: s_j, s_k               ! Current radial grid points
  real(8) :: s_j_complement         ! 1 - s_j
  real(8) :: s_k_complement         ! 1 - s_k
  
  !-----------------------------------------------------------------------------
  ! Initialize arrays to zero
  !-----------------------------------------------------------------------------
  call initialize_grid_arrays()
  
  !-----------------------------------------------------------------------------
  ! Compute radial transformation functions f2n
  !-----------------------------------------------------------------------------
  call compute_radial_transformation(f2n)
  
  !-----------------------------------------------------------------------------
  ! Compute integration kernels for metric potentials
  !-----------------------------------------------------------------------------
  call compute_metric_kernels(f2n)
  
  !-----------------------------------------------------------------------------
  ! Set boundary conditions for radial kernels
  !-----------------------------------------------------------------------------
  call set_radial_boundary_conditions()
  
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
  f_rho(:,:,:)        = 0.0d0
  f_gama(:,:,:)       = 0.0d0
  
end subroutine initialize_grid_arrays


!===============================================================================
! Subroutine: compute_radial_transformation
!
! Purpose:
!   Computes the radial transformation function f2n(n,s) used in 
!   coordinate mapping.
!
! Formula:
!   f2n(n,s) = [(1-s)/s]^(2*s_pwr*n)
!===============================================================================
subroutine compute_radial_transformation(f2n)
  use para_mod
  implicit none
  
  real(8), intent(out) :: f2n(LMAX+1, SDIV)
  integer :: n, i
  real(8) :: s_ratio
  
  f2n(:,:) = 0.0d0
  
  do n = 0, LMAX
    do i = 2, SDIV
      s_ratio = (1.0d0 - s_gp(i)) / s_gp(i)
      f2n(n+1, i) = s_ratio**(2 * s_pwr * n)
    end do
  end do
  
end subroutine compute_radial_transformation


!===============================================================================
! Subroutine: compute_metric_kernels
!
! Purpose:
!   Computes integration kernels f_rho and f_gamma for the metric potentials
!   rho and gamma across the entire computational grid.
!
! Arguments:
!   f2n - Radial transformation functions
!
! Grid indices:
!   j - source point index (observer location)
!   k - integration point index (field point)
!   n - multipole order
!===============================================================================
subroutine compute_metric_kernels(f2n)
  use para_mod
  implicit none
  
  real(8), intent(in) :: f2n(LMAX+1, SDIV)
  integer :: j, k, n
  
  ! Main interior region: both j and k away from boundaries
  do j = 2, SDIV
    do n = 1, LMAX
      do k = 2, SDIV
        if (k < j) then
          ! Integration point inside observer point
          call compute_kernel_interior_region(j, k, n, f2n)
        else
          ! Integration point outside observer point
          call compute_kernel_exterior_region(j, k, n, f2n)
        end if
      end do
    end do
  end do
  
  ! Special case: n=0 (monopole term)
  call compute_monopole_kernels(f2n)
  
end subroutine compute_metric_kernels


!===============================================================================
! Subroutine: compute_kernel_interior_region
!
! Purpose:
!   Computes f_rho and f_gamma when the integration point is interior
!   to the observer point (k < j).
!===============================================================================
subroutine compute_kernel_interior_region(j, k, n, f2n)
  use para_mod
  implicit none
  
  integer, intent(in) :: j, k, n
  real(8), intent(in) :: f2n(LMAX+1, SDIV)
  
  real(8) :: s_j, s_k, s_j_comp, s_k_comp
  real(8) :: prefactor
  
  s_j = s_gp(j)
  s_k = s_gp(k)
  s_j_comp = 1.0d0 - s_j
  s_k_comp = 1.0d0 - s_k
  
  prefactor = dble(s_pwr)
  
  ! f_rho kernel for interior region
  f_rho(j, n+1, k) = f2n(n+1, j) * (s_j_comp / s_j)**s_pwr &
                   * s_k**(s_pwr - 1) / (f2n(n+1, k) * s_k_comp**(s_pwr + 1)) &
                   * prefactor
  
  ! f_gamma kernel for interior region
  f_gama(j, n+1, k) = f2n(n+1, j) / (f2n(n+1, k) * s_k * s_k_comp) &
                    * prefactor
  
end subroutine compute_kernel_interior_region


!===============================================================================
! Subroutine: compute_kernel_exterior_region
!
! Purpose:
!   Computes f_rho and f_gamma when the integration point is exterior
!   to the observer point (k >= j).
!===============================================================================
subroutine compute_kernel_exterior_region(j, k, n, f2n)
  use para_mod
  implicit none
  
  integer, intent(in) :: j, k, n
  real(8), intent(in) :: f2n(LMAX+1, SDIV)
  
  real(8) :: s_j, s_k, s_j_comp, s_k_comp
  real(8) :: prefactor
  
  s_j = s_gp(j)
  s_k = s_gp(k)
  s_j_comp = 1.0d0 - s_j
  s_k_comp = 1.0d0 - s_k
  
  prefactor = dble(s_pwr)
  
  ! f_rho kernel for exterior region
  f_rho(j, n+1, k) = f2n(n+1, k) / (f2n(n+1, j) * s_k * s_k_comp) &
                   * prefactor
  
  ! f_gamma kernel for exterior region
  f_gama(j, n+1, k) = f2n(n+1, k) * s_j_comp**(2 * s_pwr) &
                    * s_k**(2 * s_pwr - 1) &
                    / (s_j**(2 * s_pwr) * f2n(n+1, j) * s_k_comp**(2 * s_pwr + 1)) &
                    * prefactor
  
end subroutine compute_kernel_exterior_region


!===============================================================================
! Subroutine: compute_monopole_kernels
!
! Purpose:
!   Computes special kernels for the monopole term (n=0), which requires
!   separate treatment due to its different radial dependence.
!===============================================================================
subroutine compute_monopole_kernels(f2n)
  use para_mod
  implicit none
  
  real(8), intent(in) :: f2n(LMAX+1, SDIV)
  integer :: j, k, n
  real(8) :: s_j, s_k, s_j_comp, s_k_comp
  
  n = 0
  
  ! Main grid points (j >= 2)
  do j = 2, SDIV
    do k = 2, SDIV
      s_j = s_gp(j)
      s_k = s_gp(k)
      s_j_comp = 1.0d0 - s_j
      s_k_comp = 1.0d0 - s_k
      
      if (k < j) then
        f_rho(j, n+1, k) = (s_j_comp / s_j)**s_pwr &
                         * s_k**(s_pwr - 1) / s_k_comp**(s_pwr + 1) &
                         * dble(s_pwr)
      else
        f_rho(j, n+1, k) = dble(s_pwr) / (s_k * s_k_comp)
      end if
    end do
  end do
  
end subroutine compute_monopole_kernels


!===============================================================================
! Subroutine: set_radial_boundary_conditions
!
! Purpose:
!   Sets boundary conditions for the integration kernels at s=0 and s'=0.
!===============================================================================
subroutine set_radial_boundary_conditions()
  use para_mod
  implicit none
  
  integer :: j, k, n
  real(8) :: s_k, s_k_comp
  
  !-----------------------------------------------------------------------------
  ! Boundary condition at s = 0 (center of star, j=1)
  !-----------------------------------------------------------------------------
  j = 1
  
  ! n = 0
  n = 0
  do k = 2, SDIV
    s_k = s_gp(k)
    s_k_comp = 1.0d0 - s_k
    f_rho(j, n+1, k) = dble(s_pwr) / (s_k * s_k_comp)
  end do
  
  ! n = 1
  n = 1
  do k = 2, SDIV
    s_k = s_gp(k)
    s_k_comp = 1.0d0 - s_k
    f_rho(j, n+1, k)  = 0.0d0
    f_gama(j, n+1, k) = dble(s_pwr) / (s_k * s_k_comp)
  end do
  
  ! n >= 2
  do n = 2, LMAX
    do k = 1, SDIV
      f_rho(j, n+1, k)  = 0.0d0
      f_gama(j, n+1, k) = 0.0d0
    end do
  end do
  
  !-----------------------------------------------------------------------------
  ! Boundary condition at s' = 0 (integration boundary, k=1)
  !-----------------------------------------------------------------------------
  k = 1
  
  do j = 1, SDIV
    ! n = 0
    f_rho(j, 1, k) = 0.0d0
    
    ! n >= 1
    do n = 1, LMAX
      f_rho(j, n+1, k)  = 0.0d0
      f_gama(j, n+1, k) = 0.0d0
    end do
  end do
  
end subroutine set_radial_boundary_conditions


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