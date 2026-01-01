! Gaussian quadrature module (Fortran 90)
! Gauss-Legendre nodes and weights are on [0, 1]
module gauss_quad
  implicit none
contains
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(8), intent(out) :: x(n), w(n)
    integer :: i, j, m
    real(8) :: z, z1, p1, p2, p3, pp
    real(8) :: pi
    
    pi = 4.0d0*atan(1.0d0)
    m = (n+1)/2
    
    ! Compute nodes and weights on [-1, 1]
    do i = 1, m
       z = cos(pi*(i-0.25d0)/(n+0.5d0))
       do
          p1 = 1.0d0
          p2 = 0.0d0
          do j = 1, n
             p3 = p2
             p2 = p1
             p1 = ((2.0d0*j-1.0d0)*z*p2 - (j-1.0d0)*p3)/j
          end do
          pp = n*(z*p1 - p2)/(z*z - 1.0d0)
          z1 = z
          z = z1 - p1/pp
          if (abs(z - z1) < 1.0d-14) exit
       end do
       x(i)     = -z
       x(n+1-i) =  z
       w(i)     = 2.0d0/((1.0d0 - z*z)*pp*pp)
       w(n+1-i) = w(i)
    end do
    
    ! Transform from [-1, 1] to [0, 1]
    do i = 1, n
       x(i) = (x(i) + 1.0d0)/2.0d0
       w(i) = w(i)/2.0d0
    end do
  end subroutine gauss_legendre

  ! Integrate function f over [a, b] using n-point Gauss-Legendre quadrature
  function integrate(f, a, b, n) result(integral)
    real(8), intent(in) :: a, b
    integer, intent(in) :: n
    real(8) :: integral
    interface
       function f(x)
         real(8), intent(in) :: x
         real(8) :: f
       end function f
    end interface
    real(8) :: x(n), w(n), xi
    real(8) :: scale
    integer :: i
    
    ! Get Gauss-Legendre nodes and weights on [0, 1]
    call gauss_legendre(n, x, w)
    
    ! Scale for interval [a, b]
    scale = (b - a)
    
    ! Compute integral
    integral = 0.0d0
    do i = 1, n
       xi = a + scale*x(i)
       integral = integral + w(i)*f(xi)
    end do
    integral = integral*scale
  end function integrate

  ! Integrate pre-tabulated function values at Gauss-Legendre nodes
  function integrate_tabulated(f_values, a, b, n) result(integral)
    integer, intent(in) :: n
    real(8), intent(in) :: f_values(n), a, b
    real(8) :: integral
    real(8) :: x(n), w(n), scale
    integer :: i
    
    ! Get Gauss-Legendre weights on [0, 1]
    call gauss_legendre(n, x, w)
    
    ! Scale for interval [a, b]
    scale = (b - a)
    
    ! Direct weighted sum
    integral = 0.0d0
    do i = 1, n
       integral = integral + w(i)*f_values(i)
    end do
    integral = integral*scale
  end function integrate_tabulated

  ! Compute cumulative integrals from a to each point in x_eval
  subroutine integrate_cumulative(f, a, x_eval, n_eval, n_quad, cumulative)
    real(8), intent(in) :: a
    integer, intent(in) :: n_eval, n_quad
    real(8), intent(in) :: x_eval(n_eval)
    real(8), intent(out) :: cumulative(n_eval)
    interface
       function f(x)
         real(8), intent(in) :: x
         real(8) :: f
       end function f
    end interface
    integer :: i
    
    ! Compute integral from a to each evaluation point
    do i = 1, n_eval
       cumulative(i) = integrate(f, a, x_eval(i), n_quad)
    end do
  end subroutine integrate_cumulative

  ! Compute cumulative integrals at the Gauss nodes themselves
  ! Returns: cumulative(k) = integral from a to the k-th Gauss node
  subroutine integrate_at_gauss_nodes(f, a, b, n, cumulative)
    real(8), intent(in) :: a, b
    integer, intent(in) :: n
    real(8), intent(out) :: cumulative(n)
    interface
       function f(x)
         real(8), intent(in) :: x
         real(8) :: f
       end function f
    end interface
    real(8) :: x(n), w(n), x_transformed(n)
    real(8) :: scale
    integer :: i
    
    ! Get Gauss nodes on [0, 1]
    call gauss_legendre(n, x, w)
    
    ! Transform to [a, b]
    scale = (b - a)
    do i = 1, n
       x_transformed(i) = a + scale*x(i)
    end do
    
    ! Compute cumulative integral to each Gauss node
    do i = 1, n
       cumulative(i) = integrate(f, a, x_transformed(i), n)
    end do
  end subroutine integrate_at_gauss_nodes

end module gauss_quad