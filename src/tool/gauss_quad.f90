! Gauss-Legendre nodes and weights are on [0, 1]
module gauss_quad
  use precision_mod, only: wp
  implicit none
contains
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(wp), intent(out) :: x(n), w(n)
    integer :: i, j, m
    real(wp) :: z, z1, p1, p2, p3, pp
    real(wp) :: pi

    pi = 4.e0_wp*atan(1.e0_wp)
    m = (n+1)/2

    ! Compute nodes and weights on [-1, 1]
    do i = 1, m
       z = cos(pi*(i-0.25e0_wp)/(n+0.5e0_wp))
       do
          p1 = 1.e0_wp
          p2 = 0.e0_wp
          do j = 1, n
             p3 = p2
             p2 = p1
             p1 = ((2.e0_wp*j-1.e0_wp)*z*p2 - (j-1.e0_wp)*p3)/j
          end do
          pp = n*(z*p1 - p2)/(z*z - 1.e0_wp)
          z1 = z
          z = z1 - p1/pp
          if (abs(z - z1) < 1.e-14_wp) exit
       end do
       x(i)     = -z
       x(n+1-i) =  z
       w(i)     = 2.e0_wp/((1.e0_wp - z*z)*pp*pp)
       w(n+1-i) = w(i)
    end do

    ! Transform from [-1, 1] to [0, 1]
    do i = 1, n
       x(i) = (x(i) + 1.e0_wp)/2.e0_wp
       w(i) = w(i)/2.e0_wp
    end do
  end subroutine gauss_legendre

  subroutine gauss_lobatto(n, x, w)
    ! n >= 2; nodes and weights on [0, 1]
    ! Endpoints x=0,1 are always included.
    ! Interior nodes are roots of P'_{n-1}; found via Newton on P'_{n-1}/P''_{n-1}.
    integer, intent(in)  :: n
    real(wp), intent(out) :: x(n), w(n)
    integer  :: i, j, m
    real(wp) :: z, z1, p_km2, p_km1, pk, dp, ddp, pn1
    real(wp) :: pi

    pi = 4.e0_wp * atan(1.e0_wp)
    m  = (n + 1) / 2     ! exploit symmetry

    do i = 1, m
      if (i == 1) then
        z = -1.e0_wp
      else
        z = -cos(pi * real(i-1, wp) / real(n-1, wp))   ! Chebyshev initial guess
        do
          p_km2 = 1.e0_wp
          p_km1 = z
          do j = 2, n-1
            pk    = ((2*j - 1) * z * p_km1 - (j-1) * p_km2) / j
            p_km2 = p_km1
            p_km1 = pk
          end do
          ! p_km1 = P_{n-1}(z),  p_km2 = P_{n-2}(z)
          dp  =  real(n-1, wp) * (z * p_km1 - p_km2) / (z*z - 1.e0_wp)   ! P'_{n-1}
          ddp = (2.e0_wp * z * dp - real(n*(n-1), wp) * p_km1) / (1.e0_wp - z*z)  ! P''_{n-1}
          z1  = z
          z   = z1 - dp / ddp
          if (abs(z - z1) < 1.e-14_wp) exit
        end do
        pn1 = p_km1   ! P_{n-1}(z) from last Newton step
      end if

      if (i == 1) then
        p_km2 = 1.e0_wp
        p_km1 = z
        do j = 2, n-1
          pk    = ((2*j-1)*z*p_km1 - (j-1)*p_km2) / j
          p_km2 = p_km1
          p_km1 = pk
        end do
        pn1 = p_km1
      end if

      x(i)     =  z
      x(n+1-i) = -z
      w(i)     = 2.e0_wp / (real((n-1)*n, wp) * pn1**2)
      w(n+1-i) = w(i)
    end do

    ! Transform from [-1, 1] to [0, 1]
    do i = 1, n
      x(i) = (x(i) + 1.e0_wp) / 2.e0_wp
      w(i) = w(i) / 2.e0_wp
    end do
  end subroutine gauss_lobatto

  ! Integrate pre-tabulated function values at Gauss-Legendre nodes
  function integrate_tabulated(f_values, a, b, n) result(integral)
    integer, intent(in) :: n
    real(wp), intent(in) :: f_values(n), a, b
    real(wp) :: integral
    real(wp) :: x(n), w(n), scale
    integer :: i

    ! Get Gauss-Legendre weights on [0, 1]
    call gauss_legendre(n, x, w)

    ! Scale for interval [a, b]
    scale = (b - a)

    ! Direct weighted sum
    integral = 0.e0_wp
    do i = 1, n
       integral = integral + w(i)*f_values(i)
    end do
    integral = integral*scale
  end function integrate_tabulated

  ! Compute cumulative integrals at the Gauss nodes themselves
  ! Returns: cumulative(k) = integral from a to the k-th Gauss node
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
    real(wp) :: x(n), w(n), x_transformed(n)
    real(wp) :: scale
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
