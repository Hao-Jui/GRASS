module grid_mod
  use spectral_hub_mod, only: gauss_lobatto, chebyshev_lobatto_points, clenshaw_curtis_weights, barycentric_diff_matrices
  implicit none

contains

  subroutine make_grid
    use para_mod, only: SDIV, MDIV, s_gp, mu, DS, DM, pi, wp, &
                        angular_collocation, COLLOCATION_UNI, COLLOCATION_LEG, COLLOCATION_CHEB, &
                        D_mu, D_mu_t, D2_mu, w_mu
    implicit none
    integer  :: i, n
    real(wp) :: x, x_prev, p_n, dp_n

    s_gp = [(real(i, wp), i=0, SDIV-1)] * DS

    select case (angular_collocation)
    case (COLLOCATION_UNI)
      mu = [(real(i, wp), i=0, MDIV-1)] * DM
      mu(MDIV) = 1.0_wp
      D_mu_t = transpose(D_mu)
    case (COLLOCATION_LEG)
      call gauss_lobatto(MDIV, mu, w_mu)
      call barycentric_diff_matrices(mu, D_mu, D2_mu)
      D_mu_t = transpose(D_mu)
    case (COLLOCATION_CHEB)
      call chebyshev_lobatto_points(MDIV, mu)
      call clenshaw_curtis_weights(MDIV, w_mu)
      call barycentric_diff_matrices(mu, D_mu, D2_mu)
      D_mu_t = transpose(D_mu)
    case default
      error stop "make_grid: unknown collocation mode"
    end select

    ! Note: Consider using a dedicated library for this
    if (.false.) then
      n = MDIV
      do i = 1, (n + 1) / 2
        ! Cleaned up the initial guess math
        x = cos(pi * (real(i, wp) - 0.25_wp) / (real(n, wp) + 0.5_wp)) 
        
        ! Newton-Raphson Loop
        do
          call legendre_and_deriv(n, x, p_n, dp_n)
          x_prev = x
          x = x_prev - p_n / dp_n
          if (abs(x - x_prev) < 1.0e-14_wp) exit
        end do
        
        mu(i)        = (1.0_wp - x) * 0.5_wp 
        mu(n + 1 - i) = (1.0_wp + x) * 0.5_wp
      end do
      if (mod(n, 2) == 1) mu((n + 1) / 2) = 0.5_wp
    end if
  end subroutine make_grid

  subroutine legendre_and_deriv(n, x, p, dp)
    use para_mod, only: wp
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: x
    real(wp), intent(out) :: p, dp
    integer  :: k
    real(wp) :: pkm1, pk, pkm2

    ! Early exit for low orders
    if (n == 0) then
      p = 1.0_wp; dp = 0.0_wp
      return
    end if

    pkm2 = 1.0_wp
    pkm1 = x
    pk   = x ! Default for n=1

    do k = 2, n
      pk   = (real(2*k - 1, wp) * x * pkm1 - real(k - 1, wp) * pkm2) / real(k, wp)
      pkm2 = pkm1
      pkm1 = pk
    end do
    
    p  = pk
    dp = real(n, wp) * (x * pk - pkm2) / (x**2 - 1.0_wp)
  end subroutine legendre_and_deriv

  subroutine GridTrig
    use toolkit_mod, only: legendre, plgndr
    use para_mod,    only: MDIV, LMAX, mu, sin_theta, P_2n, P1_2n_1, &
                           sin_2n_1_theta, wp
    implicit none
    real(wp) :: theta(MDIV)
    integer  :: n, i

    sin_theta = sqrt(1.0_wp - mu**2)
    theta     = asin(sin_theta)

    P_2n           = 0.0_wp
    P1_2n_1        = 0.0_wp
    sin_2n_1_theta = 0.0_wp

    ! ":" syntax to apply the function to the whole column if 'legendre' is elemental or we need do loop
    P_2n(:, 1) = [(legendre(0, mu(i)), i=1, MDIV)]

    do n = 1, LMAX
      P_2n(:, n+1)           = [(legendre(2*n, mu(i)), i=1, MDIV)]
      P1_2n_1(:, n+1)        = [(plgndr(2*n - 1, 1, mu(i)), i=1, MDIV)]
      sin_2n_1_theta(:, n)   = sin(real(2*n - 1, wp) * theta)
    end do

  end subroutine GridTrig

end module grid_mod
