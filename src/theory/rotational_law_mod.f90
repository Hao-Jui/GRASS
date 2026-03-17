module rotation_law_mod
  implicit none
  public
  ! Stored context for the Brent callback
  real(8), save :: ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h
  real(8), save :: ctx_rho_pole_h, ctx_gama_pole_h
  private :: AA_h, BB_h
contains
  subroutine diff_rotation_const_j(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
    use para_mod, only: A_diff
    real(8), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(8), intent(out):: fx
    real(8) :: LHS
    LHS =  re**2 * (g_h + rho_h - g_p - rho_p) + log(1.d0 - ((x-w_h)*exp(-re**2*rho_h))**2 )
    fx = LHS * (1.d0 - ((x-w_h)*exp(-re**2*rho_h))**2 )**2 &
        - ( (x-w_h) * exp(-2.d0*re**2*rho_h) / A_diff )**2
  end subroutine diff_rotation_const_j

  subroutine rotation_law_const_j(x, fx, re, rho_p, ww_p, sgp, mugp)
    use para_mod, only: A_diff, Omega_c
    real(8), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(8), intent(out):: fx
    fx = (Omega_c - x) * ( (1.d0-sgp)**2 - ( sgp*(x-ww_p)*exp(-re**2*rho_p) )**2 * (1.d0-mugp*mugp) ) &
        - ( (x-ww_p) * sgp * sgp * (1.d0-mugp*mugp) * exp(-2.d0 * re**2 * rho_p) ) / A_diff**2
  end subroutine rotation_law_const_j


  subroutine diff_rotation_uryu(x, fx, re, rho_h, g_h, w_h, rho_p, g_p)
  ! Compute Omega_e
    use para_mod, only: lambda2, Fmax_h
    implicit none
    real(8), intent(in) :: x, re, rho_h, g_h, w_h, rho_p, g_p
    real(8), intent(out):: fx
    real(8) :: F_e, ocre, RHS, aa, bb

    ocre = x / lambda2
    F_e  = (x - w_h) / ( exp(2.d0*re**2*rho_h) - (x-w_h)**2 )

    aa = AA_h(F_e, Fmax_h)
    bb = BB_h(F_e, Fmax_h)
    
    RHS = F_e * x - aa * ocre / 4.d0 * &
          (2.d0 * aa / bb * atan(F_e**2/aa**2) &
          - sqrt(2.d0) * ( atan(1.d0 - F_e*sqrt(2.d0)/aa) - atan(1.d0 + F_e*sqrt(2.d0)/aa) ) &
          + sqrt(2.d0) * ATANH( aa * F_e * sqrt(2.d0) / (F_e**2 + aa**2) ) )

    fx = re**2 * (g_h + rho_h - g_p - rho_p) + log( 1.d0 - ((x-w_h)*exp(-re**2*rho_h))**2 ) &
        + 2.d0 * RHS
    !write(*,"(A10,19es15.6)") "brent:", x, fx, ((x-w_e)*exp(-re**2*rho_e))**2
    !if (fx.ne.fx) stop "L536"

  end subroutine diff_rotation_uryu

  subroutine rotation_law_uryu(x, fx, re, rho_p, ww_p, sgp, mugp)
    use para_mod, only: uyru_p,uyru_q,Omega_c,F_equator_h, Fmax_h
    implicit none
    real(8), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
    real(8), intent(out):: fx
    real(8) :: Fj, aa, bb

    aa = AA_h(F_equator_h, Fmax_h)
    bb = BB_h(F_equator_h, Fmax_h)

    Fj = (x-ww_p) * sgp**2 * (1.d0-mugp**2) / ( (1.d0-sgp)**2 * exp(2.d0*re**2*rho_p) - (x-ww_p)**2 * sgp**2 * (1.d0-mugp**2) )
    fx = x / Omega_c * ( 1.d0 + (Fj / aa)**(uyru_p+uyru_q) ) - ( 1.d0 + (Fj / bb)**uyru_p )
    
    !write(*,"(10es15.6)") x,fx,tmp2,BB,AA

  end subroutine rotation_law_uryu

  elemental real(8) function intF(x, F_at_x)
    use para_mod, only: F_equator_h, Fmax_h, Omega_c
    implicit none
    real(8), intent(in) :: x, F_at_x
    real(8) :: aa, bb

    if ( x == 0.d0 .and. F_at_x == 0.d0 ) then
      intF = 0.d0
    else
      aa = AA_h( F_equator_h, Fmax_h )
      bb = BB_h( F_equator_h, Fmax_h )

      intF = F_at_x * x - aa * omega_c / 4.d0 * &
          (2.d0 * aa / bb * atan(F_at_x**2/aa**2) &
          - sqrt(2.d0) * ( atan(1.d0-F_at_x*sqrt(2.d0)/aa) - atan(1.d0+F_at_x*sqrt(2.d0)/aa) ) &
          + sqrt(2.d0) * ATANH(aa*F_at_x*sqrt(2.d0)/(F_at_x**2 + aa**2)) )
    endif

  end function intF

  pure elemental real(8) function AA_h(F_e, F_m)
  ! AA_h = A^2 * Omega_c
    use para_mod, only: lambda1, lambda2, uyru_p, uyru_q
    implicit none
    real(8), intent(in) :: F_e, F_m
    real(8) :: tmp1, tmp2
    
    tmp1 = lambda2 * F_e**uyru_q &
        - lambda1 * F_m**uyru_q
    tmp2 = (lambda1-1.d0) * F_e**uyru_p &
        - (lambda2-1.d0) * F_m**uyru_p
    AA_h = ( (F_e*F_m)**uyru_p * tmp1 / tmp2 )**(1.d0/dble(uyru_p+uyru_q))

  end function

  pure elemental real(8) function BB_h(F_e, F_m)
  ! BB_h = B^2 * Omega_c
    use para_mod, only: lambda1, lambda2, uyru_p, uyru_q
    implicit none
    real(8), intent(in) :: F_e, F_m
    real(8) :: tmp1, tmp2
    
    tmp1 = lambda2 * F_e**uyru_q &
        - lambda1 * F_m**uyru_q
    tmp2 = lambda2 * (lambda1-1.d0) * F_e**(uyru_p+uyru_q) &
        - lambda1 * (lambda2-1.d0) * F_m**(uyru_p+uyru_q)
    BB_h = F_e * F_m * ( tmp1 / tmp2 )**(1.d0/dble(uyru_p))

  end function

  subroutine set_shoot_context(rho_equator_h, gama_equator_h, ww_equator_h, rho_pole_h, gama_pole_h)
    real(8), intent(in) :: rho_equator_h, gama_equator_h, ww_equator_h, rho_pole_h, gama_pole_h
    ctx_rho_equator_h  = rho_equator_h
    ctx_gama_equator_h = gama_equator_h
    ctx_ww_equator_h   = ww_equator_h
    ctx_rho_pole_h     = rho_pole_h
    ctx_gama_pole_h    = gama_pole_h
  end subroutine set_shoot_context

  subroutine shoot_Fmax(current_Fmax_h, diff_sq)
    use para_mod, only: Fmax_h, Omega_e, Omega_c, r_e, rho, ww, F_equator_h, &
                        SDIV, s_gp, lambda1, lambda2
    use brent_mod, only: find_omege_e, zbrent_rot
    real(8), intent(in)  :: current_Fmax_h
    real(8), intent(out) :: diff_sq
    real(8) :: guess, Fa, rsm, wwsm, sgp, mum, omg_max_h, diff_Fmax
    real(8), dimension(SDIV) :: omg_mu_0
    integer :: imax, s
    real(8), parameter :: tolerance = 1.d-5

    ! Fmax_h is a module variable, so we set it here to be used by diff_rotation_uryu
    Fmax_h = current_Fmax_h

    guess = Omega_e
    call find_omege_e(guess, r_e, ctx_rho_equator_h, ctx_gama_equator_h, ctx_ww_equator_h, &
                    ctx_rho_pole_h, ctx_gama_pole_h, tolerance, Fa, diff_rotation_uryu)
    Omega_e = Fa
    F_equator_h = (Omega_e - ww(1,1)) / (exp(2.d0*r_e**2*rho(1,1)) - (Omega_e-ww(1,1))**2)
    if (F_equator_h < 0.d0) then
      diff_sq = 1.d3 ! Return a large value if F_equator_h is invalid
      return
    endif

    Omega_c = Omega_e / lambda2
    omg_mu_0(1) = Omega_c
    mum = 0.d0
    do s = 2, SDIV*2/3
        rsm = rho(s,1) ! hat
        wwsm= ww (s,1) ! hat
        sgp = s_gp(s)
        guess  = omg_mu_0(s-1)
        call zbrent_rot( guess, r_e, rsm, wwsm, sgp, mum, 1.d-5, omg_mu_0(s), rotation_law_uryu)
    enddo
    imax = maxloc(omg_mu_0(1:SDIV*2/3), 1)
    omg_max_h = omg_mu_0(imax)
    diff_Fmax = (lambda1 - omg_max_h / Omega_c)
    diff_sq = diff_Fmax**2
  end subroutine shoot_Fmax
end module rotation_law_mod
