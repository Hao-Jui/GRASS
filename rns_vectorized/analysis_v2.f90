subroutine mass_radius
  use toolkit_mod, only: deriv_s_1d, interp
  use para_mod
  implicit none
  integer :: s, m
  real(8) :: s1, s_1, s_p, r_h
  real(8) :: d_gama_s, d_rho_s, d_ww_s
  real(8) :: gama_equator, rho_equator, ww_equator, gama_pole, rho_pole
  real(8) :: J, E_tilde_m, l_tilde_m, E_tilde_p, l_tilde_p
  real(8), dimension(MDIV) :: Intm
  real(8), dimension(SDIV) :: Intr
  real(8), dimension(SDIV) :: D_m, D_m_0, D_J, D_T, D_m_p
  real(8), dimension(SDIV) :: d_o_e, d_g_e, d_r_e, d_v_e ! For Kepler properties
  real(8), dimension(SDIV) :: dd_o, dd_g, dd_r, dd_v 
  real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1
  real(8) :: n0_at_e
  real(8) :: er2
  real(8) :: sqrt_v, doe, dge, dre, dve, vek
  real(8) :: o_r, o_rr, g_r, g_rr, r_r, r_rr
  integer :: ifail
  real(8), dimension(SDIV,MDIV) :: rho_0
  character(32) :: fil1, fil2, fil3

  s_p = r_ratio**(1.d0/dble(s_pwr)) / ( 1.d0 + r_ratio**(1.d0/dble(s_pwr)) )

  if (output) then
    write(fil1,"(f6.4)") parA 
    write(fil2,"(f6.4)") parB
    write(fil3,"(f6.4)") chi
    open(21, file="./Cont/Omega.dat")
    do s = 1, 2*SDIV/3
      write(21,"(99es18.9)") (s_gp(s)/(1.d0-s_gp(s)))**s_pwr, omg(s,1)/2.d0/pi* (C/sqrt(kappa)), &
              enthalpy(s,1), F_j(s,1)
    enddo
    close(21)
  endif

  do s = 1, SDIV
    gama_mu_1(s) = gama(s,MDIV)
    rho_mu_1 (s) = rho (s,MDIV)
    gama_mu_0(s) = gama(s,1)
    rho_mu_0 (s) = rho (s,1)
    ww_mu_0  (s) = ww  (s,1)
  enddo
  call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole)
  call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator)
  call interp(s_gp,  rho_mu_1, SDIV, s_p, rho_pole)
  call interp(s_gp,  rho_mu_0, SDIV, s_e, rho_equator)
  
  ! Masses and angular momentum
  Mass   = 0.d0
  mass_0 = 0.d0
  J      = 0.d0
  T_kin  = 0.d0

  do s = 1, SDIV
    do m = 1, MDIV
      if (energy(s,m) > e_surface) then
        rho_0(s,m) = n0_at_e(energy(s,m))
        rho_0(s,m) = rho_0(s,m) * MB * KSCALE*C**2
      else
        rho_0(s,m) = 0.d0
      endif
    enddo
  enddo

  do s = 1, SDIV
    s1 = ( s_gp(s) / (1.d0 - s_gp(s)) )**s_pwr

    Intm = exp( 2.d0 * alpha(s,:) + gama(s,:) ) * &
        ( ( (energy(s,:)+pressure(s,:)) / (1.d0-velocity_sq(s,:)) ) * &
        ( 1.d0 + velocity_sq(s,:) +  2.d0 * sqrt(velocity_sq(s,:) ) * s1 &
        * sqrt(1.d0-mu(:)**2) * r_e * ww(s,:) * exp( -rho(s,:) ) ) + 2.d0*pressure(s,:) )
    call d01gaf( mu, Intm, MDIV, D_m(s), er2, ifail)

    Intm = exp( 2.d0*alpha(s,:) + (gama(s,:) - rho(s,:)) / 2.d0 ) &
          * rho_0(s,:) / sqrt(1.d0-velocity_sq(s,:))
    call d01gaf( mu, Intm, MDIV, D_m_0(s), er2, ifail)

    Intm = exp( 2.d0*alpha(s,:) + (gama(s,:) - rho(s,:)) / 2.d0 ) &
          * (energy(s,:) ) / sqrt(1.d0-velocity_sq(s,:))
    call d01gaf( mu, Intm, MDIV, D_m_p(s), er2, ifail)

    Intm = sqrt(1.d0-mu(:)**2) * exp(2.d0*alpha(s,:)+gama(s,:)-rho(s,:)) &
          * (energy(s,:) + pressure(s,:)) * sqrt(velocity_sq(s,:)) / (1.d0-velocity_sq(s,:))
    call d01gaf( mu, Intm, MDIV, D_J(s), er2, ifail)

    Intm = Intm * omg(s,:)
    call d01gaf( mu, Intm, MDIV, D_T(s), er2, ifail)
  enddo

  Intr = ( s_gp(:) / (1.d0-s_gp(:)) )**(3*s_pwr-1) / (1.d0-s_gp(:))**2 * D_m(:) * dble(s_pwr)
  call d01gaf( s_gp, Intr, SDIV, Mass, er2, ifail)

  Intr = ( s_gp(:) / (1.d0-s_gp(:)) )**(3*s_pwr-1) / (1.d0-s_gp(:))**2 * D_m_0(:) * dble(s_pwr)
  call d01gaf( s_gp, Intr, SDIV, Mass_0, er2, ifail)

  Intr = ( s_gp(:) / (1.d0-s_gp(:)) )**(3*s_pwr-1) / (1.d0-s_gp(:))**2 * D_m_p(:) * dble(s_pwr)
  call d01gaf( s_gp, Intr, SDIV, Mass_p, er2, ifail)

  Intr = ( s_gp(:) / (1.d0-s_gp(:)) )**(4*s_pwr-1) / (1.d0-s_gp(:))**2 * D_J(:) * dble(s_pwr)
  call d01gaf( s_gp, Intr, SDIV, J, er2, ifail)

  Intr = ( s_gp(:) / (1.d0-s_gp(:)) )**(4*s_pwr-1) / (1.d0-s_gp(:))**2 * D_T(:) * dble(s_pwr)
  call d01gaf( s_gp, Intr, SDIV, T_kin, er2, ifail)
  
  Mass   = Mass   * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3/G
  Mass_0 = Mass_0 * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3/G
  Mass_p = Mass_p * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3/G
  
  J = J * 4.d0 * pi * kappa * C**3 * r_e**4 / G
  T_kin = T_kin * 2.d0 * pi * sqrt(kappa) * C**2 * r_e**4 / G

  ang_mom = J*C/(G*MSUN**2)
  chi = J*C/(G*Mass**2)
  



  ! Compute the velocities of co-rotating and counter-rotating particles
  ! with respect to a ZAMO
  do s = 1+(SDIV-1)/2, SDIV
    s1 = s_gp(s) * ( 1.d0 - s_gp(s) )
    s_1 = 1.d0 - s_gp(s)

    d_rho_s  = deriv_s_1d( rho,s)
    d_gama_s = deriv_s_1d(gama,s)
    d_ww_s   = deriv_s_1d(  ww,s)

    sqrt_v = exp(-2.d0*rho(s,1)) * r_e**2 * s_gp(s)**4 * d_ww_s**2 &
             + 2.d0 * s1 * (d_gama_s+d_rho_s) + s1**2 * (d_gama_s * d_gama_s - d_rho_s * d_rho_s)

    if (sqrt_v > 0.d0) then
      sqrt_v = sqrt(sqrt_v)
    else
      sqrt_v = 0.d0
    endif

    v_plus(s) = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_ww_s + sqrt_v) / &
              ( 2.d0 + s1 * (d_gama_s - d_rho_s) )

    v_minus(s)= ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_ww_s - sqrt_v) / &
              ( 2.d0 + s1 * (d_gama_s - d_rho_s) )
  enddo


  ! Kepler angular velocity
  do s = 1 , SDIV
    d_r_e(s)   = deriv_s_1d( rho,s)
    d_g_e(s)   = deriv_s_1d(gama,s)
    d_o_e(s)   = deriv_s_1d(  ww,s)
    d_v_e(s)   = deriv_s_1d( sqrt(velocity_sq) ,s)
    ww_mu_0(s) = ww(s,1)
  enddo

  open(23, file="./Cont/velocity.dat")
  do s = 1 , SDIV
    s1  = s_gp(s) * ( 1.d0 - s_gp(s) )
    s_1 = 1.d0 - s_gp(s)
    r_h = r_e * s_gp(s) / ( 1.d0 - s_gp(s) )
    dd_r(s)   = deriv_s_1d(d_r_e,s)
    dd_g(s)   = deriv_s_1d(d_g_e,s)
    dd_o(s)   = deriv_s_1d(d_o_e,s)

    o_r  = r_e * s_gp(s)**2 / r_h**2 * d_o_e(s)
    g_r  = r_e * s_gp(s)**2 / r_h**2 * d_g_e(s)
    r_r  = r_e * s_gp(s)**2 / r_h**2 * d_r_e(s)

    o_rr = - 2.d0 * r_e / (r_e+r_h)**3 * d_o_e(s) + r_e**2 * s_gp(s)**4 / r_h**4 * dd_o(s)
    g_rr = - 2.d0 * r_e / (r_e+r_h)**3 * d_g_e(s) + r_e**2 * s_gp(s)**4 / r_h**4 * dd_g(s)
    r_rr = - 2.d0 * r_e / (r_e+r_h)**3 * d_r_e(s) + r_e**2 * s_gp(s)**4 / r_h**4 * dd_r(s)

    V_rr_p(s) = v_plus (s)**2 * ( r_h**2 * (g_rr - r_rr + g_r**2 - r_r**2) &
              + 2.d0 * ( exp(-rho_mu_0(s)) * r_h**2 * o_r )**2 + 4.d0 * r_h * r_r - 6.d0 ) &
              + 2.d0 * v_plus (s) * r_h**3 * exp(-rho_mu_0(s)) * ( 2.d0 * r_r * o_r - o_rr ) &
              - r_h**2 * ( g_rr + r_rr + g_r**2 - r_r**2 )

    V_rr_p(s) = V_rr_p(s) / r_h**2 / (1.d0-v_plus(s)**2) * exp(gama_mu_0(s))

    V_rr_m(s) = v_minus(s)**2 * ( r_h**2 * (g_rr - r_rr + g_r**2 - r_r**2) &
              + 2.d0 * ( exp(-rho_mu_0(s)) * r_h**2 * o_r )**2 + 4.d0 * r_h * r_r - 6.d0 ) &
              + 2.d0 * v_minus(s) * r_h**3 * exp(-rho_mu_0(s)) * ( 2.d0 * r_r * o_r - o_rr ) &
              - r_h**2 * ( g_rr + r_rr + g_r**2 - r_r**2 )

    V_rr_m(s) = V_rr_m(s) / r_h**2 / (1.d0-v_minus(s)**2) * exp(gama_mu_0(s))

    l_tilde_m = v_minus(s) * r_h * exp( (gama_mu_0(s) + rho_mu_0(s)) / 2.d0 ) / sqrt( 1.d0 - v_minus(s)**2 )
    E_tilde_m = exp( (gama_mu_0(s) + rho_mu_0(s)) / 2.d0 ) / sqrt( 1.d0 - v_minus(s)**2 ) + ww_mu_0(s) * l_tilde_m
    
    write(23,"(10es18.9)") s_gp(s), v_plus(s), V_rr_p(s), v_minus(s), V_rr_m(s), l_tilde_m*sqrt(KAPPA)/1.d5, E_tilde_m
  enddo
  close(23)

  call interp(s_gp, d_o_e, SDIV, s_e, doe)
  call interp(s_gp, d_g_e, SDIV, s_e, dge)
  call interp(s_gp, d_r_e, SDIV, s_e, dre)
  call interp(s_gp, d_v_e, SDIV, s_e, dve)

  vek = ( doe / ( 8.d0 + dge - dre ) ) * r_e * exp(-rho_equator) &
      + sqrt( ( (dge+dre) / ( 8.d0 + dge - dre ) ) + ( ( doe / ( 8.d0 + dge - dre ) ) * r_e * exp(-rho_equator) )**2  )

  if ( r_ratio == 1.d0 ) then
    ww_equator = 0.d0
  else
    call interp(s_gp, ww_mu_0, SDIV, s_e, ww_equator);
  endif

  Omega_K = (C/sqrt(kappa)) * ( ww_equator + vek * exp(rho_equator) / r_e )
  !write(*,"(A10,f11.5,A7,f11.5,A4)") "Keplr :",Omega_K/2/pi,"Hz = ",Omega_K,"Hz"
  
  ! Circumferential radius
  r_circ = sqrt(kappa)*r_e*exp((gama_equator-rho_equator)/2.0)
  

end  subroutine mass_radius



