subroutine mass_radius

  use para_mod
  implicit none
  integer :: s,m
  real(8) :: s1,s_1,d_gama_s,d_rho_s,d_ww_s,sqrt_v, &
            doe,dge,dre,dve,vek,r_p,s_p
  real(8) :: gama_equator,rho_equator,ww_equator, &
            gama_pole, rho_pole
  real(8) :: J
  real(8), dimension(MDIV) :: Int
  real(8), dimension(SDIV) :: D_m, D_m_0, D_J, D_T, D_m_p
  real(8), dimension(SDIV) :: d_o_e,d_g_e,d_r_e,d_v_e ! For Kepler properties
  real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1
  real(8) :: n0_at_e,deriv_s
  real(8) :: er2
  integer :: ifail
  real(8), dimension(SDIV,MDIV) :: rho_0

  r_p = r_ratio*r_e
  s_p = r_p/(r_p+r_e)

  if (output) then
    open(21, file="./Cont/Omega.dat")
    do s = 1, 2*SDIV/3
      write(21,"(99es18.9)") (s_gp(s)/(1.d0-s_gp(s))), omg(s,1)/2.d0/pi* (C/sqrt(kappa)), &
              enthalpy(s,1)
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
    Int = exp(2.d0*alpha(s,:)+gama(s,:)) * &
      ( ( (energy(s,:)+pressure(s,:))/(1.d0-velocity_sq(s,:)) ) * &
        (1.d0+velocity_sq(s,:)+(2.d0*s_gp(s)*sqrt(velocity_sq(s,:))/ &
        (1.d0-s_gp(s)))*sqrt(1.d0-mu(:)**2) * r_e * ww(s,:) * &
        exp(-rho(s,:))) + 2.d0*pressure(s,:))
    call d01gaf( mu, Int, MDIV, D_m(s), er2, ifail)

    Int = exp( 2.d0*alpha(s,:) + (gama(s,:) - rho(s,:)) / 2.d0 ) &
          * rho_0(s,:) / sqrt(1.d0-velocity_sq(s,:))
    call d01gaf( mu, Int, MDIV, D_m_0(s), er2, ifail)

    Int = exp( 2.d0*alpha(s,:) + (gama(s,:) - rho(s,:)) / 2.d0 ) &
          * (energy(s,:) + pressure(s,:)) / sqrt(1.d0-velocity_sq(s,:))
    call d01gaf( mu, Int, MDIV, D_m_p(s), er2, ifail)

    Int = sqrt(1.d0-mu(:)**2) * exp(2.d0*alpha(s,:)+gama(s,:)-rho(s,:)) &
          * (energy(s,:) + pressure(s,:)) * sqrt(velocity_sq(s,:)) / (1.d0-velocity_sq(s,:))
    call d01gaf( mu, Int, MDIV, D_J(s), er2, ifail)

    Int = Int * omg(s,:)
    call d01gaf( mu, Int, MDIV, D_T(s), er2, ifail)
  enddo

  Int = s_gp(:)**2/(1.d0-s_gp(:))**4 * D_m(:)
  call d01gaf( s_gp, Int, SDIV, Mass, er2, ifail)

  Int = s_gp(:)**2/(1.d0-s_gp(:))**4 * D_m_0(:)
  call d01gaf( s_gp, Int, SDIV, Mass_0, er2, ifail)

  Int = s_gp(:)**2/(1.d0-s_gp(:))**4 * D_m_p(:)
  call d01gaf( s_gp, Int, SDIV, Mass_p, er2, ifail)

  Int = s_gp(:)**3/(1.d0-s_gp(:))**5 * D_J(:)
  call d01gaf( s_gp, Int, SDIV, J, er2, ifail)

  Int = s_gp(:)**3/(1.d0-s_gp(:))**5 * D_T(:)
  call d01gaf( s_gp, Int, SDIV, T_kin, er2, ifail)
  
  Mass   = Mass   * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3/G
  Mass_0 = Mass_0 * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3/G
  Mass_p = Mass_p * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3/G

  if (r_ratio.eq.1.0) then
    J = 0.0
  else
    J = J * 4.d0 * pi * kappa * C**3 * r_e**4 / G
    T_kin = T_kin * 2.d0 * pi * sqrt(kappa) * C**2 * r_e**4 / G
  endif
  ang_mom = J*C/(G*MSUN**2)
  chi = J*C/(G*Mass**2)
  



  ! Compute the velocities of co-rotating and counter-rotating particles
  ! with respect to a ZAMO
  do s = 1+(SDIV-1)/2, SDIV
    s1 = s_gp(s)*(1-s_gp(s))
    s_1 = 1-s_gp(s)

    d_rho_s = deriv_s(rho,s,1)
    d_gama_s = deriv_s(gama,s,1)
    d_ww_s = deriv_s(ww,s,1)

    sqrt_v = exp(-2.0*rho(s,1))*r_e**2*s_gp(s)**4*d_ww_s**2 &
             + 2*s1*(d_gama_s+d_rho_s)+s1**2*(d_gama_s*d_gama_s-d_rho_s*d_rho_s)
    if (sqrt_v > 0.0) then
      sqrt_v = sqrt(sqrt_v)
    else
      sqrt_v = 0.0
    endif

    v_plus(s) = (exp(-rho(s,1))*r_e*s_gp(s)**2*d_ww_s + sqrt_v)/ &
              (2.0+s1*(d_gama_s-d_rho_s))
    v_minus(s) = (exp(-rho(s,1))*r_e*s_gp(s)**2*d_ww_s - sqrt_v)/ &
              (2.0+s1*(d_gama_s-d_rho_s))
  enddo

  ! Kepler angular velocity
  do s=1,SDIV
    d_r_e(s) = deriv_s(rho,s,1)
    d_g_e(s) = deriv_s(gama,s,1)
    d_o_e(s) = deriv_s(ww,s,1)
    d_v_e(s) = deriv_s( sqrt(velocity_sq) ,s,1)
    ww_mu_0(s) = ww(s,1)
  enddo

  call interp(s_gp, d_o_e, SDIV, s_e, doe)
  call interp(s_gp, d_g_e, SDIV, s_e, dge)
  call interp(s_gp, d_r_e, SDIV, s_e, dre)
  call interp(s_gp, d_v_e, SDIV, s_e, dve)

  vek = (doe/(8.0+dge-dre))*r_e*exp(-rho_equator) + sqrt(((dge+dre)/(8.0+dge &
        -dre)) + ((doe/(8.0+dge-dre))*r_e*exp(-rho_equator))**2  )

  if (r_ratio.eq.1.0) then
    ww_equator = 0.0
  else
    call interp(s_gp, ww_mu_0, SDIV, s_e, ww_equator);
  endif
  Omega_K = (C/sqrt(kappa)) * (ww_equator+vek*exp(rho_equator)/r_e)
  !write(*,"(A10,f11.5,A7,f11.5,A4)") "Keplr :",Omega_K/2/pi,"Hz = ",Omega_K,"Hz"
  
  ! Circumferential radius
  r_circ = sqrt(kappa)*r_e*exp((gama_equator-rho_equator)/2.0)
  

end  subroutine mass_radius
