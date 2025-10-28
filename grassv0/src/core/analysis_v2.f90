subroutine mass_radius
  use para_mod
  use miscellaneous_mod, only: write_eq_profile
  use toolkit_mod, only: interp, deriv_s, deriv_s_1d, integrate_profiles
  use ad_mod, only: dual, dual_var
  implicit none
  integer :: s, m
  real(8) :: s_p, r_p
  real(8) :: gama_pole, rho_pole, gama_equator, rho_equator, ww_equator, sphi_equator
  real(8) :: doe, dge, dre, dve, vek
  real(8) :: sqrt_term, mphi_local
  real(8) :: s1, s_1, r_h
  real(8), dimension(MDIV) :: scal, acoup, vphi, vel_safe
  real(8), dimension(MDIV,5) :: mu_integrand_buffer
  real(8), dimension(SDIV) :: mass_weight, ang_weight
  real(8), dimension(SDIV) :: d_m, d_m0, d_mp, d_j, d_t
  real(8), dimension(SDIV) :: d_r_e, d_g_e, d_o_e, d_v_e
  real(8), dimension(SDIV) :: dd_r_e, dd_g_e, dd_o_e
  real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  real(8), dimension(SDIV) :: sound_speed, sound_slope
  real(8), dimension(SDIV,MDIV) :: rho_0
  real(8), dimension(SDIV,5) :: integrand_buffer
  real(8), dimension(5) :: integral_results, mu_results
  real(8) :: l_minus, e_minus
  character(64) :: profile_file, mphi_str, B_str
  real(8) :: j_local
  logical :: use_scalar ! local snapshot of the flag
  real(8), external :: n0_at_e
  type(dual) :: energy_dual, pressure_dual

  interface
    function p_at_e_dual(ee) result(res)
      import :: dual
      implicit none
      type(dual), intent(in) :: ee
      type(dual) :: res
    end function p_at_e_dual
  end interface

  use_scalar = has_scalar
  s_p = r_ratio**(1.d0/dble(s_pwr)) / (1.d0 + r_ratio**(1.d0/dble(s_pwr)))

  do s = 1, SDIV
    gama_mu_1(s) = gama(s,MDIV)
    rho_mu_1(s)  = rho (s,MDIV)
    gama_mu_0(s) = gama(s,1)
    rho_mu_0(s)  = rho (s,1)
    ww_mu_0 (s)  = ww  (s,1)
    if (use_scalar) sphi_mu_0(s) = sphi(s,1)
    do m = 1, MDIV
      if (energy(s,m) > e_surface) then
        rho_0(s,m) = n0_at_e(energy(s,m)) * MB * KSCALE*C**2 ! * KSCALE*C**2 is necessary here to compute Mb
      else
        rho_0(s,m) = 0.d0
      end if
    end do
  end do

  call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole)
  call interp(s_gp,  rho_mu_1,  SDIV, s_p, rho_pole)
  call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator)
  call interp(s_gp,  rho_mu_0,  SDIV, s_e, rho_equator)
  call interp(s_gp, sphi(:,1),  SDIV, s_e, sphi_equator)

  if (output) then
    do s = 1, SDIV
      if (energy(s,1) > 0.d0) then
        energy_dual = dual_var(energy(s,1))
        pressure_dual = p_at_e_dual(energy_dual)
        sound_speed(s) = pressure_dual%der
      else
        sound_speed(s) = 0.d0
      end if
    end do
    do s = 1, SDIV 
      sound_slope(s) = deriv_s_1d(sound_speed, s)
    enddo

    write(mphi_str,"(f10.0)") mphi_goal
    write(B_str,"(f10.0)") B_goal
    profile_file = "./Cont/1dprofile_"//trim(adjustl(B_str))//"dat"
    call write_eq_profile(profile_file, 2*SDIV/3, &
         omg(:,1)/(2.d0*pi)*(C/sqrt(kappa)), F_j(:,1), &
         enthalpy(:,1), rho_0(:,1)/(KSCALE*C**2)/ n_sat, sound_speed, sound_slope, &
         sphi(:,1) )
  end if

  Mass   = 0.d0
  mass_0 = 0.d0
  mass_p = 0.d0
  T_kin  = 0.d0
  j_local = 0.d0

  if (use_scalar) then
    r_p = r_ratio * r_e
    s_p = r_p / (r_p + r_e)
    call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole)
    call interp(s_gp, rho_mu_1,  SDIV, s_p, rho_pole)
    mphi_local = mphi_r

    do s = 1, SDIV
      scal = sphi(s,:)
      acoup = exp(-scal**2 * B_coup / 4.d0)
      vphi  = mphi_local * scal**2 / 2.d0
      vel_safe = min(max(velocity_sq(s,:), 0.d0), 1.d0 - 1.d-12)

      mu_integrand_buffer(:,1) = exp(2.d0*alpha(s,:)+gama(s,:)) * &
                                 ( ((energy(s,:)+pressure(s,:))*acoup**4/(1.d0-vel_safe)) * &
                                   (1.d0 + vel_safe + 2.d0*s_gp(s)*sqrt(vel_safe)/(1.d0-s_gp(s)) * &
                                    sqrt(1.d0-mu(:)**2) * r_e * ww(s,:) * exp(-rho(s,:))) &
                                    + 2.d0*pressure(s,:)*acoup**4 - vphi/(2.d0*pi) )
      mu_integrand_buffer(:,2) = exp(2.d0*alpha(s,:) + (gama(s,:) - rho(s,:))/2.d0) * rho_0(s,:) * acoup**3 / sqrt(1.d0-vel_safe)
      mu_integrand_buffer(:,3) = exp(2.d0*alpha(s,:) + (gama(s,:) - rho(s,:))/2.d0) * energy(s,:)*acoup**4 / sqrt(1.d0-vel_safe)
      mu_integrand_buffer(:,4) = sqrt(1.d0-mu(:)**2) * exp(2.d0*alpha(s,:)+gama(s,:)-rho(s,:)) * &
                                 (energy(s,:) + pressure(s,:)) * acoup**4 * sqrt(vel_safe) / (1.d0-vel_safe)
      mu_integrand_buffer(:,5) = mu_integrand_buffer(:,4) * omg(s,:)

      call integrate_profiles(mu, mu_integrand_buffer, mu_results)

      d_m(s)  = mu_results(1)
      d_m0(s) = mu_results(2)
      d_mp(s) = mu_results(3)
      d_j(s)  = mu_results(4)
      d_t(s)  = mu_results(5)
    end do

    mass_weight = s_gp(:)**2/(1.d0-s_gp(:))**4
    ang_weight  = s_gp(:)**3/(1.d0-s_gp(:))**5

    integrand_buffer(:,1) = mass_weight * d_m(:)
    integrand_buffer(:,2) = mass_weight * d_m0(:)
    integrand_buffer(:,3) = mass_weight * d_mp(:)
    integrand_buffer(:,4) = ang_weight  * d_j(:)
    integrand_buffer(:,5) = ang_weight  * d_t(:)

    call integrate_profiles(s_gp, integrand_buffer, integral_results)

    Mass   = integral_results(1)
    mass_0 = integral_results(2)
    mass_p = integral_results(3)
    j_local = integral_results(4)
    T_kin   = integral_results(5)

    Mass   = Mass   * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3 / G
    mass_0 = mass_0 * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3 / G
    mass_p = mass_p * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3 / G

    if (r_ratio == 1.d0) then
      j_local = 0.d0
      T_kin   = 0.d0
    else
      j_local = j_local * 4.d0 * pi * kappa * C**3 * r_e**4 / G
      T_kin   = T_kin   * 2.d0 * pi * sqrt(kappa) * C**2 * r_e**4 / G
    end if
    
    call compute_dynamics(gama_equator, rho_equator, ww_equator, sphi_equator)

  else
    do s = 1, SDIV
      s1 = (s_gp(s)/(1.d0-s_gp(s)))**s_pwr
      vel_safe = min(max(velocity_sq(s,:), 0.d0), 1.d0 - 1.d-12)

      mu_integrand_buffer(:,1) = exp(2.d0*alpha(s,:)+gama(s,:)) * &
                                 ( ((energy(s,:)+pressure(s,:))/(1.d0-vel_safe)) * &
                                   (1.d0 + vel_safe + 2.d0*sqrt(vel_safe)*s1*sqrt(1.d0-mu(:)**2) * r_e * ww(s,:) * exp(-rho(s,:))) &
                                   + 2.d0*pressure(s,:) )
      mu_integrand_buffer(:,2) = exp(2.d0*alpha(s,:) + (gama(s,:) - rho(s,:))/2.d0) * rho_0(s,:) / sqrt(1.d0-vel_safe)
      mu_integrand_buffer(:,3) = exp(2.d0*alpha(s,:) + (gama(s,:) - rho(s,:))/2.d0) * energy(s,:) / sqrt(1.d0-vel_safe)
      mu_integrand_buffer(:,4) = sqrt(1.d0-mu(:)**2) * exp(2.d0*alpha(s,:)+gama(s,:)-rho(s,:)) * &
                                 (energy(s,:) + pressure(s,:)) * sqrt(vel_safe) / (1.d0-vel_safe)
      mu_integrand_buffer(:,5) = mu_integrand_buffer(:,4) * omg(s,:)

      call integrate_profiles(mu, mu_integrand_buffer, mu_results)

      d_m(s)  = mu_results(1)
      d_m0(s) = mu_results(2)
      d_mp(s) = mu_results(3)
      d_j(s)  = mu_results(4)
      d_t(s)  = mu_results(5)

    end do

    mass_weight = (s_gp(:)/(1.d0-s_gp(:)))**(3*s_pwr-1) / (1.d0-s_gp(:))**2 * dble(s_pwr)
    ang_weight  = (s_gp(:)/(1.d0-s_gp(:)))**(4*s_pwr-1) / (1.d0-s_gp(:))**2 * dble(s_pwr)

    integrand_buffer(:,1) = mass_weight * d_m(:)
    integrand_buffer(:,2) = mass_weight * d_m0(:)
    integrand_buffer(:,3) = mass_weight * d_mp(:)
    integrand_buffer(:,4) = ang_weight  * d_j(:)
    integrand_buffer(:,5) = ang_weight  * d_t(:)

    call integrate_profiles(s_gp, integrand_buffer, integral_results)

    Mass    = integral_results(1) * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3 / G
    mass_0  = integral_results(2) * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3 / G
    mass_p  = integral_results(3) * 4.d0*pi*sqrt(kappa)*C**2 * r_e**3 / G
    j_local = integral_results(4) * 4.d0 * pi * kappa * C**3 * r_e**4 / G
    T_kin   = integral_results(5) * 2.d0 * pi * sqrt(kappa) * C**2 * r_e**4 / G
    call compute_dynamics(gama_equator, rho_equator, ww_equator)
  end if

  ang_mom = j_local*C/(G*MSUN**2)
  chi     = j_local*C/(G*Mass**2)

contains

  subroutine compute_dynamics(ge, re, wwe, s_equator)
    real(8), intent(out) :: ge, re, wwe
    real(8), intent(out), optional :: s_equator
    integer :: s
    logical :: scalar_mode

    scalar_mode = present(s_equator)

    if (scalar_mode) then
      do s = 1, SDIV
        d_r_e(s) = deriv_s(rho,s,1)
        d_g_e(s) = deriv_s(gama,s,1)
        d_o_e(s) = deriv_s(ww,s,1)
        d_v_e(s) = deriv_s(sqrt(max(velocity_sq(s,:),0.d0)),s,1)
      end do
    else
      do s = 1, SDIV
        d_r_e(s) = deriv_s_1d(rho,s)
        d_g_e(s) = deriv_s_1d(gama,s)
        d_o_e(s) = deriv_s_1d(ww ,s)
        d_v_e(s) = deriv_s_1d(sqrt(max(velocity_sq(s,:),0.d0)),s)
      end do
    end if
    call interp(s_gp, d_o_e, SDIV, s_e, doe)
    call interp(s_gp, d_g_e, SDIV, s_e, dge)
    call interp(s_gp, d_r_e, SDIV, s_e, dre)
    call interp(s_gp, d_v_e, SDIV, s_e, dve)
    call interp(s_gp, gama_mu_0, SDIV, s_e, ge)
    call interp(s_gp, rho_mu_0,  SDIV, s_e, re)
    if (r_ratio.eq.1.d0) then
      wwe = 0.d0
    else
      call interp(s_gp, ww_mu_0, SDIV, s_e, wwe)
    end if
    vek = (doe/(8.d0+dge-dre))*r_e*exp(-re) + sqrt(((dge+dre)/(8.d0+dge-dre))+((doe/(8.d0+dge-dre))*r_e*exp(-re))**2)
    Omega_K = (C/sqrt(kappa)) * (wwe + vek*exp(re)/r_e)
    r_circ  = sqrt(kappa) * r_e * exp((ge-re)/2.d0)
    if (scalar_mode) then
      r_circ = r_circ * exp(-s_equator**2 * B_coup / 4.d0)
    else
      call write_velocity_table
    end if
  end subroutine compute_dynamics

  subroutine write_velocity_table
    integer :: s
    open(23,file="./Cont/velocity.dat")
    do s = 1, SDIV
      s1  = s_gp(s) * (1.d0 - s_gp(s))
      s_1 = 1.d0 - s_gp(s)
      r_h = r_e * s_gp(s) / (1.d0 - s_gp(s))
      dd_r_e(s) = deriv_s_1d(d_r_e,s)
      dd_g_e(s) = deriv_s_1d(d_g_e,s)
      dd_o_e(s) = deriv_s_1d(d_o_e,s)
      sqrt_term = exp(-2.d0*rho(s,1))*r_e**2*s_gp(s)**4*d_o_e(s)**2 + &
                  2.d0*s1*(d_g_e(s)+d_r_e(s)) + s1**2*(d_g_e(s)**2 - d_r_e(s)**2)
      sqrt_term = merge(sqrt(sqrt_term), 0.d0, sqrt_term > 0.d0)
      v_plus(s)  = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_o_e(s) + sqrt_term ) / ( 2.d0 + s1 * (d_g_e(s) - d_r_e(s)) )
      v_minus(s) = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_o_e(s) - sqrt_term ) / ( 2.d0 + s1 * (d_g_e(s) - d_r_e(s)) )
      call compute_velocity_curvature(s, r_h)
      write(23,"(10es18.9)") s_gp(s), v_plus(s), V_rr_p(s), v_minus(s), V_rr_m(s), l_minus*sqrt(KAPPA)/1.d5, e_minus
    end do
    close(23)
  end subroutine write_velocity_table

  subroutine compute_velocity_curvature(s, r_h)
    integer, intent(in) :: s
    real(8), intent(in) :: r_h
    real(8) :: o_r, g_r, r_r, o_rr, g_rr, r_rr
    o_r  = r_e * s_gp(s)**2 / r_h**2 * d_o_e(s)
    g_r  = r_e * s_gp(s)**2 / r_h**2 * d_g_e(s)
    r_r  = r_e * s_gp(s)**2 / r_h**2 * d_r_e(s)
    o_rr = -2.d0 * r_e / (r_e + r_h)**3 * d_o_e(s) + r_e**2 * s_gp(s)**4 / r_h**4 * dd_o_e(s)
    g_rr = -2.d0 * r_e / (r_e + r_h)**3 * d_g_e(s) + r_e**2 * s_gp(s)**4 / r_h**4 * dd_g_e(s)
    r_rr = -2.d0 * r_e / (r_e + r_h)**3 * d_r_e(s) + r_e**2 * s_gp(s)**4 / r_h**4 * dd_r_e(s)
    V_rr_p(s) = v_plus(s)**2 * ( r_h**2 * (g_rr - r_rr + g_r**2 - r_r**2) + &
                   2.d0 * ( exp(-rho_mu_0(s)) * r_h**2 * o_r )**2 + 4.d0 * r_h * r_r - 6.d0 ) + &
                 2.d0 * v_plus(s) * r_h**3 * exp(-rho_mu_0(s)) * ( 2.d0 * r_r * o_r - o_rr ) - &
                 r_h**2 * ( g_rr + r_rr + g_r**2 - r_r**2 )
    V_rr_p(s) = V_rr_p(s) / r_h**2 / (1.d0 - v_plus(s)**2) * exp(gama_mu_0(s))
    V_rr_m(s) = v_minus(s)**2 * ( r_h**2 * (g_rr - r_rr + g_r**2 - r_r**2) + &
                   2.d0 * ( exp(-rho_mu_0(s)) * r_h**2 * o_r )**2 + 4.d0 * r_h * r_r - 6.d0 ) + &
                 2.d0 * v_minus(s) * r_h**3 * exp(-rho_mu_0(s)) * ( 2.d0 * r_r * o_r - o_rr ) - &
                 r_h**2 * ( g_rr + r_rr + g_r**2 - r_r**2 )
    V_rr_m(s) = V_rr_m(s) / r_h**2 / (1.d0 - v_minus(s)**2) * exp(gama_mu_0(s))
    l_minus = v_minus(s) * r_h * exp( (gama_mu_0(s) + rho_mu_0(s))/2.d0 ) / sqrt(1.d0 - v_minus(s)**2)
    e_minus = exp( (gama_mu_0(s) + rho_mu_0(s))/2.d0 ) / sqrt(1.d0 - v_minus(s)**2) + ww_mu_0(s) * l_minus
  end subroutine compute_velocity_curvature

end subroutine mass_radius
