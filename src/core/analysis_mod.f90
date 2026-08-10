module analysis_mod
contains

subroutine mass_radius()
  use para_mod, only: wp, SDIV, MDIV, s_gp, s_pwr, s_e, mu, &
                      gama, rho, alpha, ww, omg, sphi, energy, pressure, velocity_sq, &
                      r_e, r_ratio, r_circ, &
                      mass, mass_0, mass_s, mass_p, T_kin, ang_mom, chi, Omega_K, &
                      B_coup, mphi_r, pi, KAPPA, C, G, MSUN, &
                      v_plus, v_minus, V_rr_p, V_rr_m, output
  use toolkit_mod, only: interp, deriv_s_1d, integrate_profiles
  implicit none
  integer :: s, j
  integer, dimension(MDIV) :: s_peak_idx
  real(wp) :: s1
  real(wp), dimension(MDIV,6) :: mu_integrand_buffer
  real(wp), dimension(SDIV) :: mass_weight, ang_weight
  real(wp), dimension(SDIV) :: d_m, d_m0, d_ms, d_mp, d_j, d_t
  real(wp), dimension(SDIV,MDIV) :: rho_0
  real(wp), dimension(SDIV,6) :: integrand_buffer
  real(wp), dimension(6) :: integral_results, mu_results
  real(wp), dimension(MDIV) :: sqr1mmu2
  real(wp) :: j_local
  logical :: use_scalar ! local snapshot of the flag
  real(wp), allocatable :: alpha_t(:,:), energy_t(:,:), gama_t(:,:), omg_t(:,:), pressure_t(:,:)
  real(wp), allocatable :: rho_0_t(:,:), rho_t(:,:), sphi_t(:,:), velocity_sq_t(:,:), ww_t(:,:)

  real(wp), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  real(wp), dimension(SDIV) :: d_r_e, d_g_e, d_o_e
  real(wp), dimension(SDIV) :: dd_r_e, dd_g_e, dd_o_e
  real(wp) :: l_minus, e_minus

  call prepare_common_data(rho_0, gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0, use_scalar)
  Mass   = 0.0_wp
  mass_0 = 0.0_wp
  mass_s = 0.0_wp
  mass_p = 0.0_wp
  T_kin  = 0.0_wp
  j_local = 0.0_wp
  if (.not. use_scalar) sphi = 0.0_wp
  sqr1mmu2 = sqrt(max(0.0_wp, 1.0_wp - mu**2))

  allocate(alpha_t(MDIV,SDIV), energy_t(MDIV,SDIV), gama_t(MDIV,SDIV), omg_t(MDIV,SDIV), &
           pressure_t(MDIV,SDIV), rho_0_t(MDIV,SDIV), rho_t(MDIV,SDIV), sphi_t(MDIV,SDIV), &
           velocity_sq_t(MDIV,SDIV), ww_t(MDIV,SDIV))
  alpha_t = transpose(alpha)
  energy_t = transpose(energy)
  gama_t = transpose(gama)
  omg_t = transpose(omg)
  pressure_t = transpose(pressure)
  rho_0_t = transpose(rho_0)
  rho_t = transpose(rho)
  sphi_t = transpose(sphi)
  velocity_sq_t = transpose(velocity_sq)
  ww_t = transpose(ww)

  ! Find radial index where sphi peaks for each angular direction
  do j = 1, MDIV
    s_peak_idx(j) = maxloc(sphi(:,j), dim=1)
  enddo

  do s = 1, SDIV
    block
      real(wp), dimension(MDIV) :: acoup, acoup3, acoup4, vphi, vel_safe, e2ag, e_mr, ehgr, invlf, sqrlf, mask_sphi
      acoup = exp(-sphi_t(:,s)**2 * B_coup / 4.0_wp)
      acoup3 = acoup**3
      acoup4 = acoup**4
      vphi  = mphi_r * sphi_t(:,s)**2 / 2.0_wp
      s1 = (s_gp(s)/(1.0_wp-s_gp(s)))**s_pwr
      vel_safe = min(max(velocity_sq_t(:,s), 0.0_wp), 1.0_wp - 1.e-12_wp)

      e2ag = exp(2.0_wp*alpha_t(:,s)+gama_t(:,s))
      e_mr = exp(-rho_t(:,s))
      ehgr = exp(2.0_wp*alpha_t(:,s) + (gama_t(:,s) - rho_t(:,s))/2.0_wp)
      invlf = 1.0_wp / (1.0_wp - vel_safe)
      sqrlf = sqrt(vel_safe)

      ! Conditional mask: only compute baryon mass where:
      ! 1) sphi < 1e-3 AND
      ! 2) s is inside the radius where sphi peaks (for each mu direction)
      mask_sphi = merge(1.0_wp, 0.0_wp, &
                       (abs(sphi_t(:,s)) < 1.0e-4_wp) .and. (s <= s_peak_idx(:)))

      mu_integrand_buffer(:,1) = e2ag * &
                                ( ( ( energy_t(:,s) + pressure_t(:,s) ) * acoup4 * invlf ) * &
                                  ( 1.0_wp + vel_safe + 2.0_wp * sqrlf * s1 * sqr1mmu2 * r_e * ww_t(:,s) * e_mr ) &
                                  + 2.0_wp*pressure_t(:,s) * acoup4 - vphi / (2.0_wp*pi) )
      mu_integrand_buffer(:,2) = ehgr * rho_0_t(:,s) * acoup3 / sqrt(1.0_wp-vel_safe)
      mu_integrand_buffer(:,3) = ehgr * energy_t(:,s) * acoup4 / sqrt(1.0_wp-vel_safe)
      mu_integrand_buffer(:,4) = sqr1mmu2 * exp(2.0_wp*alpha_t(:,s) + gama_t(:,s) - rho_t(:,s)) * &
                                (energy_t(:,s) + pressure_t(:,s)) * acoup4 * sqrlf * invlf
      mu_integrand_buffer(:,5) = mu_integrand_buffer(:,4) * omg_t(:,s)
      mu_integrand_buffer(:,6) = mu_integrand_buffer(:,2) * mask_sphi
      call integrate_profiles(mu, mu_integrand_buffer, mu_results)
      d_m(s)  = mu_results(1)
      d_m0(s) = mu_results(2)
      d_ms(s) = mu_results(6)
      d_mp(s) = mu_results(3)
      d_j(s)  = mu_results(4)
      d_t(s)  = mu_results(5)
    end block
  enddo

  mass_weight = (s_gp(:)/(1.0_wp-s_gp(:)))**(3*s_pwr-1) / (1.0_wp-s_gp(:))**2 * real(s_pwr, wp)
  ang_weight  = (s_gp(:)/(1.0_wp-s_gp(:)))**(4*s_pwr-1) / (1.0_wp-s_gp(:))**2 * real(s_pwr, wp)

  integrand_buffer(:,1) = mass_weight * d_m(:)
  integrand_buffer(:,2) = mass_weight * d_m0(:)
  integrand_buffer(:,3) = mass_weight * d_ms(:)
  integrand_buffer(:,4) = mass_weight * d_mp(:)
  integrand_buffer(:,5) = ang_weight  * d_j(:)
  integrand_buffer(:,6) = ang_weight  * d_t(:)

  call integrate_profiles(s_gp, integrand_buffer, integral_results)

  Mass    = integral_results(1) * 4.0_wp * pi * sqrt(kappa)*C**2 * r_e**3 / G
  mass_0  = integral_results(2) * 4.0_wp * pi * sqrt(kappa)*C**2 * r_e**3 / G
  mass_s  = integral_results(3) * 4.0_wp * pi * sqrt(kappa)*C**2 * r_e**3 / G
  mass_p  = integral_results(4) * 4.0_wp * pi * sqrt(kappa)*C**2 * r_e**3 / G
  j_local = integral_results(5) * 4.0_wp * pi * kappa * C**3 * r_e**4 / G
  T_kin   = integral_results(6) * 2.0_wp * pi * sqrt(kappa) * C**2 * r_e**4 / G
  Omega_K = Kepler()

  ang_mom = j_local * C / (G*MSUN**2)
  chi     = j_local * C / (G*Mass**2)
contains

  real(wp) function Kepler() result(val)
    integer :: s 
    real(wp) :: doe, dge, dre, vek
    real(wp) :: s_p, gama_pole, rho_pole, gama_equator, rho_equator, sphi_equator, wwe
    s_p = r_ratio**(1.0_wp/real(s_pwr, wp)) / (1.0_wp + r_ratio**(1.0_wp/real(s_pwr, wp)))
    do s = 1, SDIV
      d_r_e(s) = deriv_s_1d(rho(:,1),s)
      d_g_e(s) = deriv_s_1d(gama(:,1),s)
      d_o_e(s) = deriv_s_1d(ww(:,1) ,s)
    end do
    call interp(s_gp, gama_mu_1,  SDIV, s_p, gama_pole)
    call interp(s_gp,  rho_mu_1,  SDIV, s_p, rho_pole)
    call interp(s_gp, gama_mu_0,  SDIV, s_e, gama_equator)
    call interp(s_gp,  rho_mu_0,  SDIV, s_e, rho_equator)
    call interp(s_gp, sphi_mu_0,  SDIV, s_e, sphi_equator)
    if (r_ratio >= 1.0_wp) then
      wwe = 0.0_wp
    else
      call interp(s_gp, ww_mu_0, SDIV, s_e, wwe)
    end if
    call interp(s_gp, d_o_e, SDIV, s_e, doe)
    call interp(s_gp, d_g_e, SDIV, s_e, dge)
    call interp(s_gp, d_r_e, SDIV, s_e, dre)
    vek = ( doe / ( 8.0_wp + dge - dre ) ) * r_e * exp(-rho_equator) + sqrt( ( (dge+dre) / (8.0_wp+dge-dre) ) &
        +(( doe / ( 8.0_wp + dge - dre ) ) * r_e * exp(-rho_equator) )**2 )
    val = (C/sqrt(kappa)) * (wwe + vek*exp(rho_equator)/r_e)
    r_circ  = sqrt(kappa) * r_e * exp((gama_equator-rho_equator)/2.0_wp) * exp(-sphi_equator**2 * B_coup / 4.0_wp)

    if (output) call write_velocity_table
  end function Kepler

  subroutine write_velocity_table
    integer :: s, unit, ios
    real(wp) :: sqrt_term, r_h
    open(newunit=unit,file="./Cont/velocity.dat",status='unknown',action='write',iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening file Cont/velocity.dat, IOSTAT=', ios
      return
    end if
    do s = 1, SDIV
      s1  = s_gp(s) * (1.0_wp - s_gp(s))
      r_h = r_e * s_gp(s) / (1.0_wp - s_gp(s))
      dd_r_e(s) = deriv_s_1d(d_r_e, s)
      dd_g_e(s) = deriv_s_1d(d_g_e, s)
      dd_o_e(s) = deriv_s_1d(d_o_e, s)
      sqrt_term = exp(-2.0_wp*rho(s,1))*r_e**2*s_gp(s)**4*d_o_e(s)**2 + &
                  2.0_wp*s1*(d_g_e(s)+d_r_e(s)) + s1**2*(d_g_e(s)**2 - d_r_e(s)**2)
      sqrt_term = merge(sqrt(sqrt_term), 0.0_wp, sqrt_term > 0.0_wp)
      v_plus(s)  = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_o_e(s) + sqrt_term ) / ( 2.0_wp + s1 * (d_g_e(s) - d_r_e(s)) )
      v_minus(s) = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_o_e(s) - sqrt_term ) / ( 2.0_wp + s1 * (d_g_e(s) - d_r_e(s)) )
      call compute_velocity_curvature(s, r_h)
      write(unit,"(10es18.9)") s_gp(s), v_plus(s), V_rr_p(s), v_minus(s), V_rr_m(s), l_minus*sqrt(KAPPA)/1.e5_wp, e_minus
    end do
    close(unit)
  end subroutine write_velocity_table

  subroutine compute_velocity_curvature(s, r_h)
    integer, intent(in) :: s
    real(wp), intent(in) :: r_h
    real(wp) :: o_r, g_r, r_r, o_rr, g_rr, r_rr, term1, term2, term3
    real(wp) :: r_h_inv, s_gp_s, r_e_inv_rh_cubed, r_e2_s_gp4_rh4

    r_h_inv = 1.0_wp / r_h
    s_gp_s = s_gp(s)
    r_e_inv_rh_cubed = -2.0_wp * r_e / (r_e + r_h)**3
    r_e2_s_gp4_rh4 = (r_e * s_gp_s**2 * r_h_inv**2)**2

    o_r  = r_e * s_gp_s**2 * r_h_inv**2 * d_o_e(s)
    g_r  = r_e * s_gp_s**2 * r_h_inv**2 * d_g_e(s)
    r_r  = r_e * s_gp_s**2 * r_h_inv**2 * d_r_e(s)
    o_rr = r_e_inv_rh_cubed * d_o_e(s) + r_e2_s_gp4_rh4 * dd_o_e(s)
    g_rr = r_e_inv_rh_cubed * d_g_e(s) + r_e2_s_gp4_rh4 * dd_g_e(s)
    r_rr = r_e_inv_rh_cubed * d_r_e(s) + r_e2_s_gp4_rh4 * dd_r_e(s)

    term1 = r_h**2 * (g_rr - r_rr + g_r**2 - r_r**2) + 2.0_wp * ( exp(-rho_mu_0(s)) * r_h**2 * o_r )**2 + 4.0_wp * r_h * r_r - 6.0_wp
    term2 = r_h**3 * exp(-rho_mu_0(s)) * ( 2.0_wp * r_r * o_r - o_rr )
    term3 = r_h**2 * ( g_rr + r_rr + g_r**2 - r_r**2 )

    V_rr_p(s) = (v_plus(s)**2 * term1 + 2.0_wp * v_plus(s) * term2 - term3) * r_h_inv**2 / (1.0_wp - v_plus(s)**2) * exp(gama_mu_0(s))
    V_rr_m(s) = (v_minus(s)**2 * term1 + 2.0_wp * v_minus(s) * term2 - term3) * r_h_inv**2 / (1.0_wp - v_minus(s)**2) * exp(gama_mu_0(s))
    l_minus = v_minus(s) * r_h * exp( (gama_mu_0(s) + rho_mu_0(s))/2.0_wp ) / sqrt(1.0_wp - v_minus(s)**2)
    e_minus = exp( (gama_mu_0(s) + rho_mu_0(s))/2.0_wp ) / sqrt(1.0_wp - v_minus(s)**2) + ww_mu_0(s) * l_minus
  end subroutine compute_velocity_curvature
end subroutine mass_radius

subroutine solution_properties()
  use eos_mod, only: p_at_e_dual, pressure_derivative_n
  use para_mod, only: wp, SDIV, MDIV, s_gp, s_e, &
                      gama, rho, alpha, ww, omg, sphi, enthalpy, &
                      energy, pressure, sound_speed, velocity_sq, &
                      r_e, r_circ, mass, MSUN, l_uni, I_inertia, Love2, &
                      pi, mphi_r, B_goal, mphi_goal, h_center, sphi_m, &
                      B_coup, KAPPA, C, n_sat, KSCALE, output
  use cheb_mod, only: cheb_diff_matrix, cheb_std_base, cheb_get_val_point
  use miscellaneous_mod, only: write_eq_profile
  use exporter_mod, only: initial_data_for_spec
  use toolkit_mod, only: interp, interp_dual, deriv_s_1d, integrate_profiles
  use ad_mod, only: dual, dual_var
  implicit none
  integer :: s, ifail, susceptibility_idx_hint
  real(wp), dimension(SDIV) :: d_r_e, d_g_e
  real(wp), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  real(wp), dimension(SDIV) :: sound_slope, sphi_deriv, pres_deriv
  real(wp), dimension(SDIV) :: effective_pressure, effective_energy, scal_op
  real(wp), dimension(SDIV) :: pressure_slope, energy_slope, effective_cs, traceT
  real(wp), dimension(SDIV) :: susceptibility, suscep_slope
  real(wp), dimension(SDIV,MDIV) :: rho_0
  character(128) :: profile_file, mphi_str, B_str, M_str, sdiv_str, mdiv_str
  real(wp) :: moi_love(2), cc, yy, dom
  logical :: use_scalar ! local snapshot of the flag
  real(wp), dimension(SDIV)  :: gamj, schwarz, gg, BV
  real(wp) :: delt = 0.005
  type(dual) :: energy_dual, pressure_dual, sphi_dual

  call prepare_common_data(rho_0, gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0, use_scalar)

  traceT = - energy(:,1) + 3.0_wp * pressure(:,1)

  if (use_scalar) then
    do s = 1, SDIV
      call interp_dual(s_gp, sphi(:,1), SDIV, dual_var(s_gp(s)), sphi_dual)
      sphi_deriv(s) = sphi_dual%der
    end do
    effective_pressure = sphi_deriv**2 / (8.0_wp * pi) - mphi_r * sphi_mu_0**2 / (4.0_wp * pi)
    effective_energy   = sphi_deriv**2 / (8.0_wp * pi) + mphi_r * sphi_mu_0**2 / (4.0_wp * pi)
    scal_op = ( 2.0_wp * pi * B_goal * traceT * exp( - B_goal * sphi_mu_0**2 / 2.0_wp ) + mphi_r ) &
            * exp( 2.0_wp * alpha(:,1) )
  else
    effective_pressure = 0.0_wp; effective_energy = 0.0_wp; scal_op = 0.0_wp; sphi_deriv = 0.0_wp
  end if

  susceptibility_idx_hint = 1
  do s = 1, SDIV
    if (energy(s,1) > 0.0_wp) then
      energy_dual = dual_var(energy(s,1))
      pressure_dual = p_at_e_dual(energy_dual)
      pres_deriv(s) = pressure_dual%der
      sound_speed(s) = pressure_dual%der
      pressure_slope(s) = deriv_s_1d(effective_pressure, s)
      energy_slope(s)   = deriv_s_1d(energy(:,1),        s)
      effective_cs(s)   = pressure_slope(s) / energy_slope(s)
      call pressure_derivative_n(energy(s,1), 2, susceptibility(s), ifail, susceptibility_idx_hint)
    else
      sound_speed(s) = 0.0_wp
      pressure_slope(s) = 0.0_wp
      energy_slope(s) = 0.0_wp
    end if
  end do
  do s = 1, SDIV
    sound_slope (s) = deriv_s_1d(sound_speed, s)
    suscep_slope(s) = deriv_s_1d(susceptibility, s)
  enddo

  if (output) call radial_configuration()
  !call to_alexis()
  !call to_sizeng()
  call mass_radius()
  moi_love = moment_inertia()
  I_inertia = moi_love(1) / (mass/MSUN*l_uni)**3

  cc = (mass/MSUN*l_uni) / (r_circ/1.e5_wp)
  yy = moi_love(2)
  dom = 2.0_wp * cc * (6.0_wp - 3.0_wp * yy + 3.0_wp * cc * (5.0_wp * yy - 8.0_wp)) &
      + 4.0_wp * cc**3 * (13.0_wp - 11.0_wp * yy + cc * (3.0_wp * yy - 2.0_wp) + 2.0_wp * cc**2 * (1.0_wp + yy)) &
      + 3.0_wp * (1.0_wp - 2.0_wp * cc)**2 * (2.0_wp * cc * (yy - 1.0_wp) - yy + 2.0_wp) * log(1.0_wp - 2.0_wp * cc)
  Love2 = 8.0_wp / 5.0_wp * cc**5 * (1.0_wp - 2.0_wp * cc)**2 * (2.0_wp * cc * (yy - 1.0_wp) - yy + 2.0_wp) &
        / dom * 2.0_wp / cc**(2*2+1) / real(2*2-1, wp)
contains
  subroutine radial_configuration()
    use para_mod, only: mass_0, eos_file
    integer :: i = 2
    select case(i)
    case(1)
        write(profile_file, '(A, A, "_", I0, "_B", ES0.2, "_mphi", ES0.2,"_M",F0.2,"_Mb",F0.4,".dat")') &
          "/Users/horay/Data4Projects/crazy/Map/muses/1dprofile_", &
          trim(eos_file), SDIV, B_goal, mphi_goal, mass / MSUN, mass_0 / MSUN
    case(2)
        profile_file = "/Users/horay/ptmp/GRASS/Res/1dprofile.dat"
    end select
    call write_eq_profile(profile_file,(SDIV-1)/2, &
          gama(:,1), rho(:,1), alpha(:,1),         &
          ww(:,1), omg(:,1),                       &
          enthalpy(:,1),                           & ! 7
          rho_0(:,1)/(KSCALE*C**2) / n_sat,        &
          energy(:,1)/(C*C*KSCALE),                &
          pressure(:,1)/KSCALE,                    &
          sphi(:,1) * sqrt(B_coup),                & ! 11
          sphi_deriv * sqrt(B_coup),               &
          sound_speed,                             & ! 13
          effective_cs,                            &
          effective_pressure/KSCALE,               & 
          effective_energy/(C*C*KSCALE)  )

    if (.false.) then ! Debug: Chebyshev fit of gama over the stellar interior [s_gp(1), s_gp(res+1)]
      block
        integer, parameter :: n_cheb = 39
        integer :: res, i
        real(wp) :: coeffs(0:n_cheb)
        real(wp), allocatable :: Y_cheb(:), Y(:)
        res = (SDIV-1)/2
        allocate(Y_cheb(res+1), Y(res+1))
        Y = alpha(1:res+1,1)
        call cheb_std_base(n_cheb, res+1, s_gp(1:res+1), Y, s_gp(1), s_gp(res+1), coeffs)
        do i = 1, res+1
          Y_cheb(i) = cheb_get_val_point(n_cheb, coeffs, s_gp(1), s_gp(res+1), s_gp(i))
        end do
        write(*,'(A,2es12.4)') 'cheb err (abs, rel):', &
            maxval(abs(Y_cheb - Y)), &
            maxval(abs(Y_cheb/Y - 1.0_wp))
        block
          integer :: uid
          open(newunit=uid, file='./Cont/cheb_check.dat', status='replace', action='write')
          do i = 1, res+1
            write(uid,'(3es25.16)') s_gp(i), Y_cheb(i), Y(i)
          end do
          close(uid)
        end block
        deallocate(Y_cheb)
      end block
    end if
  end subroutine radial_configuration

  subroutine to_alexis()
    profile_file = "./Cont/alexis.dat"
    gamj  = sound_speed * (energy(:,1)/pressure(:,1) - 1.0_wp) / (1.0_wp+delt)
    schwarz  = pres_deriv / pressure(:,1)/ gamj * delt / (1.0_wp+delt)
    gg    = -exp(2.0_wp*rho(:,1)) * pres_deriv / ( energy(:,1) + pressure(:,1) )
    BV    = sqrt( - gg * schwarz ) / r_e * ( 1.0_wp - s_gp )**2 * sqrt(KAPPA) / C
    call write_eq_profile(profile_file, (SDIV+1)/2, &
        gama(:,1), rho(:,1), alpha(:,1), & ! 1-3
        enthalpy(:,1), rho_0(:,1)/(KSCALE*C**2), pressure(:,1)/KSCALE, & ! 4-6
        sound_speed, BV) ! 7-8
  end subroutine to_alexis  

  subroutine to_sizeng()
    real(wp), parameter :: RHO_TO_KM = 1.e12_wp * 6.67408e-20_wp / (2.99792458e5_wp)**2
    real(wp), parameter :: K_KM = 218.04217865726338_wp
    profile_file = "./Cont/sizeng.dat"
    call initial_data_for_spec( profile_file, rho_0 / (KSCALE*C**2) * rho_to_km * K_km, &
        alpha, rho, gama, ww * ( sqrt(K_km) / sqrt(KAPPA) ), sqrt(velocity_sq) )
  end subroutine to_sizeng  

end subroutine solution_properties

subroutine prepare_common_data(rho_0, gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0, use_scalar)
  use eos_mod, only: n0_at_e
  use para_mod, only: wp, SDIV, MDIV, has_scalar, &
                      gama, rho, ww, sphi, energy, e_surface, MB, KSCALE, C
  implicit none
  real(wp), dimension(SDIV,MDIV), intent(out) :: rho_0
  real(wp), dimension(SDIV), intent(out) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  logical, intent(out) :: use_scalar

  use_scalar = has_scalar

  ! Vectorized extraction of boundary values
  gama_mu_1(:) = gama(:,MDIV)
  rho_mu_1(:)  = rho (:,MDIV)
  gama_mu_0(:) = gama(:,1)
  rho_mu_0(:)  = rho (:,1)
  ww_mu_0 (:)  = ww  (:,1)
  sphi_mu_0(:) = merge( sphi(:,1), 0.0_wp, use_scalar )

  block
    integer :: s_, m_
    real(wp) :: n0_val
    do m_ = 1, MDIV; do s_ = 1, SDIV
        if (energy(s_,m_) > e_surface) then
          n0_val = n0_at_e(energy(s_,m_))
          rho_0(s_,m_) = n0_val * MB * KSCALE * C**2
        else
          rho_0(s_,m_) = 0.0_wp
        end if
    end do; end do
  end block
end subroutine prepare_common_data
  
function moment_inertia() result(val)
  use iso_fortran_env, only: error_unit
  use para_mod, only: wp, SDIV, s_gp, r_e, KAPPA, energy, pressure, sound_speed, &
                      C, KSCALE, rho_uni, prs_uni, pi
  use nag_compat_mod, only: d02pcf
  use toolkit_mod, only: interp
  implicit none
  real(wp) :: r_in, r_surf, r_start ! in km
  integer, parameter :: neqn = 4
  real(wp) :: relerr = 1.e-8_wp, abserr = 1.e-8_wp
  real(wp) :: val(2), ec, pc, y(neqn), yp(neqn)
  real(wp) :: s_fail, e_fail, p_fail, vs2_fail
  integer :: flag, step_count
  logical :: debug = .false.

  r_surf = r_e * sqrt(KAPPA) / 1.e5_wp
  r_in = r_surf * s_gp(2) / ( 1.0_wp - s_gp(2) )
  r_start = r_in
  ec = energy(1,1) / (C*C*KSCALE) * rho_uni
  pc = pressure(1,1) / KSCALE * prs_uni

  y(1) = r_in - pi * ec * r_in**3
  y(2) = 8.0_wp * pi / 15.0_wp * y(1)**5 * (ec+pc)
  y(3) = 4.0_wp * pi * ec * y(1)**2 * (1.0_wp + 4.0_wp * pi * ec * y(1)**2 / 3.0_wp)
  y(4) = 2.0_wp
  yp(1)= 1 - 3.0_wp * pi * ec * r_in**2
  yp(2)= 8.0_wp * pi / 3.0_wp * y(1)**4 * (ec+pc)
  yp(3)= 8.0_wp * pi * ec * y(1) * (1.0_wp + 2.0_wp * pi * ec * y(1)**2)
  yp(4)= 0.0_wp

  flag = 1
  call d02pcf(deriv, neqn, y, yp, r_in, r_surf, relerr, abserr, flag, step_count, debug)
  if (abs(flag) /= 2) then
    write(error_unit,'(A,I0)') "d02pcf failed with flag = ", flag
    select case (abs(flag))
    case (4)
      write(error_unit,'(A)') "  reason    : maximum step count reached"
    case (6)
      write(error_unit,'(A)') "  reason    : step size fell below HMIN"
    case (8)
      write(error_unit,'(A)') "  reason    : bad input"
    case default
      write(error_unit,'(A)') "  reason    : unknown failure"
    end select
    write(error_unit,'(A,1X,ES22.14)') "  t_start   :", r_start
    write(error_unit,'(A,1X,ES22.14)') "  t_fail    :", r_in
    write(error_unit,'(A,1X,ES22.14)') "  t_out     :", r_surf
    write(error_unit,'(A,1X,I0)')      "  steps     :", step_count
    write(error_unit,'(A,1X,ES22.14)') "  ec        :", ec
    write(error_unit,'(A,1X,ES22.14)') "  pc        :", pc
    write(error_unit,'(A,1X,4(ES22.14,1X))') "  y         :", y
    write(error_unit,'(A,1X,4(ES22.14,1X))') "  yp        :", yp
    s_fail = r_in / (r_in + r_surf)
    call interp(s_gp, energy(:,1),   SDIV, s_fail, e_fail)
    call interp(s_gp, pressure(:,1), SDIV, s_fail, p_fail)
    call interp(s_gp, sound_speed,   SDIV, s_fail, vs2_fail)
    e_fail = e_fail / (C*C*KSCALE) * rho_uni
    p_fail = p_fail / KSCALE * prs_uni
    write(error_unit,'(A,1X,ES22.14)') "  s_fail    :", s_fail
    write(error_unit,'(A,1X,ES22.14)') "  e_local   :", e_fail
    write(error_unit,'(A,1X,ES22.14)') "  p_local   :", p_fail
    write(error_unit,'(A,1X,ES22.14)') "  vs2_local :", vs2_fail
    write(error_unit,'(A,1X,ES22.14)') "  (e+p)/vs2 :", (e_fail + p_fail) / vs2_fail
    error stop "moment_inertia: d02pcf failed"
  end if
  val(1) = y(2)
  val(2) = y(4)
  !write(*,"('Double check Schwarzschild radius:',es18.9,'  ADM mass:',es18.9, '  RK45 steps :', i5)") &
  !  abs(1.0_wp-r_circ/1.e5_wp/y(1)), abs(1.0_wp-y(3)/l_uni/(mass/MSUN)), step_count

end function moment_inertia

subroutine deriv(t, y, yp)
  use precision_mod, only: wp
  use para_mod, only: SDIV, KAPPA, r_e, s_gp, &
                      gama, rho, alpha, sphi, energy, pressure, sound_speed, &
                      C, KSCALE, rho_uni, prs_uni, pi, has_scalar
  use ad_mod, only: dual, dual_var
  use toolkit_mod, only: interp, interp_dual
  implicit none
  real(wp), intent(in) :: t, y(:)
  real(wp), intent(out) :: yp(:)
  real(wp) :: gama_val, rho_val, s_h, e, p, elm, dpdr, QQ
  real(wp) :: dalphads, logP, dsphids, vs2, r_surf
  type(dual) :: s_d, alpha_d, sphi_d

  r_surf = r_e * sqrt(KAPPA) / 1.e5_wp
  s_h = t / ( t + r_surf )
  call interp(s_gp, gama(:,1)     , SDIV, s_h, gama_val)
  call interp(s_gp, rho(:,1)      , SDIV, s_h,  rho_val)
  call interp(s_gp, energy(:,1)   , SDIV, s_h,        e)
  call interp(s_gp, pressure(:,1) , SDIV, s_h,        p)
  call interp(s_gp, sound_speed   , SDIV, s_h,      vs2)
  s_d = dual_var(s_h)
  call interp_dual(s_gp, alpha(:,1), SDIV, s_d, alpha_d)
  call interp_dual(s_gp, sphi(:,1),  SDIV, s_d,  sphi_d)
  dalphads = alpha_d%der
  dsphids  = sphi_d%der

  logP = ( 4.0_wp * alpha_d%val + gama_val - rho_val ) / 12.0_wp

  e = e / (C*C*KSCALE) * rho_uni
  p = p / KSCALE * prs_uni

  elm = 1.0_wp / ( 1.0_wp + s_h * (1.0_wp - s_h) * dalphads )**2
  dpdr= - (e+p) * (y(3) + 4.0_wp * pi * y(1)**3 * p) / (y(1) * (y(1) - 2.0_wp * y(3)))

  yp(1) = exp(2.0_wp * logP) * ( 1.0_wp + s_h * (1.0_wp - s_h) * dalphads )
  yp(2) = 8.0_wp / 3.0_wp * pi * y(1)**4 * (e+p) * ( 1.0_wp - 5.0_wp * y(2) / 2.0_wp / y(1)**3 + y(2)**2 / y(1)**6 ) * elm
  yp(3) = 4.0_wp * pi * y(1)**2 * e

  QQ  = -real((1+1)*(1+2), wp) * elm / y(1)**2 - dpdr**2 &
      + 4.0_wp * pi * elm * (5.0_wp * e + 9.0_wp * p + merge((e+p) / vs2, 0._wp, vs2 > 1.e-6_wp))
  yp(4) = -y(4)**2 / y(1) - y(4) * elm / y(1) * (1.0_wp + 4.0_wp * pi * y(1)**2 * (p-e)) - QQ * y(1)

  if (has_scalar) then
    yp(2) = yp(2) + y(2) * y(1) * ( 1.0_wp - 2.0_wp * y(2) / y(1)**3 ) / yp(1)**2 &
          * ( (1.0_wp - s_h)**2 / r_surf * dsphids )**2
  endif

  yp(2) = yp(2) * yp(1)
  yp(3) = yp(3) * yp(1)
  yp(4) = yp(4) * yp(1)
end subroutine deriv

end module analysis_mod
