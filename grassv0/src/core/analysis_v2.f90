subroutine mass_radius()
  use para_mod
  use toolkit_mod, only: interp, deriv_s_1d, integrate_profiles
  implicit none
  integer :: s
  real(8) :: s1, mphi_local
  real(8), dimension(MDIV) :: scal, acoup, vphi, vel_safe
  real(8), dimension(MDIV,5) :: mu_integrand_buffer
  real(8), dimension(SDIV) :: mass_weight, ang_weight
  real(8), dimension(SDIV) :: d_m, d_m0, d_mp, d_j, d_t
  real(8), dimension(SDIV,MDIV) :: rho_0
  real(8), dimension(SDIV,5) :: integrand_buffer
  real(8), dimension(5) :: integral_results, mu_results
  real(8) :: j_local
  logical :: use_scalar ! local snapshot of the flag

  real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  real(8), dimension(SDIV) :: d_r_e, d_g_e, d_o_e
  real(8), dimension(SDIV) :: dd_r_e, dd_g_e, dd_o_e
  real(8) :: l_minus, e_minus

  call prepare_common_data(rho_0, gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0, use_scalar)
  Mass   = 0.d0
  mass_0 = 0.d0
  mass_p = 0.d0
  T_kin  = 0.d0
  j_local = 0.d0
  if (.not. use_scalar) sphi = 0.d0

  do s = 1, SDIV
    acoup = exp(-scal**2 * B_coup / 4.d0)
    vphi  = mphi_local * scal**2 / 2.d0
    s1 = (s_gp(s)/(1.d0-s_gp(s)))**s_pwr
    vel_safe = min(max(velocity_sq(s,:), 0.d0), 1.d0 - 1.d-12)

    mu_integrand_buffer(:,1) = exp(2.d0*alpha(s,:)+gama(s,:)) * &
                               ( ( ( energy(s,:) + pressure(s,:) ) * acoup**4 / (1.d0 - vel_safe) ) * &
                                 ( 1.d0 + vel_safe + 2.d0 * sqrt(vel_safe) * s1 * sqrt(1.d0-mu(:)**2) * r_e * ww(s,:) * exp(-rho(s,:)) ) &
                                 + 2.d0*pressure(s,:) * acoup**4 - vphi / (2.d0*pi) )
    mu_integrand_buffer(:,2) = exp(2.d0*alpha(s,:) + (gama(s,:) - rho(s,:))/2.d0) * rho_0(s,:) * acoup**3  / sqrt(1.d0-vel_safe)
    mu_integrand_buffer(:,3) = exp(2.d0*alpha(s,:) + (gama(s,:) - rho(s,:))/2.d0) * energy(s,:)* acoup**4  / sqrt(1.d0-vel_safe)
    mu_integrand_buffer(:,4) = sqrt(1.d0-mu(:)**2) * exp( 2.d0*alpha(s,:) + gama(s,:) - rho(s,:) ) * &
                               (energy(s,:) + pressure(s,:)) * acoup**4 * sqrt(vel_safe) / (1.d0-vel_safe)
    mu_integrand_buffer(:,5) = mu_integrand_buffer(:,4) * omg(s,:)
    call integrate_profiles(mu, mu_integrand_buffer, mu_results)
    d_m(s)  = mu_results(1)
    d_m0(s) = mu_results(2)
    d_mp(s) = mu_results(3)
    d_j(s)  = mu_results(4)
    d_t(s)  = mu_results(5)
  enddo

  mass_weight = (s_gp(:)/(1.d0-s_gp(:)))**(3*s_pwr-1) / (1.d0-s_gp(:))**2 * dble(s_pwr)
  ang_weight  = (s_gp(:)/(1.d0-s_gp(:)))**(4*s_pwr-1) / (1.d0-s_gp(:))**2 * dble(s_pwr)

  integrand_buffer(:,1) = mass_weight * d_m(:)
  integrand_buffer(:,2) = mass_weight * d_m0(:)
  integrand_buffer(:,3) = mass_weight * d_mp(:)
  integrand_buffer(:,4) = ang_weight  * d_j(:)
  integrand_buffer(:,5) = ang_weight  * d_t(:)

  call integrate_profiles(s_gp, integrand_buffer, integral_results)

  Mass    = integral_results(1) * 4.d0 * pi * sqrt(kappa)*C**2 * r_e**3 / G
  mass_0  = integral_results(2) * 4.d0 * pi * sqrt(kappa)*C**2 * r_e**3 / G
  mass_p  = integral_results(3) * 4.d0 * pi * sqrt(kappa)*C**2 * r_e**3 / G
  j_local = integral_results(4) * 4.d0 * pi * kappa * C**3 * r_e**4 / G
  T_kin   = integral_results(5) * 2.d0 * pi * sqrt(kappa) * C**2 * r_e**4 / G
  Omega_K = Kepler()

  ang_mom = j_local*C/(G*MSUN**2)
  chi     = j_local*C/(G*Mass**2)
contains

  real(8) function Kepler() result(val)
    integer :: s 
    real(8) :: doe, dge, dre, vek
    real(8) :: s_p, gama_pole, rho_pole, gama_equator, rho_equator, ww_equator, sphi_equator, wwe
    s_p = r_ratio**(1.d0/dble(s_pwr)) / (1.d0 + r_ratio**(1.d0/dble(s_pwr)))
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
    if (r_ratio.eq.1.d0) then
      wwe = 0.d0
    else
      call interp(s_gp, ww_mu_0, SDIV, s_e, wwe)
    end if
    call interp(s_gp, d_o_e, SDIV, s_e, doe)
    call interp(s_gp, d_g_e, SDIV, s_e, dge)
    call interp(s_gp, d_r_e, SDIV, s_e, dre)
    vek = ( doe / ( 8.d0 + dge - dre ) ) * r_e * exp(-rho_equator) + sqrt( ( (dge+dre) / (8.d0+dge-dre) ) &
        +(( doe / ( 8.d0 + dge - dre ) ) * r_e * exp(-rho_equator) )**2 )
    val = (C/sqrt(kappa)) * (wwe + vek*exp(rho_equator)/r_e)
    r_circ  = sqrt(kappa) * r_e * exp((gama_equator-rho_equator)/2.d0) * exp(-sphi_equator**2 * B_coup / 4.d0)
    if (output) call write_velocity_table
  end function Kepler

  subroutine write_velocity_table
    integer :: s, unit, ios
    real(8) :: sqrt_term, r_h
    open(newunit=unit,file="./Cont/velocity.dat",status='unknown',action='write',iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening file Cont/velocity.dat, IOSTAT=', ios
    end if
    do s = 1, SDIV
      s1  = s_gp(s) * (1.d0 - s_gp(s))
      r_h = r_e * s_gp(s) / (1.d0 - s_gp(s))
      dd_r_e(s) = deriv_s_1d(d_r_e, s)
      dd_g_e(s) = deriv_s_1d(d_g_e, s)
      dd_o_e(s) = deriv_s_1d(d_o_e, s)
      sqrt_term = exp(-2.d0*rho(s,1))*r_e**2*s_gp(s)**4*d_o_e(s)**2 + &
                  2.d0*s1*(d_g_e(s)+d_r_e(s)) + s1**2*(d_g_e(s)**2 - d_r_e(s)**2)
      sqrt_term = merge(sqrt(sqrt_term), 0.d0, sqrt_term > 0.d0)
      v_plus(s)  = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_o_e(s) + sqrt_term ) / ( 2.d0 + s1 * (d_g_e(s) - d_r_e(s)) )
      v_minus(s) = ( exp(-rho(s,1)) * r_e * s_gp(s)**2 * d_o_e(s) - sqrt_term ) / ( 2.d0 + s1 * (d_g_e(s) - d_r_e(s)) )
      call compute_velocity_curvature(s, r_h)
      write(unit,"(10es18.9)") s_gp(s), v_plus(s), V_rr_p(s), v_minus(s), V_rr_m(s), l_minus*sqrt(KAPPA)/1.d5, e_minus
    end do
    close(unit)
  end subroutine write_velocity_table

  subroutine compute_velocity_curvature(s, r_h)
    integer, intent(in) :: s
    real(8), intent(in) :: r_h
    real(8) :: o_r, g_r, r_r, o_rr, g_rr, r_rr, term1, term2, term3
    real(8) :: r_h_inv, s_gp_s, r_e_inv_rh_cubed, r_e2_s_gp4_rh4

    r_h_inv = 1.d0 / r_h
    s_gp_s = s_gp(s)
    r_e_inv_rh_cubed = -2.d0 * r_e / (r_e + r_h)**3
    r_e2_s_gp4_rh4 = (r_e * s_gp_s**2 * r_h_inv**2)**2

    o_r  = r_e * s_gp_s**2 * r_h_inv**2 * d_o_e(s)
    g_r  = r_e * s_gp_s**2 * r_h_inv**2 * d_g_e(s)
    r_r  = r_e * s_gp_s**2 * r_h_inv**2 * d_r_e(s)
    o_rr = r_e_inv_rh_cubed * d_o_e(s) + r_e2_s_gp4_rh4 * dd_o_e(s)
    g_rr = r_e_inv_rh_cubed * d_g_e(s) + r_e2_s_gp4_rh4 * dd_g_e(s)
    r_rr = r_e_inv_rh_cubed * d_r_e(s) + r_e2_s_gp4_rh4 * dd_r_e(s)

    term1 = r_h**2 * (g_rr - r_rr + g_r**2 - r_r**2) + 2.d0 * ( exp(-rho_mu_0(s)) * r_h**2 * o_r )**2 + 4.d0 * r_h * r_r - 6.d0
    term2 = r_h**3 * exp(-rho_mu_0(s)) * ( 2.d0 * r_r * o_r - o_rr )
    term3 = r_h**2 * ( g_rr + r_rr + g_r**2 - r_r**2 )

    V_rr_p(s) = (v_plus(s)**2 * term1 + 2.d0 * v_plus(s) * term2 - term3) * r_h_inv**2 / (1.d0 - v_plus(s)**2) * exp(gama_mu_0(s))
    V_rr_m(s) = (v_minus(s)**2 * term1 + 2.d0 * v_minus(s) * term2 - term3) * r_h_inv**2 / (1.d0 - v_minus(s)**2) * exp(gama_mu_0(s))
    l_minus = v_minus(s) * r_h * exp( (gama_mu_0(s) + rho_mu_0(s))/2.d0 ) / sqrt(1.d0 - v_minus(s)**2)
    e_minus = exp( (gama_mu_0(s) + rho_mu_0(s))/2.d0 ) / sqrt(1.d0 - v_minus(s)**2) + ww_mu_0(s) * l_minus
  end subroutine compute_velocity_curvature
end subroutine mass_radius

subroutine solution_properties()
  use para_mod
  use cheb_mod, only: cheb_diff_matrix
  use miscellaneous_mod, only: write_eq_profile
  use toolkit_mod, only: interp, interp_dual, deriv_s, deriv_s_1d, integrate_profiles
  use ad_mod, only: dual, dual_var
  implicit none
  integer :: s, m, ifail
  real(8) :: s_p, r_p
  real(8) :: gama_pole, rho_pole, gama_equator, rho_equator, ww_equator, sphi_equator
  real(8) :: doe, dge, dre, vek
  real(8) :: sqrt_term, mphi_local
  real(8) :: s1, r_h
  real(8), dimension(MDIV) :: scal, acoup, vphi, vel_safe
  real(8), dimension(MDIV,5) :: mu_integrand_buffer
  real(8), dimension(SDIV) :: mass_weight, ang_weight
  real(8), dimension(SDIV) :: d_m, d_m0, d_mp, d_j, d_t
  real(8), dimension(SDIV) :: d_r_e, d_g_e, d_o_e
  real(8), dimension(SDIV) :: dd_r_e, dd_g_e, dd_o_e
  real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  real(8), dimension(SDIV) :: sound_slope, sphi_deriv, pres_deriv
  real(8), dimension(SDIV) :: effective_pressure, scal_potential
  real(8), dimension(SDIV) :: pressure_slope, energy_slope, effective_cs, T_trace
  real(8), dimension(SDIV) :: susceptibility, suscep_slope
  real(8), dimension(SDIV,MDIV) :: rho_0
  real(8), dimension(SDIV,5) :: integrand_buffer
  real(8), dimension(5) :: integral_results, mu_results
  real(8) :: l_minus, e_minus
  character(128) :: profile_file, mphi_str, B_str, Mb_str
  real(8) :: j_local
  logical :: use_scalar ! local snapshot of the flag
  real(8), dimension(SDIV)  :: gamj, schwarz, gg, BV
  real(8) :: delt = 0.005
  type(dual) :: energy_dual, pressure_dual, sphi_dual

  interface
    function p_at_e_dual(ee) result(res)
      import :: dual
      implicit none
      type(dual), intent(in) :: ee
      type(dual) :: res
    end function p_at_e_dual
  end interface

  if (.not. output) return
  call prepare_common_data(rho_0, gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0, use_scalar)

  effective_pressure =  - mphi_r * sphi_mu_0**2 / (4.d0 * pi)
  T_trace = - energy(:,1) + 3.d0 * pressure(:,1)
  scal_potential = ( 2.d0 * pi * B_goal * T_trace &
          * exp( - B_goal * sphi_mu_0**2 / 2.d0 ) + mphi_r ) * exp( gama(:,1) ) &
          + (1.d0-s_gp)**3 / max(s_gp, 1.d-10) / r_e**2 * exp( rho(:,1) ) * d_r_e

  do s = 1, SDIV
    if (energy(s,1) > 0.d0) then
      energy_dual = dual_var(energy(s,1)) ! seeded with dual_var(energy)
      pressure_dual = p_at_e_dual(energy_dual)
      pres_deriv(s) = pressure_dual%der
      ! The code builds energy_dual = dual_var(energy(s,1)), and dual_var fixes its derivative to 1. 
      ! When p_at_e_dual is evaluated on that dual number, the derivative component it returns 
      ! (pressure_dual%der) is literally (dp/de) * energy_dual%der = dp/de. 
      ! Dividing by energy_dual%der would just give the same result, since it is 1. 
      sound_speed(s) = pressure_dual%der !/ energy_dual%der

      pressure_slope(s) = deriv_s_1d(effective_pressure, s)
      energy_slope(s)   = deriv_s_1d(energy(:,1),        s)
      effective_cs(s)   = pressure_slope(s) / energy_slope(s)
      call pressure_derivative_n(energy(s,1), 2, susceptibility(s), ifail)
    else
      sound_speed(s) = 0.d0
      pressure_slope(s) = 0.d0
      energy_slope(s) = 0.d0
    end if
  end do
  do s = 1, SDIV 
    sound_slope (s) = deriv_s_1d(sound_speed, s)
    suscep_slope(s) = deriv_s_1d(susceptibility, s)
  enddo

  if (use_scalar) then
    do s = 1, SDIV
      call interp_dual(s_gp, sphi(:,1), SDIV, dual_var(s_gp(s)), sphi_dual)
      ! seeded with dual_var(s_gp)
      sphi_deriv(s) = sphi_dual%der
    end do
  else
    sphi_deriv = 0.d0
  end if

  call radial_configuration()
  !call to_alexis()
  call mass_radius()
  !call spectral_analysis()
contains

  subroutine spectral_analysis()
    real(8), allocatable, dimension(:) :: x, s_collocate, dgama_coll, src_coll
    real(8), allocatable, dimension(:,:) :: D, op
    real(8), allocatable :: A(:,:), wr(:), wi(:)
    real(8), allocatable :: WORK(:)
    real(8) :: scale
    real(8) :: vl_dummy(1,1), vr_dummy(1,1)
    integer :: N, i, info, unit, ios, LWORK
    
    scale = 2.d0 / max(s_e, 1.d-12)

    write(mphi_str,"(f10.0)") mphi_goal
    profile_file = './Cont/eigenvalues_mphi'//trim(adjustl(mphi_str))//'log'
    open(newunit=unit,file=profile_file,status='unknown',position='append',iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Error opening file eigenvalues.log, IOSTAT=', ios
    end if
    
    do N = 31, 55, 8
      LWORK = 3*(N+1) + 10
      allocate( x(0:N), s_collocate(0:N), dgama_coll(0:N), src_coll(0:N) )
      allocate( D(0:N,0:N), op(0:N,0:N) )
      allocate( A(N+1,N+1), wr(N+1), wi(N+1) )
      allocate( WORK(LWORK) )
      !---------------------------------------------------------------
      ! Chebyshev grid in enthalpy s ∈ [0, s_e]
      !---------------------------------------------------------------
      do i = 0, N
        x(i) = cos( pi * dble(i) / dble(N) )
      end do
      call cheb_diff_matrix(N, D)
      D = scale * D
      s_collocate = 0.5d0 * s_e * (x + 1.d0)
      do i = 0, N
        call interp(s_gp, d_g_e, SDIV, s_collocate(i), dgama_coll(i))
        call interp(s_gp, scal_potential, SDIV, s_collocate(i), src_coll(i))
      enddo
      
      op = matmul(D, D)
      do i = 0, N
        op(i, :) = op(i, :) * ( (1.d0-s_collocate(i))**2 / r_e )**2 &
                  + (dgama_coll(i)*( (1.d0-s_collocate(i))**2 / r_e )**2 &
                  - 2.d0*( (1.d0-s_collocate(i)) / r_e )**2) * D(i, :)
        op(i, i) = op(i, i) - src_coll(i)
      end do

      A = op
      info = 0
      
      call dgeev('N', 'N', N+1, A, N+1, wr, wi, vl_dummy, 1, vr_dummy, 1, WORK, LWORK, info)

      write(unit,"(99es15.6e3)",Advance='NO') maxval(wr), minval(wr)
      deallocate(x, s_collocate, dgama_coll, src_coll)
      deallocate(D, op, A, wr, wi)
      deallocate(WORK)
    enddo
    write(unit,"(10es15.6e3)") h_center, mphi_r, sphi_m
    close(unit)
  end subroutine spectral_analysis

  subroutine radial_configuration()
    write(mphi_str,"(es8.2)") mphi_goal
    write(B_str,"(es8.2)") B_goal
    write(Mb_str,"(f8.0)") Mb_goal*1.d2
    profile_file = "/Users/horay/Data4Projects/crazy/Dat/1dprofile_" &
              // trim(adjustl(B_str)) //"_"// trim(adjustl(mphi_str)) //"_" &
              // trim(adjustl(Mb_str)) // "d"
    call write_eq_profile(profile_file,(SDIV-1)/2,&
         enthalpy(:,1),                           &
         rho_0(:,1)/(KSCALE*C**2)/ n_sat,         &
         sound_speed,                             &
         sound_slope,                             &
         effective_cs,                            &
         suscep_slope,                            &
         sphi(:,1) * sqrt(B_coup),                &
         sphi_deriv * sqrt(B_coup),               &
         effective_pressure,                      & 
         scal_potential  )
  end subroutine radial_configuration

  subroutine to_alexis()
    profile_file = "./Cont/alexis.dat"
    gamj  = sound_speed * (energy(:,1)/pressure(:,1) - 1.d0) / (1.d0+delt)
    schwarz  = pres_deriv / pressure(:,1)/ gamj * delt / (1.d0+delt)
    gg    = -exp(2.d0*rho(:,1)) * pres_deriv / ( energy(:,1) + pressure(:,1) )
    BV    = sqrt( - gg * schwarz ) / r_e * ( 1.d0 - s_gp )**2 * sqrt(KAPPA) / C
    call write_eq_profile(profile_file, (SDIV+1)/2, &
         gama(:,1), rho(:,1), alpha(:,1), & ! 1-3
         enthalpy(:,1), rho_0(:,1)/(KSCALE*C**2), pressure(:,1)/KSCALE, & ! 4-6
         sound_speed, BV) ! 7-8
  end subroutine to_alexis  

end subroutine solution_properties

subroutine prepare_common_data(rho_0, gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0, use_scalar)
  use para_mod
  implicit none
  real(8), dimension(SDIV,MDIV), intent(out) :: rho_0
  real(8), dimension(SDIV), intent(out) :: gama_mu_0, rho_mu_0, ww_mu_0, gama_mu_1, rho_mu_1, sphi_mu_0
  logical, intent(out) :: use_scalar

  interface
    pure elemental function n0_at_e(ee)
      implicit none
      real(8), intent(in) :: ee
      real(8) :: n0_at_e
    end function n0_at_e
  end interface

  use_scalar = (active_theory /= THEORY_GR)

  ! Vectorized extraction of boundary values
  gama_mu_1(:) = gama(:,MDIV)
  rho_mu_1(:)  = rho (:,MDIV)
  gama_mu_0(:) = gama(:,1)
  rho_mu_0(:)  = rho (:,1)
  ww_mu_0 (:)  = ww  (:,1)
  sphi_mu_0(:) = merge( sphi(:,1), 0.d0, use_scalar )

  where (energy > e_surface)
    rho_0 = n0_at_e(energy) * MB * KSCALE * C**2
  elsewhere
    rho_0 = 0.d0
  end where
end subroutine prepare_common_data

