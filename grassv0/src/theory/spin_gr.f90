module rotation_gr
  implicit none
contains

subroutine spin_gr
  use toolkit_mod
  use simpson_mod
  use nag_compat_mod, only: d01gaf
  use para_mod
  use miscellaneous_mod, only: spectral_tail_fit, composite_richardson
  implicit none
  integer :: m, s, n, k, n_of_it, ifail
  integer :: tail_unit, tail_start, tail_idx
  real(8) :: r_p, s_p, r_e_old, dif
  real(8) :: r_e_new, r_e_new_sq, grgr, term_in_Omega_h
  real(8) :: gama_pole_h, gama_center_h, gama_equator_h
  real(8) :: rho_pole_h, rho_center_h, rho_equator_h, ww_equator_h
  real(8) :: er2, rho_0, r_inf
  real(8) :: c0, c1
  logical :: valid(MDIV), diverg, tail_ok
  character(32) :: fil1, fil2, fil3, fil4, fil5, fil6
  real(8), parameter :: cf = 1.d0

  ! externs used in precompute
  real(8) :: n0_at_e, e_at_p, p_at_h, e_at_h

  ! --- precompute caches ---
  real(8), allocatable :: e_g_half_cache(:,:), e_mhalf_cache(:,:), e_mrho_cache(:,:), e2alpha_r2_cache(:,:)
  real(8), allocatable :: s1_cache(:), s2_cache(:)
  real(8), allocatable :: dr_s_cache(:,:), dr_m_cache(:,:), dg_s_cache(:,:), dg_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
  real(8), allocatable :: dg_ss_cache(:,:), dg_mm_cache(:,:)

  ! --- Main arrays (given by modules) ---
  real(8) :: sum_rho, sum_gama, sum_omega
  real(8), dimension(SDIV) :: gama_mu_1, gama_mu_0, rho_mu_1, rho_mu_0, ww_mu_0
  real(8), dimension(SDIV,MDIV) :: da_dm
  real(8), dimension(SDIV,MDIV) :: S_metric_rho, S_metric_gama, S_metric_omega
  real(8), dimension(LMAX+1,SDIV) :: D1_metric_rho, D1_metric_gama, D1_metric_omega
  real(8), dimension(SDIV,LMAX+1) :: D2_metric_rho, D2_metric_gama, D2_metric_omega
  real(8), dimension(MDIV,SDIV) :: Int_m
  real(8), dimension(SDIV) :: Int_s

  ! local temps
  real(8) :: sgp, mum, s_1, s1, s2, m1, ea
  real(8) :: gsm, rsm, wwsm, esm, psm, v2sm, e_gsm, e_rsm
  real(8) :: d_gama_s, d_gama_m, d_gama_ss, d_gama_mm, d_gama_sm
  real(8) :: d_rho_s, d_rho_m, d_ww_s, d_ww_m
  real(8) :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
  real(8), dimension(MDIV) :: row_rho, row_gama, row_ww, row_energy, row_pressure
  real(8), dimension(MDIV) :: row_v2, row_e_g, row_e_r, row_ea, row_mum, row_m1
  real(8), dimension(MDIV) :: row_dg_s, row_dg_m, row_dr_s, row_dr_m, row_dww_s, row_dww_m
  real(8), dimension(MDIV) :: row_diff, row_term_a, row_term_b, row_term_c, row_term_d, row_term_e, row_term_f
  real(8), dimension(MDIV) :: row_term_g, row_term_h, row_term_i, row_term_j, row_term_k
  real(8), dimension(MDIV) :: row_term_l, row_term_m, row_term_n, row_term_o, row_term_p, row_term_q, row_term_r, row_term_s

  real(8) :: D2_rho_row(LMAX+1)
  real(8) :: D2_gama_row(LMAX+1)
  real(8) :: D2_omega_row(LMAX+1)
  real(8) :: nvec(LMAX)
  real(8), dimension(LMAX) :: P1_term,  sin_term

  ! compute alpha
  real(8) :: alpha_end, shift, inv_r
  ! ---------------------------------------------------------------

  dif = 1.d0
  n_of_it =0
  r_e_new = r_e
  r_e_new_sq = r_e_new**2
  sphi = 0.d0
  s_p = r_ratio**(1.d0/dble(s_pwr)) / ( 1.d0 + r_ratio**(1.d0/dble(s_pwr)) )

  if (.not. allocated(e_g_half_cache)) then
     allocate(e_g_half_cache(SDIV,MDIV), e_mhalf_cache(SDIV,MDIV), e_mrho_cache(SDIV,MDIV), e2alpha_r2_cache(SDIV,MDIV))
     allocate(s1_cache(SDIV), s2_cache(SDIV))
     allocate(dr_s_cache(SDIV,MDIV), dr_m_cache(SDIV,MDIV), dg_s_cache(SDIV,MDIV))
     allocate(dg_m_cache(SDIV,MDIV), dww_s_cache(SDIV,MDIV), dww_m_cache(SDIV,MDIV))
     allocate(dg_ss_cache(SDIV,MDIV), dg_mm_cache(SDIV,MDIV))
  end if
  
  do while( dif > 1.d-7 .or. n_of_it <2 )
      !call cpu_time(start)
      rho    = rho / r_e_new_sq
      gama   = gama / r_e_new_sq
      alpha  = alpha / r_e_new_sq
      ww     = ww * r_e_new
      do s = 1, SDIV
        ! along x-axis
        rho_mu_0 (s) = rho (s,1) ! hat
        gama_mu_0(s) = gama(s,1) ! hat
        ww_mu_0  (s) = ww  (s,1) ! hat
        ! along z-axis
        rho_mu_1 (s) = rho (s,MDIV) ! hat
        gama_mu_1(s) = gama(s,MDIV) ! hat
      enddo
    
      r_e_old = r_e_new ! only to compute dif

      ! --- Compute r_e ---
      call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole_h)
      call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator_h)
      gama_center_h = gama(1,1) ! hat

      call interp(s_gp,  rho_mu_1, SDIV, s_p, rho_pole_h)
      call interp(s_gp,  rho_mu_0, SDIV, s_e, rho_equator_h)
      rho_center_h = rho(1,1) ! hat
    
      grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h ! hat
      r_e_new_sq = 2.d0 * ( h_center - enthalpy_min ) / grgr
      !write(*,"(es15.6,A10,es15.6)") r_e_new, "--->", sqrt(r_e_new_sq)

      r_e_new = sqrt( r_e_new_sq )
      if (r_e_new/r_e_old > 2) stop "r_e changed too much"
      if (r_e_new .ne. r_e_new) stop "nan in r_e_new"

      ! --- Angular velocity (rigid rotation) ---
      if(r_ratio == 1.d0) then
          Omega_c = 0.d0
      else
          call interp( s_gp, ww_mu_0, SDIV, s_e, ww_equator_h )
          grgr = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h ! hat
          term_in_Omega_h = 1.d0 - exp( r_e_new_sq * grgr )
          if (term_in_Omega_h >= 0.d0) then
            Omega_c = ww_equator_h + exp(r_e_new_sq*rho_equator_h) * sqrt(term_in_Omega_h) ! hat
          else
            write(*,*) Omega_c,term_in_Omega_h; stop "Omega can't be found; L88 in spin"
        endif
      endif

      ! > Compute velocity, energy density and pressure
      do s = 1, SDIV
        sgp = s_gp(s)
        velocity_sq(s,:) = merge(0.d0, ((Omega_c - ww(s,:)) * (sgp / (1.d0 - sgp)) * &
                            sin_theta(:) * exp(-rho(s,:) * r_e_new_sq))**2, r_ratio == 1.d0)

        where (velocity_sq(s,:) > 1.d0) velocity_sq(s,:) = 0.d0

        ! enthalpy_min is the assumed small value for the value at the pole
        enthalpy(s,:) = enthalpy_min + 5.d-1 * ( &
                r_e_new_sq * ( gama_pole_h + rho_pole_h - gama(s,:) - rho(s,:) ) - log(1.d0-velocity_sq(s,:))  )

        valid = (enthalpy(s,:) > enthalpy_min) .and. (sgp <= s_e)
        where (.not. valid)
          enthalpy(s,:) = enthalpy_min
          pressure(s,:) = 0.d0
          energy(s,:)   = 0.d0
        elsewhere
          pressure(s,:) = exp( interp_log_h_to_p( log(enthalpy(s,:)) ) )
          energy(s,:)   = exp( interp_log_p_to_e( log(pressure(s,:)) ) )
        end where
      enddo

      rho  =  rho * r_e_new_sq
      gama = gama * r_e_new_sq
      alpha=alpha * r_e_new_sq
      
      if ( disk_present ) call set_disk(r_e_new)

      !------------------------------------------------------------
      ! PRECOMPUTE derivatives, exponentials, and s-factors (on rescaled fields)
      !------------------------------------------------------------
      do s = 1, SDIV
        sgp = s_gp(s)
        s_1 = 1.d0 - sgp
        s1_cache(s) = sgp * s_1 / dble(s_pwr)
        s2_cache(s) = ( sgp / s_1 )**(2*s_pwr)
        do m = 1, MDIV
          e_g_half_cache(s,m) = exp( gama(s,m) / 2.d0 )
          e_mhalf_cache(s,m)  = exp(-gama(s,m) / 2.d0 )
          e_mrho_cache(s,m)   = exp(-rho(s,m))
          e2alpha_r2_cache(s,m)= exp(2.d0*alpha(s,m)) * r_e_new_sq
        end do
      end do

      do s = 1, SDIV
        do m = 1, MDIV
          dg_s_cache(s,m)  = deriv_s(gama, s, m)
          dr_s_cache(s,m)  = deriv_s(rho,  s, m)
          dww_s_cache(s,m) = deriv_s(ww,   s, m)
          dg_m_cache(s,m)  = deriv_m(gama, s, m)
          dr_m_cache(s,m)  = deriv_m(rho,  s, m)
          dww_m_cache(s,m) = deriv_m(ww,   s, m)
          dg_ss_cache(s,m) = deriv_ss(gama, s, m)
          dg_mm_cache(s,m) = deriv_mm(gama, s, m)
        end do
      end do

      !> Compute metric potentials (use caches)
      S_metric_rho   = 0.d0
      S_metric_gama  = 0.d0
      S_metric_omega = 0.d0

      do s = 1, SDIV-1
        sgp = s_gp(s)
        s1  = s1_cache(s)
        s2  = s2_cache(s)
        do m = 1, MDIV
          mum = mu(m)
          m1  = 1.d0 - mum**2
          
          temp1 = 1.d0 / max(1.d-14, 1.d0 - velocity_sq(s,m))
          temp2 = energy(s,m) + pressure(s,m)
          temp3 = (s1 * dww_s_cache(s,m))**2 + m1 * dww_m_cache(s,m)**2
          temp4 = s1 * dg_s_cache(s,m) - mum * dg_m_cache(s,m)
          temp5 = 16.d0 * pi * e2alpha_r2_cache(s,m) * pressure(s,m) * s2
          temp6 = s1 * dg_s_cache(s,m)
          temp7 = m1 * 0.5d0 * dg_m_cache(s,m) - mum
          temp8 = rho(s,m) * 0.5d0 * ( temp5 - temp6 * (0.5d0 * temp6 + 1.d0) - dg_m_cache(s,m) * temp7 )
          
          S_metric_rho(s,m) = e_g_half_cache(s,m) * ( &
               16.d0 * pi * e2alpha_r2_cache(s,m) * 0.5d0 * temp2 * s2 * (1.d0 + velocity_sq(s,m)) * temp1 &
               + s2 * m1 * e_mrho_cache(s,m)**2 * temp3 &
               + temp4 + temp8 )

          temp2 = 0.5d0 * ( (s1 * dg_s_cache(s,m))**2 + m1 * dg_m_cache(s,m)**2 )
          S_metric_gama(s,m) = e_g_half_cache(s,m) * ( temp5 + gama(s,m) * 0.5d0 * ( temp5 - temp2 ) )

          temp2 = energy(s,m) + pressure(s,m)
          temp3 = -16.d0 * pi * e2alpha_r2_cache(s,m) * (Omega_c - ww(s,m)) * temp2 * s2 * temp1
          temp4 = -0.5d0 * 16.d0 * pi * e2alpha_r2_cache(s,m) * s2 * ( &
                  ((1.d0 + velocity_sq(s,m)) * energy(s,m) + 2.d0 * velocity_sq(s,m) * pressure(s,m)) * temp1 )
          temp5 = - s1 * ( 2.d0 * dr_s_cache(s,m) + 0.5d0 * dg_s_cache(s,m) )
          temp6 =   mum * ( 2.d0 * dr_m_cache(s,m) + 0.5d0 * dg_m_cache(s,m) )
          temp7 = 0.25d0 * s1**2 * ( 4.d0 * dr_s_cache(s,m)**2 - dg_s_cache(s,m)**2 )
          temp8 = 0.25d0 * m1 * ( 4.d0 * dr_m_cache(s,m)**2 - dg_m_cache(s,m)**2 )
          temp1 = - m1 * e_mrho_cache(s,m)**2 * s2 * ( (s1 * dww_s_cache(s,m))**2 + m1 * dww_m_cache(s,m)**2 )
          
          S_metric_omega(s,m) = e_g_half_cache(s,m) * e_mrho_cache(s,m) * ( temp3 + ww(s,m) * &
                                ( temp4 + temp5 + temp6 + temp7 + temp8 + temp1 ) )
        enddo
      enddo


      ! ---------------------------------------------------------------
      ! ANGULAR INTEGRATION
      ! ---------------------------------------------------------------
      call integrate_mu(S_metric_rho, S_metric_gama, S_metric_omega, &
                        D1_metric_rho, D1_metric_gama, D1_metric_omega, Int_m)

      ! ---------------------------------------------------------------
      ! RADIAL INTEGRATION
      ! ---------------------------------------------------------------
      ifail = 0
      er2   = 0.d0
      n = 0
      do s = 1, SDIV
        do k = 1, SDIV
          Int_s(k) = f_rho(s,n+1,k) * D1_metric_rho(n+1,k)
        enddo
        call d01gaf(s_gp, Int_s, SDIV, D2_metric_rho(s,n+1), er2, ifail)
      enddo
      D2_metric_gama(:,n+1)  = 0.d0
      D2_metric_omega(:,n+1) = 0.d0

      do n = 1, LMAX
        do s = 1, SDIV
          do k = 1, SDIV
            Int_s(k) = f_rho(s,n+1,k) * D1_metric_rho(n+1,k)
          enddo
          call d01gaf(s_gp, Int_s, SDIV, D2_metric_rho(s,n+1), er2, ifail)

          do k = 1, SDIV
            Int_s(k) = f_gama(s,n+1,k) * D1_metric_gama(n+1,k)
          enddo
          call d01gaf(s_gp, Int_s, SDIV, D2_metric_gama(s,n+1), er2, ifail)

          do k = 1, SDIV
            Int_s(k) = merge( f_rho(s,n+1,k)*D1_metric_omega(n+1,k), &
                f_gama(s,n+1,k)*D1_metric_omega(n+1,k), k < s )
          enddo
          call d01gaf(s_gp, Int_s, SDIV, D2_metric_omega(s,n+1), er2, ifail)
        enddo
      enddo
      
      sum_rho   = 0.d0
      sum_gama  = 0.d0
      sum_omega = 0.d0
      ! ---------------------------------------------------------------
      ! SUMMATION OF COEFFICIENTS & UPDATE
      ! ---------------------------------------------------------------
      do s = 1, SDIV
        ! Pre-load arrays for D2_metric at s
        D2_rho_row   = D2_metric_rho  (s,1:LMAX+1)
        D2_gama_row  = D2_metric_gama (s,1:LMAX+1)
        D2_omega_row = D2_metric_omega(s,1:LMAX+1)

        do m = 1, MDIV
          gsm   = gama(s,m)
          rsm   = rho (s,m)
          wwsm  = ww  (s,m)

          e_gsm = e_mhalf_cache(s,m)
          e_rsm = exp(rsm)
          temp1 = sin_theta(m)
          
          ! > Vectorized sum over n (0:LMAX)
          sum_rho = - e_gsm * sum( P_2n(m,1:LMAX+1) * D2_rho_row )

          if (m == MDIV) then
            ! m = MDIV: simpler formulas without n-dependence except D2 values
            sum_gama  = - (2.d0/pi) * e_gsm * sum( D2_gama_row(2:LMAX+1) )
            sum_omega = e_rsm * e_gsm * sum( D2_omega_row(2:LMAX+1) ) * 0.5d0
          else
            ! Build n-dependent factors for gamma and omega
            nvec = [(n, n=1, LMAX)]

            sin_term = sin_2n_1_theta(m,1:LMAX) / ( (2.d0*nvec - 1.d0) * temp1 )
            P1_term  = P1_2n_1(m,2:LMAX+1) / ( 2.d0*nvec*(2.d0*nvec - 1.d0)*temp1 )

            sum_gama  = - (2.d0/pi) * e_gsm * sum( sin_term * D2_gama_row(2:LMAX+1) )
            sum_omega = - e_rsm * e_gsm * sum( P1_term  * D2_omega_row(2:LMAX+1) )
          endif

          ! > Update fields (vectorizable)
          rho(s,m)  = rsm  + 1.2 * cf * (sum_rho  - rsm)
          gama(s,m) = gsm  + cf * (sum_gama - gsm)
          ww(s,m)   = wwsm + cf * (sum_omega- wwsm)
        enddo
      enddo

      ! check for divergence
      if (abs(rho(2,1))>100.d0 .or. abs(gama(2,1))>300.d0 .or. abs(ww(2,1))>100.d0) then
        write(*,"(3es18.9)") rho(2,1), gama(2,1), ww(2,1)
        write(*,"(3es18.9)") e_at_h (h_center)/(C * C * KSCALE), h_center, r_ratio
        stop "Line 300 in spin"
      endif
#if defined(Fishbone)
#else
      if (r_ratio == 1.d0) then
        rho(s,:)       = rho(s,1)
        gama(s,:)      = gama(s,1)
        alpha(s,:)     = alpha(s,1)
        ww(s,:)        = 0.d0
        pressure(s,:)  = pressure(s,1)
        energy(s,:)    = energy(s,1)
        enthalpy(s,:)  = enthalpy(s,1)
        velocity_sq(s,:)= 0.d0
      endif
#endif   
    ! --- the fourth equation
    
    ! compute first order derivatives of gama
      ! alpha
      alpha(:,:) = 0.d0
      da_dm(:,:) = 0.0

      do s = 2, SDIV
        do m = 1, MDIV
            da_dm(1,m) = 0.d0
            sgp = s_gp(s)
            s1  = s1_cache(s) ! c_1 in the paper
            mum = mu(m)
            m1  = 1.d0 - mum**2
            d_gama_s  = dg_s_cache(s,m)
            d_gama_m  = dg_m_cache(s,m)
            d_rho_s   = dr_s_cache(s,m)
            d_rho_m   = dr_m_cache(s,m)
            d_gama_sm = deriv_sm(gama,s,m)
                        
            d_ww_s    = dww_s_cache(s,m)
            d_ww_m    = dww_m_cache(s,m)
                        
            d_gama_ss = s1 * dg_ss_cache(s,m) &
                      + d_gama_s * ( 1.d0 - 2.d0 * sgp ) / dble(s_pwr)
            d_gama_mm = m1 * dg_mm_cache(s,m) - 2.d0 * mum * d_gama_m


            temp2 = 1.d0/( m1 * (1.d0+s1*d_gama_s)**2 + (-mum+m1*d_gama_m)**2 )

            temp3 = s1 * d_gama_ss + ( s1 * d_gama_s )**2 ! d_gama_ss = (c_1 gama_s),s

            temp4 = d_gama_m * (-mum + m1 * d_gama_m)

            temp5 = ( ( s1 * (d_rho_s+d_gama_s) )**2 - m1*(d_rho_m+d_gama_m)**2 ) &
              * (-mum + m1*d_gama_m)

            temp6 = s1 * m1 * (  (d_rho_s+d_gama_s) * (d_rho_m+d_gama_m) / 2.d0 &
              + d_gama_sm + d_gama_s * d_gama_m ) * ( 1.d0 + s1 * d_gama_s )

            temp7 = s1 * mum * d_gama_s * ( 1.d0 + s1 * d_gama_s )

            temp8 = m1 * exp( -2.d0 * rho(s,m) ) * (sgp / (1.d0-sgp) )**(2*s_pwr)

            temp1 = 2.d0 * s1 * m1 * d_ww_s * d_ww_m * ( 1.d0 + s1 * d_gama_s ) &
              - ( (s1 * d_ww_s)**2 - m1 * d_ww_m**2 ) * (-mum + m1*d_gama_m)
                  
            da_dm(s,m) = - ( d_rho_m + d_gama_m ) * 0.5d0 &
              - temp2 * ( (temp3 - d_gama_mm - temp4) * (-mum + m1 * d_gama_m) * 0.5d0 & ! OK
              + temp5 * 0.25d0 - temp6  + temp7 + temp8 * temp1 * 0.25d0 )
        enddo
      enddo

      do s = 1, SDIV-1
        alpha(s,1) = 0.d0
        do m = 1, MDIV-1
          alpha(s,m+1) = alpha(s,m) + dm * ( da_dm(s,m+1) + da_dm(s,m) ) * 0.5d0
        enddo
      enddo
      alpha(SDIV,:) = 0.d0
      diverg = .false.
      inv_r = 1.d0 / r_e_new

      do s = 1, SDIV
        alpha_end = alpha(s,MDIV)
        shift = ( gama(s,MDIV) - rho(s,MDIV) ) * 0.5d0
        do m = 1, MDIV
          alpha(s,m) = alpha(s,m) - alpha_end + shift
          if ( .not. diverg .and. alpha(s,m) >= 300.d0) then
            diverg = .true.
          end if
          ww(s,m) = ww(s,m) * inv_r
        enddo
      enddo
      if (diverg) stop "L428, alpha fails"
      dif = abs(r_e_old-r_e_new) * inv_r
      n_of_it = n_of_it + 1; n_of_relaxation_steps = n_of_relaxation_steps + 1
      if (n_of_it == 100) stop "Cannot converge; L404 in spin"
      !call cpu_time(finish); write(*,*) finish-start; stop
  enddo
! --- End of iteration

  ! compute omega
  do s = 1, SDIV
    if ( s_gp(s) < s_inner ) then 
      omg(s,:) = Omega_c / r_e_new
    else   
      omg(s,:) = omg(s,:) / r_e_new
    endif
  enddo
  Omega_c  = Omega_c / r_e_new
  Omega_e  = Omega_c
  r_e      = r_e_new
  rho_0    = n0_at_e( energy(1,1) ) * MB

  if (output) call output_helper

contains

  subroutine integrate_mu(src_rho, src_gama, src_omega, dst_rho, dst_gama, dst_omega, scratch)
    real(8), intent(in)    :: src_rho(SDIV,MDIV), src_gama(SDIV,MDIV), src_omega(SDIV,MDIV)
    real(8), intent(out)   :: dst_rho(LMAX+1,SDIV), dst_gama(LMAX+1,SDIV), dst_omega(LMAX+1,SDIV)
    real(8), intent(inout) :: scratch(MDIV,SDIV)
    integer :: n, k, m
    real(8) :: rho_base(MDIV,SDIV), gama_base(MDIV,SDIV), omega_base(MDIV,SDIV)
    real(8) :: weight_vec(MDIV)
    integer :: status_dummy

    rho_base  = transpose(src_rho)
    gama_base = transpose(src_gama)
    omega_base= transpose(src_omega)

    scratch = rho_base
    weight_vec = P_2n(:,1)
    do concurrent (m = 1:MDIV)
      scratch(m,:) = scratch(m,:) * weight_vec(m)
    end do
    do k = 1, SDIV
      call integrate_column_spline(scratch(:,k), mu, dst_rho(1,k), status_dummy)
    end do
    dst_gama(1,:)  = 0.d0
    dst_omega(1,:) = 0.d0

    do n = 1, LMAX
      scratch = rho_base
      weight_vec = P_2n(:,n+1)
      do concurrent (m = 1:MDIV)
        scratch(m,:) = scratch(m,:) * weight_vec(m)
      end do
      do k = 1, SDIV
        call integrate_column_spline(scratch(:,k), mu, dst_rho(n+1,k), status_dummy)
      end do

      scratch = gama_base
      weight_vec = sin_2n_1_theta(:,n)
      do concurrent (m = 1:MDIV)
        scratch(m,:) = scratch(m,:) * weight_vec(m)
      end do
      do k = 1, SDIV
        call integrate_column_spline(scratch(:,k), mu, dst_gama(n+1,k), status_dummy)
      end do

      scratch = omega_base
      weight_vec = sin_theta(:) * P1_2n_1(:,n+1)
      do concurrent (m = 1:MDIV)
        scratch(m,:) = scratch(m,:) * weight_vec(m)
      end do
      do k = 1, SDIV
        call integrate_column_spline(scratch(:,k), mu, dst_omega(n+1,k), status_dummy)
      end do
    end do
  end subroutine integrate_mu

  subroutine output_helper()
    call mass_radius
    r_inf = r_e_new * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1.d0 - s_gp(SDIV - 1) ))**s_pwr

    tail_start = max(2, SDIV - 24)
    open(newunit=tail_unit, file="tail_samples.dat", status='replace', action='write')
    do tail_idx = tail_start, SDIV
      write(tail_unit,'(6es25.16)') s_gp(tail_idx), D2_metric_rho(tail_idx,2), D2_metric_omega(tail_idx,3), &
                                    D2_metric_rho(tail_idx,3), D2_metric_omega(tail_idx,4), D2_metric_rho(tail_idx,4)
    end do
    close(tail_unit)

    call spectral_tail_fit(D2_metric_rho(:,1+1), 2, r_e_new, c0, c1, tail_ok)
    if (tail_ok) then
      M2 = (-0.5d0 * c0) * ( C * C / G / Mass )**3
      write(*,'(A)') "Multipole note: spectral tail fit accepted for M2."
    else
      call composite_richardson(D2_metric_rho(:,1+1), 2, r_e_new, c0, c1, tail_ok)
      if (tail_ok) then
        M2 = (-0.5d0 * c0) * ( C * C / G / Mass )**3
        write(*,'(A)') "Multipole note: Composite Richardson fallback accepted for M2."
      else
        M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2.d0 * r_inf**3 * ( C * C / G / Mass )**3
      end if
    end if

    call spectral_tail_fit(D2_metric_omega(:,2+1), 4, r_e_new, c0, c1, tail_ok)
    if (tail_ok) then
      S3 = (-0.5d0 * c0) * ( C * C / G / Mass )**4 / sqrt(KAPPA)
      write(*,'(A)') "Multipole note: spectral tail fit accepted for S3."
    else
      call composite_richardson(D2_metric_omega(:,2+1), 4, r_e_new, c0, c1, tail_ok)
      if (tail_ok) then
        S3 = (-0.5d0 * c0) * ( C * C / G / Mass )**4 / sqrt(KAPPA)
        write(*,'(A)') "Multipole note: Composite Richardson fallback accepted for S3."
      else
        S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2.d0 * r_inf**5 * ( C * C / G / Mass )**4 / sqrt(KAPPA)
      end if
    end if

    call spectral_tail_fit(D2_metric_rho(:,2+1), 4, r_e_new, c0, c1, tail_ok)
    if (tail_ok) then
      M4 = (-0.5d0 * c0) * ( C * C / G / Mass )**5
      write(*,'(A)') "Multipole note: spectral tail fit accepted for M4."
    else
      call composite_richardson(D2_metric_rho(:,2+1), 4, r_e_new, c0, c1, tail_ok)
      if (tail_ok) then
        M4 = (-0.5d0 * c0) * ( C * C / G / Mass )**5
        write(*,'(A)') "Multipole note: Composite Richardson fallback accepted for M4."
      else
        M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2.d0 * r_inf**5 * ( C * C / G / Mass )**5
      end if
    end if

    !write(fil1,"(f6.2)") ang_mom
    !write(fil2,"(f16.5)") mass_0/MSUN
    !write(fil3,"(es15.3)") rho_0
    !write(fil4,"(i6)") SDIV
    !open(98,file="./Cont/"//trim(adjustl(eos_file))//"_J_"//trim(adjustl(fil1))//&
    !    "_Mb"//trim(adjustl(fil2))//"_rhoc"//trim(adjustl(fil3))//".dat")
    !write(98,"(3i5,99es27.17)") SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.d5, &
    !  energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)) , Omega_c* (C/sqrt(kappa)) 
    !do s = 1, SDIV
    !  do m = 1, MDIV
    !    ! r, \theta, \apha, \gamma, \rho, \omega, \phi, \varepsilon, \rho_0, p
    !    if (enthalpy(s,m) > enthalpy_min) then 
    !      rho_0 = n0_at_e( energy(s,m) ) * MB
    !    else 
    !      rho_0 = 0.d0
    !    endif
    !    write(98,"(99es27.17)") s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), ww(s,m) * (C/sqrt(kappa)), & ! 1-6
    !      pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), enthalpy(s,m), rho_0, & ! 7-10
    !      velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)) ! 11-12
    !  enddo
    !enddo
    !close(98)

    open(99,file="./Res/res.dat")
    write(99,"(3i5,99es27.17)") SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.d5, &
      energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)) , Omega_c* (C/sqrt(kappa))
    if (r_ratio == 1.d0) then
      do s = 1, SDIV
        if (enthalpy(s,1) > enthalpy_min) then
          rho_0 = n0_at_e( energy(s,1) ) * MB
        else
          rho_0 = 0.d0
        end if
        do m = 1, MDIV
          write(99,"(99es27.17)") s_gp(s), mu(m), alpha(s,1), gama(s,1), rho(s,1), 0.d0, &
            pressure(s,1)/KSCALE, energy(s,1)/(C*C*KSCALE), enthalpy(s,1), rho_0, &
            0.d0, omg(s,1) * (C/sqrt(kappa))
        end do
      end do
    else
      do s = 1, SDIV
        do m = 1, MDIV
          if (enthalpy(s,m) > enthalpy_min) then 
            rho_0 = n0_at_e( energy(s,m) ) * MB
          else 
            rho_0 = 0.d0
          endif
          write(99,"(99es27.17)") s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), ww(s,m) * (C/sqrt(kappa)), &
            pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), enthalpy(s,m), rho_0, &
            velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa))
        enddo
      enddo
    end if
    close(99)
  end subroutine output_helper

end subroutine spin_gr

end module rotation_gr
