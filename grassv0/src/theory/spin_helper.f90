module spin_helper
  use para_mod
  use toolkit_mod
  use simpson_mod
  use nag_compat_mod, only : d01gaf
  implicit none

  real(8) :: dif
  interface
    subroutine DPOSV(UPLO, N, NRHS, A, LDA, B, LDB, INFO)
      character, intent(in) :: UPLO
      integer, intent(in) :: N, NRHS, LDA, LDB
      integer, intent(out) :: INFO
      real(8), intent(inout) :: A(LDA,N), B(LDB,NRHS)
    end subroutine DPOSV

    subroutine DTRSV(UPLO, TRANS, DIAG, N, A, LDA, X, INCX)
      character, intent(in) :: UPLO, TRANS, DIAG
      integer, intent(in) :: N, LDA, INCX
      real(8), intent(in) :: A(LDA,N)
      real(8), intent(inout) :: X(N)
    end subroutine DTRSV
  end interface

  ! --- Workspace arrays, moved to module scope ---
  real(8), allocatable, target :: e_gsm_cache(:,:), e_rsm_cache(:,:), e2alpha_r2_cache(:,:)
  real(8), allocatable, target :: Acoup_cache(:,:), Acoup4_cache(:,:)
  real(8), allocatable, target :: dg_s_cache(:,:), dg_m_cache(:,:), d2g_ss_cache(:,:), d2g_mm_cache(:,:)
  real(8), allocatable, target :: dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
  real(8), allocatable, target :: ds_s_cache(:,:), ds_m_cache(:,:) 
  real(8), allocatable, target :: mr_cache(:), besseli_cache(:,:), besselk_cache(:,:), wfac_cache(:)
  real(8), allocatable, target :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
  real(8), allocatable, target :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
  real(8), allocatable, target :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
  real(8), allocatable, target :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
  real(8), allocatable, target :: target_field_perturbed(:,:)

contains

  subroutine update_equatorial_radius(r_e_old, r_e_new, dif, sphi_pole_h, gama_pole_h, rho_pole_h, &
                                      gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h, &
                                      sphi_center_h, gama_center_h, rho_center_h)
    real(8), intent(in)    :: r_e_old
    real(8), intent(out)   :: r_e_new, dif
    real(8), intent(out)   :: sphi_pole_h, gama_pole_h, rho_pole_h
    real(8), intent(out)   :: gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h
    real(8), intent(out)   :: sphi_center_h, gama_center_h, rho_center_h
    real(8) :: s_p, r_e_new_sq, grgr
    real(8), dimension(SDIV) :: gama_mu_1, gama_mu_0, rho_mu_1, rho_mu_0, ww_mu_0, sphi_mu_0, sphi_mu_1
    integer :: s

    rho_mu_0 = rho(:,1)
    gama_mu_0 = gama(:,1)
    ww_mu_0  = ww(:,1)
    rho_mu_1 = rho(:,MDIV)
    gama_mu_1 = gama(:,MDIV)
    sphi_mu_0 = sphi(:,1)
    sphi_mu_1 = sphi(:,MDIV)

    s_p     = r_ratio**(1.d0/dble(s_pwr)) / ( 1.d0 + r_ratio**(1.d0/dble(s_pwr)) ) 
    call interp(s_gp, sphi_mu_1, SDIV, s_p, sphi_pole_h   )
    call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole_h   )
    call interp(s_gp, rho_mu_1,  SDIV, s_p, rho_pole_h    )
    call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator_h)
    call interp(s_gp, rho_mu_0,  SDIV, s_e, rho_equator_h )
    call interp(s_gp, ww_mu_0,   SDIV, s_e, ww_equator_h  )
    call interp(s_gp, sphi_mu_0, SDIV, s_e, sphi_equator_h)
    sphi_center_h = sphi(1,1)
    gama_center_h = gama(1,1)
    rho_center_h  = rho(1,1)

    grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h
    if (active_theory /= THEORY_GR) then
      grgr = grgr + B_coup / 2.d0 * ( sphi_center_h**2 - sphi_pole_h**2 )
    end if

    r_e_new_sq = ( 2.d0 * ( h_center - enthalpy_min ) ) / grgr
    r_e_new = sqrt( r_e_new_sq )
    dif = abs(r_e_old - r_e_new) / r_e_new
    
    if (r_e_new / r_e_old > 2) stop 'r_e changed too much'
    if (r_e_new /= r_e_new) stop 'nan in r_e_new'
  end subroutine update_equatorial_radius

  subroutine update_angular_velocity(r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, &
                                     sphi_pole_h, sphi_equator_h, ww_equator_h)
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h
    real(8), intent(in) :: sphi_pole_h, sphi_equator_h, ww_equator_h
    real(8) :: grgr, term_in_Omega_h
    if (r_ratio == 1.d0) then
      Omega_c = 0.d0
    else
      grgr = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h + &
             B_coup/2.d0 * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.d0 - exp( r_e_new**2 * grgr )
      if (term_in_Omega_h >= 0.d0) then
        Omega_c = ww_equator_h + exp(r_e_new**2 * rho_equator_h) * sqrt(term_in_Omega_h)
      else
        write(*,"(10A15)") "gama_pole", "rho_pole", "gama_equator", "rho_equator", "sphi_pole", "sphi_equator"
        write(*,"(10es15.3)") gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, sphi_pole_h, sphi_equator_h
        stop "Omega can't be found; Line 99 of spin helper"
      end if
    end if
  end subroutine update_angular_velocity

  subroutine update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h)
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    integer :: s
    real(8) :: sgp
    logical, dimension(MDIV) :: valid
    do s = 1, SDIV
      sgp = s_gp(s)
      velocity_sq(s,:) = merge(0.d0, ((Omega_c - ww(s,:)) * (sgp / (1.d0 - sgp)) * &
                            sin_theta(:) * exp(-rho(s,:) * r_e_new**2))**2, r_ratio == 1.d0)

      where (velocity_sq(s,:) > 1.d0) velocity_sq(s,:) = 0.d0

      enthalpy(s,:) = enthalpy_min + 5.d-1 * ( &
            r_e_new**2 * ( gama_pole_h + rho_pole_h - gama(s,:) - rho(s,:) &
            + ( sphi(s,:)**2 - sphi_pole_h**2 )*B_coup/2.d0 &
            ) - log( max(1.d-300, 1.d0-velocity_sq(s,:)) )  )

      valid = (enthalpy(s,:) > enthalpy_min) .and. (sgp <= s_e)

      where (.not. valid)
        enthalpy(s,:) = enthalpy_min
        pressure(s,:) = 0.d0
        energy(s,:)   = 0.d0
      elsewhere
        pressure(s,:) = exp( interp_log_h_to_p( log(enthalpy(s,:)) ) )
        energy(s,:)   = exp( interp_log_p_to_e( log(pressure(s,:)) ) )
      end where

      ! Rescale back metric potentials (except omega)
      rho  (s,:) = rho  (s,:) * r_e_new**2
      gama (s,:) = gama (s,:) * r_e_new**2
      alpha(s,:) = alpha(s,:) * r_e_new**2
      sphi (s,:) = sphi (s,:) * r_e_new
    end do
  end subroutine update_eos_and_velocity

  subroutine precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, wfac_cache, besseli_cache, besselk_cache, &
                                                dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                                e2alpha_r2_cache, Acoup_cache, Acoup4_cache)
    real(8), intent(in)  :: r_e_new
    real(8), intent(in)  :: root_mphi_re
    real(8), intent(out) :: mr_cache(:), wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    real(8), intent(out) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(8), intent(out) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(8), intent(out) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(8), intent(out) :: e2alpha_r2_cache(:,:), Acoup_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m, n
    real(8) :: sgp, s1, mum, m1

    do s = 1, SDIV
      mr_cache(s)    = root_mphi_re * s_gp(s) / (1.d0 - s_gp(s))
      wfac_cache(s)  = 1.d0 / (1.d0 - s_gp(s))**2
    end do

    do n = 0, LMAX
      do s = 1, SDIV
        besseli_cache(n+1,s) = besseli(2*n, mr_cache(s))
        besselk_cache(n+1,s) = besselk(2*n, mr_cache(s))
      end do
    end do

    do m = 1, MDIV
      do s = 1, SDIV
        dg_s_cache(s,m) = deriv_s(gama,s,m)
        dg_m_cache(s,m) = deriv_m(gama,s,m)
        dr_s_cache(s,m) = deriv_s(rho ,s,m)
        dr_m_cache(s,m) = deriv_m(rho ,s,m)
        dww_s_cache(s,m)= deriv_s(ww  ,s,m)
        dww_m_cache(s,m)= deriv_m(ww  ,s,m)
        ds_s_cache(s,m) = deriv_s(sphi,s,m)
        ds_m_cache(s,m) = deriv_m(sphi,s,m)
      end do
    end do

    do m = 1, MDIV
      do s = 1, SDIV
        sgp = s_gp(s)
        s1 = sgp*(1.d0 - sgp)
        mum = mu(m)
        m1 = 1.d0 - mum*mum
        d2g_ss_cache(s,m) = s1*deriv_s(dg_s_cache,s,m) + (1.d0-2.d0*sgp)*dg_s_cache(s,m)
        d2g_mm_cache(s,m) = m1*deriv_m(dg_m_cache,s,m) - 2.d0*mum*dg_m_cache(s,m)
        e_gsm_cache(s,m)      = exp(  0.5d0 * gama(s,m) )
        e_rsm_cache(s,m)      = exp( -       rho(s,m) )
        e2alpha_r2_cache(s,m) = exp(  2.d0 * alpha(s,m) ) * r_e_new**2
        Acoup_cache(s,m)      = exp( - sphi(s,m)**2 * B_coup / 4.d0 )
        Acoup4_cache(s,m)     = Acoup_cache(s,m)**4
      end do
    end do
  end subroutine precompute_derivatives_and_bessels

  subroutine build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, &
                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                e2alpha_r2_cache, Acoup4_cache)
    real(8), intent(in)  :: r_e_new
    real(8), intent(out) :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(8), intent(in)  :: dr_s_cache(:,:), dr_m_cache(:,:), dg_s_cache(:,:), dg_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(8), intent(in)  :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(8), intent(in)  :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(8), intent(in)  :: e2alpha_r2_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m
    real(8) :: sgp, mum, s1, s2, m1, gsm, rsm, wwsm, esm, psm, v2sm, scal_p, Vphi
    real(8) :: d_rho_s,d_rho_m,d_gama_s,d_gama_m,d_ww_s,d_ww_m,d_sphi_s,d_sphi_m,d_gama_ss,d_gama_mm
    real(8), dimension(SDIV) :: s1_cache, s2_cache, sgp_cache

    S_metric_rho   = 0.d0
    S_metric_sphi  = 0.d0
    S_metric_gama  = 0.d0
    S_metric_omega = 0.d0

    ! Precompute s-dependent quantities outside m loop
    do s = 1, SDIV
      sgp_cache(s) = s_gp(s)
      s1_cache(s)  = sgp_cache(s) * (1.d0 - sgp_cache(s))
      s2_cache(s) = (sgp_cache(s) / (1.d0 - sgp_cache(s)))**2
    end do

    do m = 1, MDIV
      mum = mu(m)
      m1  = 1.d0 - mum*mum
      do s = 1, SDIV
        sgp  = sgp_cache(s)
        s1   = s1_cache(s)
        s2   = s2_cache(s)
        gsm  = gama(s,m)
        rsm  = rho (s,m)
        wwsm = ww  (s,m)
        esm  = energy(s,m)  * Acoup4_cache(s,m)
        psm  = pressure(s,m)* Acoup4_cache(s,m)
        v2sm = velocity_sq(s,m)
        scal_p = sphi(s,m)*e_gsm_cache(s,m)
        Vphi  = sphi(s,m)**2 * mphi_r * 0.5d0 * e2alpha_r2_cache(s,m)

        d_rho_s   = dr_s_cache(s,m)
        d_rho_m   = dr_m_cache(s,m)
        d_gama_s  = dg_s_cache(s,m)
        d_gama_m  = dg_m_cache(s,m)
        d_ww_s    = dww_s_cache(s,m)
        d_ww_m    = dww_m_cache(s,m)
        d_sphi_s  = ds_s_cache(s,m)
        d_sphi_m  = ds_m_cache(s,m)
        d_gama_ss = d2g_ss_cache(s,m)
        d_gama_mm = d2g_mm_cache(s,m)

        S_metric_rho(s,m) = e_gsm_cache(s,m)*( 16.d0*pi*e2alpha_r2_cache(s,m)/2.d0 * (esm + psm) * s2 * (1.d0+v2sm)/(1.d0-v2sm) &
          + s2 * m1 * e_rsm_cache(s,m)**2 * ( (s1*d_ww_s)**2 + m1*d_ww_m**2 ) &
          + s1 * d_gama_s - mum * d_gama_m &
          + rsm / 2.d0 * ( (16.d0*pi*e2alpha_r2_cache(s,m)*psm-4.d0*Vphi) * s2 - s1*d_gama_s*(s1/2.d0*d_gama_s+1.d0) &
          - d_gama_m*(m1/2.d0*d_gama_m-mum)) )

        S_metric_gama(s,m) = e_gsm_cache(s,m) * ( (16.d0*pi*e2alpha_r2_cache(s,m)*psm-4.d0*Vphi) * s2 &
          + gsm/2.d0 * ( (16.d0*pi*e2alpha_r2_cache(s,m)*psm-4.d0*Vphi) * s2 - (s1*d_gama_s)**2/2.d0 - m1*d_gama_m**2/2.d0 ) )

        S_metric_omega(s,m) = e_gsm_cache(s,m) * e_rsm_cache(s,m) * &
            ( -16.d0*pi*e2alpha_r2_cache(s,m)*(Omega_c-wwsm)*(esm+psm)*s2/(1.d0-v2sm) &
          + wwsm * ( -0.5d0*16.d0*pi*e2alpha_r2_cache(s,m)*s2* (((1.d0+v2sm)*esm + 2.d0*v2sm*psm)/(1.d0-v2sm)) &
          - s1 *(2.d0*d_rho_s+0.5d0*d_gama_s) + mum*(2.d0*d_rho_m+0.5d0*d_gama_m) &
          + 0.25d0*s1**2*(4.d0*d_rho_s**2 - d_gama_s**2) + 0.25d0*m1*(4.d0*d_rho_m**2 - d_gama_m**2) &
          - m1*e_rsm_cache(s,m)**2*(sgp**4*d_ww_s**2 + s2*m1*d_ww_m**2) - 2.d0 * Vphi * s2 ) )

        S_metric_sphi(s,m) = - r_e_new**2 * s2 * scal_p * mphi_r &
          + r_e_new**2 * s2 * scal_p * ( e2alpha_r2_cache(s,m)/r_e_new**2 ) * ( -2.d0 * pi * B_coup * (esm-3.d0*psm) + mphi_r) &
          + scal_p * ( - 2.d0 * s1 * sgp * d_gama_s + s1**2 * d_gama_ss + s1**2 * d_gama_s**2 / 4.d0 &
          + s1 * d_gama_s + m1 * d_gama_mm / 2.d0 + d_gama_m**2 * r_e_new**2 * s2 / 4.d0 - 2.d0 * mum * d_gama_m )
      end do
    end do
  end subroutine build_source_terms

  subroutine angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                 D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    real(8), intent(in)  :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(8), intent(out) :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    integer :: n, k, ifail
    real(8) :: sum_rho, sum_sphi, sum_gama, sum_omega, er2
    real(8), dimension(MDIV) :: Int_m

    n = 0
    do k = 1, SDIV
      Int_m(:) = P_2n(:,n+1) * S_metric_rho(k,:)
      call d01gaf(mu, Int_m, MDIV, sum_rho, er2, ifail)

      Int_m(:) = P_2n(:,n+1) * S_metric_sphi(k,:)
      call d01gaf(mu, Int_m, MDIV, sum_sphi, er2, ifail)

      D1_metric_rho  (n+1,k) = sum_rho
      D1_metric_sphi (n+1,k) = sum_sphi
      D1_metric_gama (n+1,k) = 0.d0
      D1_metric_omega(n+1,k) = 0.d0
    enddo

    do n = 1, LMAX
      do k = 1, SDIV
        Int_m(:) = P_2n(:,n+1) * S_metric_rho(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_rho, er2, ifail)

        Int_m(:) = P_2n(:,n+1) * S_metric_sphi(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_sphi, er2, ifail)

        Int_m(:) = sin_2n_1_theta(:,n) * S_metric_gama(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_gama, er2, ifail)

        Int_m(:) = sin_theta(:) * P1_2n_1(:,n+1) * S_metric_omega(k,:)
        call d01gaf(mu, Int_m, MDIV, sum_omega, er2, ifail)

        D1_metric_rho  (n+1,k) = sum_rho
        D1_metric_sphi (n+1,k) = sum_sphi
        D1_metric_gama (n+1,k) = sum_gama
        D1_metric_omega(n+1,k) = sum_omega
      enddo
    enddo
  end subroutine angular_integration

  subroutine radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
                                D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, &
                                root_mphi_re, wfac_cache, besseli_cache, besselk_cache)
    real(8), intent(in)  :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    real(8), intent(out) :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
    real(8), intent(in)  :: root_mphi_re
    real(8), intent(in)  :: wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    integer :: s, n, k, ifail
    real(8) :: sum_rho, sum_sphi, sum_gama, sum_omega, er2, f_sphi
    real(8), dimension(SDIV) :: Int_s

    n = 0
    do s = 1, SDIV
      do k = 1, SDIV
        Int_s(k) = f_rho(s,n+1,k) * D1_metric_rho(n+1,k)
      end do
      call d01gaf(s_gp, Int_s, SDIV, sum_rho, er2, ifail)

      ! Optimize merge: split into two loops for better branch prediction
      do k = 1, s-1
        f_sphi = besseli_cache(n+1,s) * besselk_cache(n+1,k)
        Int_s(k) = f_sphi * wfac_cache(k) * root_mphi_re * D1_metric_sphi(n+1,k)
      end do
      do k = s, SDIV
        f_sphi = besseli_cache(n+1,k) * besselk_cache(n+1,s)
        Int_s(k) = f_sphi * wfac_cache(k) * root_mphi_re * D1_metric_sphi(n+1,k)
      end do
      call d01gaf(s_gp, Int_s, SDIV, sum_sphi, er2, ifail)

      D2_metric_rho  (s,n+1) = sum_rho
      D2_metric_sphi (s,n+1) = sum_sphi
      D2_metric_gama (s,n+1) = 0.d0
      D2_metric_omega(s,n+1) = 0.d0
    end do

    do s = 1, SDIV
      do n = 1, LMAX
        do k = 1, SDIV
          Int_s(k) = f_rho(s,n+1,k) * D1_metric_rho(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_rho, er2, ifail)

        ! Optimize merge: split into two loops for better branch prediction
        do k = 1, s-1
          f_sphi = besseli_cache(n+1,s) * besselk_cache(n+1,k)
          Int_s(k) = f_sphi * wfac_cache(k) * root_mphi_re * D1_metric_sphi(n+1,k)
        end do
        do k = s, SDIV
          f_sphi = besseli_cache(n+1,k) * besselk_cache(n+1,s)
          Int_s(k) = f_sphi * wfac_cache(k) * root_mphi_re * D1_metric_sphi(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_sphi, er2, ifail)

        do k = 1, SDIV
          Int_s(k) = f_gama(s,n+1,k) * D1_metric_gama(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_gama, er2, ifail)

        ! Optimize merge: split into two loops for better branch prediction
        do k = 1, s-1
          Int_s(k) = f_rho(s,n+1,k) * D1_metric_omega(n+1,k)
        end do
        do k = s, SDIV
          Int_s(k) = f_gama(s,n+1,k) * D1_metric_omega(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_omega, er2, ifail)

        D2_metric_rho  (s,n+1) = sum_rho
        D2_metric_sphi (s,n+1) = sum_sphi
        D2_metric_gama (s,n+1) = sum_gama
        D2_metric_omega(s,n+1) = sum_omega
      end do
    end do
  end subroutine radial_integration

  subroutine sum_coefficients_and_get_targets(target_rho, target_gama, target_ww, target_sphi, &
                                              D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    real(8), intent(out) :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
    real(8), intent(in)  :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
    integer :: s, m, n
    real(8) :: gsm, rsm, temp1, exp_mhalf_gsm, exp_rsm_mhalf_gsm
    real(8) :: sum_rho, sum_sphi, sum_gama, sum_omega

    target_rho  = 0.d0
    target_gama = 0.d0
    target_ww   = 0.d0
    target_sphi = 0.d0
    do s = 1, SDIV
      do m = 1, MDIV
        gsm   = gama(s,m)
        rsm   = rho (s,m)
        temp1 = sin_theta(m)

        exp_mhalf_gsm = exp(-0.5d0*gsm)
        exp_rsm_mhalf_gsm = exp(rsm - 0.5d0*gsm)

        sum_rho = -exp_mhalf_gsm * P_2n(m,1) * D2_metric_rho(s,1)
        sum_sphi= -exp_mhalf_gsm * P_2n(m,1) * D2_metric_sphi(s,1)
        sum_gama = 0.d0
        sum_omega = 0.d0

        do n = 1, LMAX
          sum_rho = sum_rho - exp_mhalf_gsm * P_2n(m,n+1) * D2_metric_rho(s,n+1)
          sum_sphi= sum_sphi- exp_mhalf_gsm * (2.d0*dble(n)+1.d0) * P_2n(m,n+1) * D2_metric_sphi(s,n+1)
          if (m == MDIV) then
            sum_gama  = sum_gama  - 2.d0 / pi * exp_mhalf_gsm * D2_metric_gama(s,n+1)
            sum_omega = sum_omega + exp_rsm_mhalf_gsm*D2_metric_omega(s,n+1)/2.d0
          else
            sum_gama  = sum_gama - 2.d0/pi*exp_mhalf_gsm*(sin_2n_1_theta(m,n)/((2.d0*n-1)*temp1)) * D2_metric_gama(s,n+1)
            sum_omega = sum_omega - exp_rsm_mhalf_gsm*(P1_2n_1(m,n+1)/(2.d0*n*(2.d0*n-1)*temp1)) * D2_metric_omega(s,n+1)
          end if
        end do

        target_rho (s,m) = sum_rho
        target_gama(s,m) = sum_gama
        target_ww  (s,m) = sum_omega
        target_sphi(s,m) = sum_sphi
      end do
    end do
  end subroutine sum_coefficients_and_get_targets

  subroutine update_alpha_potential(r_e_new, dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                    ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_rsm_cache)
    real(8), intent(in) :: r_e_new
    real(8), intent(in) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(8), intent(in) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(8), intent(in) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_rsm_cache(:,:)
    integer :: s, m
    real(8) :: sgp, s1, mum, m1, d_gama_s, d_gama_m, d_rho_s, d_rho_m, d_sphi_s, d_sphi_m
    real(8) :: d_gama_sm, d_ww_s, d_ww_m, d_gama_ss, d_gama_mm
    real(8) :: temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9
    real(8) :: sum_even, sum_odd, sgp_ratio, adj_const
    real(8), dimension(SDIV,MDIV) :: da_dm
    real(8), dimension(SDIV) :: s1_cache, sgp_ratio_cache

    alpha(:,:) = 0.d0
    if (r_ratio == 1.d0) then
      da_dm(:,:) = 0.0
    else
      ! Precompute s-dependent quantities outside m loop
      do s = 2, SDIV
        sgp = s_gp(s)
        s1_cache(s) = sgp * (1.d0-sgp)
        sgp_ratio_cache(s) = sgp / (1.d0-sgp)
      end do
      da_dm(1,:) = 0.0
      do s = 2, SDIV
        do m = 1, MDIV
          sgp = s_gp(s)
          s1  = s1_cache(s)
          mum = mu(m)
          m1  = 1.d0 - mum**2
          d_gama_s  = dg_s_cache(s,m)
          d_gama_m  = dg_m_cache(s,m)
          d_rho_s   = dr_s_cache(s,m)
          d_rho_m   = dr_m_cache(s,m)
          d_sphi_s  = ds_s_cache(s,m)
          d_sphi_m  = ds_m_cache(s,m)
          d_gama_sm = deriv_sm(gama,s,m)
          d_ww_s    = dww_s_cache(s,m)
          d_ww_m    = dww_m_cache(s,m)
          d_gama_ss = d2g_ss_cache(s,m)
          d_gama_mm = d2g_mm_cache(s,m)

          sgp_ratio = sgp_ratio_cache(s)
          temp1 = 2.d0 * sgp**2 * sgp_ratio * m1 * d_ww_s * d_ww_m * (1.d0+s1*d_gama_s) &
            - ( (sgp**2 * d_ww_s)**2 - (sgp*d_ww_m*sgp_ratio)**2*m1 ) * (-mum + m1*d_gama_m)
          temp2 = 1.d0/( m1 * (1.d0+s1*d_gama_s)**2 + (-mum+m1*d_gama_m)**2 )
          temp3 = s1*d_gama_ss + (s1*d_gama_s)**2
          temp4 = d_gama_m*(-mum + m1*d_gama_m)
          temp5 = ( (s1*(d_rho_s+d_gama_s))**2 - m1*(d_rho_m+d_gama_m)**2 ) * (-mum + m1*d_gama_m)
          temp6 = s1*m1*(  (d_rho_s+d_gama_s)*(d_rho_m+d_gama_m)/2.d0 + d_gama_sm + d_gama_s*d_gama_m  ) * (1.d0 + s1*d_gama_s)
          temp7 = s1 * mum * d_gama_s * ( 1.d0 + s1 * d_gama_s )
          temp8 = m1 * (e_rsm_cache(s,m)**2)
          temp9 = -temp2*(-mum+m1*d_gama_m)*( (s1*d_sphi_s)**2 - m1*d_sphi_m**2 ) &
                - m1 * s1 * ( 1.d0 + s1 * d_gama_s ) * 2.d0 * d_sphi_m * d_sphi_s

          da_dm(s,m) = - (d_rho_m+d_gama_m)/2.d0 &
            - temp2*( (temp3 - d_gama_mm - temp4)*(-mum+m1*d_gama_m)/2.d0 &
            + temp5/4.d0 - temp6  + temp7 + temp8*temp1/4.d0 ) + temp9
        end do
      end do
    end if

    do s = 1, SDIV
      alpha(s,1) = 0.d0
      sum_even = 0.d0
      sum_odd  = 0.d0
      do m = 2, MDIV
        if (mod(m,2) == 0) then
          sum_even = sum_even + da_dm(s,m)
          alpha(s,m) = alpha(s,m-1) + 0.5d0*dm*(da_dm(s,m-1) + da_dm(s,m))
        else
          alpha(s,m) = (dm/3.d0) * ( da_dm(s,1) + da_dm(s,m) + 4.d0*sum_even + 2.d0*sum_odd )
          sum_odd = sum_odd + da_dm(s,m)
        end if
      end do
    end do

    do s = 1, SDIV
      adj_const = alpha(s,MDIV) - ( gama(s,MDIV) - rho(s,MDIV) )/2.d0
      alpha(s,:) = alpha(s,:) - adj_const
      if (any(alpha(s,:) .ge. 300.0)) then
          print *, "Error: Alpha fails in row s=", s
          stop "alpha fails"
      end if
      ww(s,:) = ww(s,:) / r_e_new
    end do
  end subroutine update_alpha_potential

  subroutine get_all_targets(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, root_mphi_re, &
                             out_target_rho, out_target_gama, out_target_ww, out_target_sphi)
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re
    real(8), intent(out) :: out_target_rho(SDIV,MDIV), out_target_gama(SDIV,MDIV), out_target_ww(SDIV,MDIV), out_target_sphi(SDIV,MDIV)
   
    call precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, wfac_cache, besseli_cache, besselk_cache, &
         dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
         d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup_cache, Acoup4_cache)

    call build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
         dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
         d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup4_cache)

    call angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
         D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)

    call radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
         D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, root_mphi_re, wfac_cache, &
         besseli_cache, besselk_cache)

    call sum_coefficients_and_get_targets(out_target_rho, out_target_gama, out_target_ww, out_target_sphi, &
         D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    
  end subroutine get_all_targets

  subroutine anderson_acceleration(current_field, target_field, history_f, iter, m_hist)
    real(8), dimension(SDIV,MDIV), intent(inout) :: current_field
    real(8), dimension(SDIV,MDIV), intent(in)    :: target_field
    real(8), dimension(SDIV,MDIV,m_hist), intent(inout) :: history_f
    integer, intent(in) :: iter, m_hist
    real(8), dimension(m_hist) :: gamma, rhs
    real(8), dimension(m_hist, m_hist) :: F_mat
    real(8) :: residual(SDIV,MDIV), df(SDIV,MDIV), accel_field(SDIV,MDIV)
    integer :: k, i, j, info, idx_curr, idx_i, idx_j
    integer, dimension(m_hist) :: idx_array 
    integer, dimension(m_hist) :: ipiv
    real(8), parameter :: blend = 0.3d0
    real(8) :: acc_fac

    residual = current_field - target_field
    ! Optimize modulo: use direct indexing when possible
    idx_curr = modulo(iter, m_hist) + 1
    history_f(:,:,idx_curr) = residual
    k = min(iter, m_hist)
    
    ! Precompute indices to avoid repeated modulo operations
    do i = 1, k
        idx_array(i) = modulo(iter - i, m_hist) + 1
    end do
    
    do i = 1, k
        idx_i = idx_array(i)
        do j = i, k
            idx_j = idx_array(j)
            F_mat(i,j) = sum((history_f(:,:,idx_curr) - history_f(:,:,idx_i)) * &
                             (history_f(:,:,idx_curr) - history_f(:,:,idx_j)))
            if (i /= j) F_mat(j,i) = F_mat(i,j)
        end do
    end do
    do i = 1, k
        idx_i = idx_array(i)
        df = history_f(:,:,idx_curr) - history_f(:,:,idx_i)
        rhs(i) = sum(df * history_f(:,:,idx_curr))
    end do
    gamma(1:k) = rhs(1:k)
    call DPOSV('U', k, 1, F_mat, m_hist, gamma, m_hist, info)

    if (info == 0) then
        accel_field = target_field
        do i = 1, k
            idx_i = idx_array(i)
            accel_field = accel_field - gamma(i) * history_f(:,:,idx_i)
        end do
        acc_fac = 1.d0
        current_field = (1.d0-acc_fac) * current_field + acc_fac * accel_field
    else
        ! Picard damping: cheap but doesn’t use any history
        current_field = blend * current_field + (1.d0-blend) * target_field
    end if
  end subroutine anderson_acceleration

  subroutine relaxation(r_e_new, target_rho, target_gama, target_ww, target_sphi, root_mphi_re, n_of_it)
    real(8), intent(in) :: r_e_new
    real(8), intent(in) :: target_rho(SDIV,MDIV), target_gama(SDIV,MDIV)
    real(8), intent(in) :: target_ww(SDIV,MDIV),  target_sphi(SDIV,MDIV)
    real(8), intent(in) :: root_mphi_re
    integer, intent(in) :: n_of_it
    integer :: s, m
    integer, parameter :: m_hist = 3
    real(8), allocatable, save :: hist_f_rho(:,:,:)
    real(8), allocatable, save :: hist_f_gama(:,:,:)
    real(8), allocatable, save :: hist_f_ww(:,:,:)
    real(8), allocatable, save :: hist_f_sphi(:,:,:)

    select case (trim(relaxation_scheme))
    case ('anderson')
        if (.not. allocated(hist_f_rho)) then
          allocate(hist_f_rho(SDIV,MDIV,m_hist))
          allocate(hist_f_gama(SDIV,MDIV,m_hist))
          allocate(hist_f_ww(SDIV,MDIV,m_hist))
          allocate(hist_f_sphi(SDIV,MDIV,m_hist))
        end if
        call anderson_acceleration(rho, target_rho, hist_f_rho, n_of_it, m_hist)
        call anderson_acceleration(gama, target_gama, hist_f_gama, n_of_it, m_hist)
        call anderson_acceleration(ww, target_ww, hist_f_ww, n_of_it, m_hist)
        call anderson_acceleration(sphi, target_sphi, hist_f_sphi, n_of_it, m_hist)
    case ('newton')
        call newton_krylov('rho',  rho,  n_of_it, r_e_new, root_mphi_re, gama(1,MDIV), rho(1,MDIV), sphi(1,1))
        call newton_krylov('gama', gama, n_of_it, r_e_new, root_mphi_re, gama(1,MDIV), rho(1,MDIV), sphi(1,1))
        call newton_krylov('ww',   ww,   n_of_it, r_e_new, root_mphi_re, gama(1,MDIV), rho(1,MDIV), sphi(1,1))
        call newton_krylov('sphi', sphi, n_of_it, r_e_new, root_mphi_re, gama(1,MDIV), rho(1,MDIV), sphi(1,1))
    case default
        stop "Unknown relaxation scheme specified"
    end select

    where(sphi .ne. sphi) sphi = 0.d0
    do m=1,MDIV
      do s=2,SDIV
        if (sphi(s,m) < 0.d0) then
          sphi(s,m) = sphi(s-1,m) * exp(root_mphi_re * (s_gp(s-1)/(1.d0-s_gp(s-1)) - s_gp(s)/(1.d0-s_gp(s))))
        end if
      end do
    end do
  end subroutine relaxation

  subroutine allocate_workspace
    allocate(e_gsm_cache(SDIV,MDIV), e_rsm_cache(SDIV,MDIV), e2alpha_r2_cache(SDIV,MDIV))
    allocate(Acoup_cache(SDIV,MDIV), Acoup4_cache(SDIV,MDIV))
    allocate(dg_s_cache(SDIV,MDIV), dg_m_cache(SDIV,MDIV), d2g_ss_cache(SDIV,MDIV), d2g_mm_cache(SDIV,MDIV))
    allocate(dr_s_cache(SDIV,MDIV), dr_m_cache(SDIV,MDIV), dww_s_cache(SDIV,MDIV), dww_m_cache(SDIV,MDIV))
    allocate(ds_s_cache(SDIV,MDIV), ds_m_cache(SDIV,MDIV)) 
    allocate(mr_cache(SDIV), besseli_cache(LMAX+1,SDIV), besselk_cache(LMAX+1,SDIV), wfac_cache(SDIV))
    allocate(target_rho(SDIV,MDIV), target_gama(SDIV,MDIV), target_ww(SDIV,MDIV), target_sphi(SDIV,MDIV))
    allocate(S_metric_rho(SDIV,MDIV), S_metric_gama(SDIV,MDIV), S_metric_omega(SDIV,MDIV), S_metric_sphi(SDIV,MDIV))
    allocate(D1_metric_rho(LMAX+1,SDIV), D1_metric_gama(LMAX+1,SDIV), D1_metric_omega(LMAX+1,SDIV), D1_metric_sphi(LMAX+1,SDIV))
    allocate(D2_metric_rho(SDIV,LMAX+1), D2_metric_gama(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1), D2_metric_sphi(SDIV,LMAX+1))
    allocate(target_field_perturbed(SDIV,MDIV))
  end subroutine allocate_workspace

  subroutine deallocate_workspace
    deallocate(e_gsm_cache, e_rsm_cache, e2alpha_r2_cache)
    deallocate(Acoup_cache, Acoup4_cache)
    deallocate(dg_s_cache, dg_m_cache, d2g_ss_cache, d2g_mm_cache)
    deallocate(dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache)
    deallocate(ds_s_cache, ds_m_cache) 
    deallocate(mr_cache, besseli_cache, besselk_cache, wfac_cache)
    deallocate(target_rho, target_gama, target_ww, target_sphi)
    deallocate(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi)
    deallocate(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    deallocate(D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    deallocate(target_field_perturbed)
  end subroutine deallocate_workspace

  subroutine output_helper(D2_metric_rho, D2_metric_omega)
    real(8), intent(in) :: D2_metric_rho(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1)
    real(8) :: r_inf, M2, S3, M4, rho_0
    character(32) :: fil1, fil2, fil3, fil4, fil5, fil6
    integer :: s, m
    external :: mass_radius
    real(8), external :: n0_at_e

    call mass_radius

    r_inf = r_e * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1.d0 - s_gp(SDIV - 1) ))**s_pwr
    M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2.d0 * r_inf**3 * ( C**2 / G / Mass )**3
    S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2.d0 * r_inf**5 * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2.d0 * r_inf**5 * ( C**2 / G / Mass )**5
    
    write(fil1,"(f6.2)") ang_mom
    write(fil2,"(f16.5)") mass_0/MSUN
    write(fil3,"(es15.2)") B_coup
    write(fil4,"(es15.2)") sqrt(mphi_r*1.d10/KAPPA)*l_uni
    rho_0 = n0_at_e( energy(1,1) ) * MB
    write(fil5,"(es15.3)") rho_0
    write(fil6,"(f15.3)") sphi_m

    ! ------------------------------------------------------------------
    !  all the fields are computed and exported as in Einstein frame
    !  except for SPHI due to restart reason.
    ! ------------------------------------------------------------------
    open(98,file="./Cont/"//trim(adjustl(eos_file))//"_J_"//trim(adjustl(fil1))//&
      "_Mb"//trim(adjustl(fil2))//"_B"//trim(adjustl(fil3))//&
      "_mphi"//trim(adjustl(fil4))//"_rhoc"//trim(adjustl(fil5))//"_sphim"//trim(adjustl(fil6))//".dat")
    write(98,"(3i5,99es27.17e3)") SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.d5, &
            energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)) , Omega_c* (C/sqrt(kappa))
    do s = 1, SDIV
      do m = 1, MDIV
        if (enthalpy(s,m) > enthalpy_min) then
          rho_0 = n0_at_e( energy(s,m) ) * MB
        else
          rho_0 = 0.d0
        end if
        write(98,"(99es27.17e3)") s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), ww(s,m) * (C/sqrt(kappa)), &
          pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), enthalpy(s,m), rho_0, &
          velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), sphi(s,m)*sqrt(B_coup)
      end do
    end do
    close(98)

    open(99,file="./Res/res.dat")
    write(99,"(3i5,99es27.17e3)") SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.d5, &
            energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)) , Omega_c* (C/sqrt(kappa))
    do s = 1, SDIV
      do m = 1, MDIV
        if (enthalpy(s,m) > enthalpy_min) then
          rho_0 = n0_at_e( energy(s,m) ) * MB
        else
          rho_0 = 0.d0
        end if
        write(99,"(99es27.17e3)") s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), ww(s,m) * (C/sqrt(kappa)), &
          pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), enthalpy(s,m), rho_0, &
          velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), sphi(s,m)*sqrt(B_coup)
      end do
    end do
    close(99)
  end subroutine output_helper

  ! ---------------------------

  subroutine residual(field_name, x, F_x, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(in) :: x
    real(8), dimension(SDIV,MDIV), intent(out) :: F_x
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re
    real(8), dimension(SDIV,MDIV) :: G_x

    call evaluate_field_target(field_name, x, G_x, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)

    F_x = G_x - x
  end subroutine residual

  subroutine jacobian_vector_product(field_name, x, F_x, v, Jv, eps, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(in) :: x
    real(8), dimension(SDIV,MDIV), intent(in) :: F_x
    real(8), dimension(SDIV*MDIV), intent(in) :: v
    real(8), dimension(SDIV*MDIV), intent(out) :: Jv
    real(8), intent(in) :: eps
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re

    real(8), dimension(SDIV,MDIV) :: F_x_plus_eps_v
    real(8), dimension(SDIV,MDIV) :: x_plus_eps_v
    real(8) :: norm_v, h

    norm_v = sqrt(sum(v**2))
    if (norm_v < 1.d-12) then
        Jv = 0.d0
        return
    end if

    h = eps / norm_v

    ! F(x + h*v)
    x_plus_eps_v = x + reshape(v, [SDIV, MDIV]) * h
    call residual(field_name, x_plus_eps_v, F_x_plus_eps_v, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)

    ! J*v
    Jv = reshape(F_x_plus_eps_v - F_x, [SDIV*MDIV]) / h
  end subroutine jacobian_vector_product

  subroutine newton_krylov(field_name, current_field, iter, r_e_new, root_mphi_re, &
                           gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(inout) :: current_field
    integer, intent(in) :: iter
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h
    real(8), intent(in) :: root_mphi_re

    integer, parameter :: gmres_m = 3
    integer, parameter :: max_iter = 50
    real(8), parameter :: epsilon = 1.d-6
    real(8), parameter :: step_clip = 0.5d0
    real(8), parameter :: armijo_c1 = 1.d-4
    real(8), parameter :: gmres_tol_min = 1.d-6
    real(8), parameter :: gmres_tol_scale = 1.d-2
    integer, parameter :: conservative_steps = 3

    integer :: N, restart_iter, i, j, k
    real(8), allocatable, save :: residual_vec(:), v(:), w(:), delta_x(:), s(:)
    real(8), allocatable, save :: V_basis(:,:), H(:,:), g(:), cs(:), sn(:)
    real(8), dimension(SDIV,MDIV) :: F_x
    real(8) :: norm_residual, beta, r_norm, step_scale
    real(8) :: phi_old, phi_new, alpha
    real(8) :: local_gmres_tol, gmres_target
    real(8), dimension(SDIV,MDIV) :: target_current
    real(8), dimension(SDIV,MDIV) :: trial_field, trial_residual_reshaped

    N = SDIV * MDIV
    if (.not. allocated(residual_vec) .or. size(residual_vec) /= N) then
      if (allocated(residual_vec)) then
        deallocate(residual_vec, v, w, delta_x, s, V_basis, H, g, cs, sn)
      end if
      allocate(residual_vec(N), v(N), w(N), delta_x(N), s(gmres_m))
      allocate(V_basis(N, gmres_m+1))
      allocate(H(gmres_m+1, gmres_m))
      allocate(g(gmres_m+1), cs(gmres_m), sn(gmres_m))
    end if

    if (iter < conservative_steps) then
      call evaluate_field_target(field_name, current_field, target_current, r_e_new, root_mphi_re, &
           gama_pole_h, rho_pole_h, sphi_pole_h)
      current_field = 0.9d0 * current_field + 0.1d0 * target_current
      return
    end if

    ! 1. Compute residual F(x) = G(x) - x
    call residual(field_name, current_field, F_x, r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
    residual_vec = reshape(F_x, [N])
    norm_residual = sqrt(sum(residual_vec**2))
    local_gmres_tol = max(gmres_tol_min, gmres_tol_scale * norm_residual)
    gmres_target = local_gmres_tol * norm_residual
    phi_old = 0.5d0 * norm_residual**2

    if (norm_residual < 1.d-12) return

    ! 2. Solve J * delta_x = -F_x using GMRES
    delta_x = 0.d0
    r_norm = norm_residual
    
    outer_gmres: do restart_iter = 1, max_iter
      k = gmres_m
      V_basis(:,1) = -residual_vec / r_norm
      g = 0.d0
      g(1) = r_norm
      H = 0.d0

      inner_gmres: do j = 1, gmres_m
        v = V_basis(:,j)
        call jacobian_vector_product(field_name, current_field, F_x, v, w, epsilon, &
             r_e_new, root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)

        do i = 1, j
          H(i,j) = dot_product(w, V_basis(:,i))
          w = w - H(i,j) * V_basis(:,i)
        end do

        beta = sqrt(sum(w**2))
        H(j+1,j) = beta
        if (beta < 1.d-14) then
          k = j
          exit inner_gmres
        end if
        V_basis(:,j+1) = w / H(j+1,j)

        do i = 1, j-1
          call apply_givens(cs(i), sn(i), H(i,j), H(i+1,j))
        end do

        call generate_givens(H(j,j), beta, cs(j), sn(j))
        H(j,j) = cs(j)*H(j,j) + sn(j)*beta
        g(j+1) = -sn(j)*g(j)
        g(j)   =  cs(j)*g(j)

        r_norm = abs(g(j+1))
        if (r_norm < gmres_target) then
          k = j
          exit inner_gmres
        end if
      end do inner_gmres

      s(1:k) = g(1:k)
      call DTRSV('U', 'N', 'N', k, H, gmres_m+1, s, 1)
      delta_x = delta_x + matmul(V_basis(:,1:k), s(1:k))

      if (r_norm < gmres_target) exit outer_gmres
    end do outer_gmres

    ! 3. Line search and update (Armijo backtracking on 0.5*||F||^2)
    step_scale = min(1.d0, step_clip / max(1.d-12, sqrt(sum(delta_x**2)/dble(N))))
    delta_x = step_scale * delta_x
    alpha = 1.d0

    do i = 1, 3
      trial_field = current_field + reshape(delta_x, [SDIV,MDIV])
      call residual(field_name, trial_field, trial_residual_reshaped, r_e_new, &
           root_mphi_re, gama_pole_h, rho_pole_h, sphi_pole_h)
     residual_vec = reshape(trial_residual_reshaped, [N])
      phi_new = 0.5d0 * sum(residual_vec**2)
      if (phi_new <= phi_old - armijo_c1 * alpha * norm_residual**2) then
        current_field = trial_field
        return
      else
        delta_x = 0.5d0 * delta_x
        alpha = 0.5d0 * alpha
      end if
    end do

    ! Fallback: simple mixing
    call evaluate_field_target(field_name, current_field, target_current, r_e_new, root_mphi_re, &
         gama_pole_h, rho_pole_h, sphi_pole_h)
    current_field = 0.75d0 * current_field + 0.25d0 * target_current
  end subroutine newton_krylov

  subroutine apply_givens(c, s, v1, v2)
    real(8), intent(in) :: c, s
    real(8), intent(inout) :: v1, v2
    real(8) :: temp
    temp =  c*v1 + s*v2
    v2   = -s*v1 + c*v2
    v1   = temp
  end subroutine apply_givens

  subroutine generate_givens(a, b, c, s)
    real(8), intent(in) :: a, b
    real(8), intent(out) :: c, s
    real(8) :: r
    if (b == 0.d0) then
      c = 1.d0
      s = 0.d0
    else if (abs(b) > abs(a)) then
      r = a / b
      s = 1.d0 / sqrt(1.d0 + r*r)
      c = r * s
    else
      r = b / a
      c = 1.d0 / sqrt(1.d0 + r*r)
      s = r * c
    end if
  end subroutine generate_givens

  subroutine evaluate_field_target(field_name, trial_field, target_out, r_e_new, root_mphi_re, &
                                   gama_pole_h, rho_pole_h, sphi_pole_h)
    character(len=*), intent(in) :: field_name
    real(8), dimension(SDIV,MDIV), intent(in) :: trial_field
    real(8), dimension(SDIV,MDIV), intent(out) :: target_out
    real(8), intent(in) :: r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h 
    real(8), intent(in) :: root_mphi_re
    real(8), allocatable :: backup(:,:), temp_target_rho(:,:), temp_target_gama(:,:), temp_target_ww(:,:), temp_target_sphi(:,:)

    allocate(backup(SDIV,MDIV))

    select case (trim(field_name))
    case ('rho')
      backup = rho
      rho = trial_field
    case ('gama')
      backup = gama
      gama = trial_field
    case ('ww')
      backup = ww
      ww = trial_field
    case ('sphi')
      backup = sphi
      sphi = trial_field
    case default
      stop 'evaluate_field_target: unknown field name'
    end select

    allocate(temp_target_rho(SDIV,MDIV), temp_target_gama(SDIV,MDIV))
    allocate(temp_target_ww(SDIV,MDIV), temp_target_sphi(SDIV,MDIV))
    
    call get_all_targets(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, root_mphi_re, &
         temp_target_rho, temp_target_gama, temp_target_ww, temp_target_sphi)

    select case (trim(field_name))
    case ('rho')
      target_out = temp_target_rho
      rho = backup
    case ('gama')
      target_out = temp_target_gama
      gama = backup
    case ('ww')
      target_out = temp_target_ww
      ww = backup
    case ('sphi')
      target_out = temp_target_sphi
      sphi = backup
    end select

    deallocate(backup, temp_target_rho, temp_target_gama, temp_target_ww, temp_target_sphi)
  end subroutine evaluate_field_target
end module spin_helper
