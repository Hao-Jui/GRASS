module rotation_massless
  implicit none

  interface
    subroutine DGESV(N, NRHS, A, LDA, IPIV, B, LDB, INFO)
      integer, intent(in) :: N, NRHS, LDA, LDB
      integer, intent(out) :: INFO
      integer, intent(out) :: IPIV(N)
      real(8), intent(inout) :: A(LDA,N), B(LDB,NRHS)
    end subroutine DGESV
  end interface
contains

subroutine spin_massless
#include "option_macro.h"
  use toolkit_mod
  use para_mod
  use simpson_mod
  use nag_compat_mod, only : d01gaf
  implicit none
  integer :: m, s, n, k, n_of_it, ifail
  real(8) :: r_p, s_p, r_e_old, dif
  real(8) :: r_e_new, r_e_new_sq, grgr, term_in_Omega_h
  real(8) :: gama_pole_h, gama_center_h, gama_equator_h
  real(8) :: rho_pole_h, rho_center_h, rho_equator_h, ww_equator_h
  real(8) :: sphi_pole_h, sphi_center_h, sphi_equator_h
  real(8) :: er2, rho_0
  logical :: valid(MDIV)
  character(32) :: fil1, fil2, fil3, fil4, fil5, fil6

  ! externs used in precompute
  real(8) :: n0_at_e

  ! --- Cached helpers ------------------------------------------------------
  real(8), allocatable :: e_gsm_cache(:,:), e_rsm_cache(:,:), e2alpha_r2_cache(:,:)
  real(8), allocatable :: Acoup_cache(:,:), Acoup4_cache(:,:)
  real(8), allocatable :: dg_s_cache(:,:), dg_m_cache(:,:), d2g_ss_cache(:,:), d2g_mm_cache(:,:)
  real(8), allocatable :: dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
  real(8), allocatable :: ds_s_cache(:,:), ds_m_cache(:,:)

  ! --- Main arrays ---------------------------------------------------------
  real(8) :: sum_rho, sum_gama, sum_omega, sum_sphi
  real(8), dimension(SDIV)      :: gama_mu_1, gama_mu_0, rho_mu_1, rho_mu_0, ww_mu_0
  real(8), dimension(SDIV,MDIV) :: da_dm
  real(8), dimension(SDIV,MDIV) :: S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi
  real(8), dimension(LMAX+1,SDIV) :: D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi
  real(8), dimension(SDIV,LMAX+1) :: D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi
  real(8), allocatable :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
  real(8), dimension(MDIV) :: Int_m
  real(8), dimension(SDIV) :: Int_s

  ! --- Local temps ---------------------------------------------------------
  real(8) :: sgp, mum, s1, s2, m1, gsm, rsm, wwsm, esm, psm, v2sm, scal_p
  real(8) :: d_gama_s, d_gama_m, d_gama_ss, d_gama_mm, d_gama_sm
  real(8) :: d_rho_s, d_rho_m, d_ww_s, d_ww_m, d_sphi_s, d_sphi_m
  real(8) :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8,temp9
  real(8) :: ea, start_local, finish_local

  real(8) :: r_ratio_const

  dif = 1.d0
  n_of_it = 0
  r_e_new = r_e
  r_e_new_sq = r_e_new**2
  r_ratio_const = r_ratio

  if (.not. allocated(e_gsm_cache)) then
    allocate(e_gsm_cache(SDIV,MDIV), e_rsm_cache(SDIV,MDIV), e2alpha_r2_cache(SDIV,MDIV))
    allocate(Acoup_cache(SDIV,MDIV), Acoup4_cache(SDIV,MDIV))
    allocate(dg_s_cache(SDIV,MDIV), dg_m_cache(SDIV,MDIV), d2g_ss_cache(SDIV,MDIV), d2g_mm_cache(SDIV,MDIV))
    allocate(dr_s_cache(SDIV,MDIV), dr_m_cache(SDIV,MDIV), dww_s_cache(SDIV,MDIV), dww_m_cache(SDIV,MDIV))
    allocate(ds_s_cache(SDIV,MDIV), ds_m_cache(SDIV,MDIV))
  end if

  if ( maxval(sphi*sqrt(B_coup)) < 1.d-3 ) sphi = sphi * 10.d0
  if ( any(isnan(sphi)) ) stop "spin_massless: NaN found in sphi"

  call cpu_time(start_local)
  do while( dif > 1.d-7 .or. n_of_it < 2 )

    ! --- Rescale metric potentials ----------------------------------------
    do s = 1, SDIV
      do m = 1, MDIV
        rho  (s,m) = rho  (s,m) / r_e_new_sq
        sphi (s,m) = sphi (s,m) / r_e_new
        gama (s,m) = gama (s,m) / r_e_new_sq
        alpha(s,m) = alpha(s,m) / r_e_new_sq
        ww   (s,m) = ww   (s,m) * r_e_new
      end do
      rho_mu_0 (s) = rho (s,1)
      gama_mu_0(s) = gama(s,1)
      ww_mu_0  (s) = ww  (s,1)
      rho_mu_1 (s) = rho (s,MDIV)
      gama_mu_1(s) = gama(s,MDIV)
    end do
    
    ! --- Equatorial radius update ----------------------------------------
    r_e_old = r_e_new
    r_p     = r_ratio_const * r_e_new
    s_p     = r_p / (r_p + r_e_new)

    call interp(s_gp, sphi(:,MDIV), SDIV, s_p, sphi_pole_h)
    call interp(s_gp, gama_mu_1,    SDIV, s_p, gama_pole_h)
    call interp(s_gp, rho_mu_1,     SDIV, s_p, rho_pole_h)
    call interp(s_gp, gama_mu_0,    SDIV, s_e, gama_equator_h)
    call interp(s_gp, rho_mu_0,     SDIV, s_e, rho_equator_h)
    call interp(s_gp, ww_mu_0,      SDIV, s_e, ww_equator_h)
    call interp(s_gp, sphi(:,1),    SDIV, s_e, sphi_equator_h)

    sphi_center_h = sphi(1,1)
    gama_center_h = gama(1,1)
    rho_center_h  = rho(1,1)

    grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h + &
            B_coup/2.d0 * ( sphi_center_h**2 - sphi_pole_h**2 )

    r_e_new_sq = ( 2.d0 * ( h_center - enthalpy_min ) ) / grgr
    dif = abs(r_e_old - sqrt(r_e_new_sq)) / sqrt(r_e_new_sq)
    r_e_new = sqrt(r_e_new_sq)

    if (r_e_new/r_e_old > 2.d0) stop "spin_massless: r_e changed too much"
    if (r_e_new .ne. r_e_new)  stop "spin_massless: NaN in r_e_new"

    ! --- Angular velocity -------------------------------------------------
    if (r_ratio_const == 1.d0) then
      Omega_c = 0.d0
      ww_equator_h = 0.d0
    else
      grgr = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h + &
             B_coup/2.d0 * ( sphi_equator_h**2 - sphi_pole_h**2 )
      term_in_Omega_h = 1.d0 - exp( r_e_new_sq * grgr )
      if (term_in_Omega_h >= 0.d0) then
        Omega_c = ww_equator_h + exp(r_e_new_sq*rho_equator_h) * sqrt(term_in_Omega_h)
      else
        ! Unphysical state. Set Omega_c to a reasonable value and let the solver try to recover.
        Omega_c = ww_equator_h
      end if
    end if

    ! --- Velocity, enthalpy, EoS -----------------------------------------
    do s = 1, SDIV
      sgp = s_gp(s)
      velocity_sq(s,:) = merge(0.d0, ((Omega_c - ww(s,:)) * (sgp / (1.d0 - sgp)) * &
                            sin_theta(:) * exp(-rho(s,:) * r_e_new_sq))**2, r_ratio_const == 1.d0)
      where (velocity_sq(s,:) > 1.d0) velocity_sq(s,:) = 0.d0

      enthalpy(s,:) = enthalpy_min + 5.d-1 * ( &
            r_e_new_sq * ( gama_pole_h + rho_pole_h - gama(s,:) - rho(s,:) &
            + ( sphi(s,:)**2 - sphi_pole_h**2 )*B_coup/2.d0 ) &
            - log( max(1.d-300, 1.d0-velocity_sq(s,:)) ) )

      valid = (enthalpy(s,:) > enthalpy_min) .and. (sgp <= s_e)

      where (.not. valid)
        enthalpy(s,:) = enthalpy_min
        pressure(s,:) = 0.d0
        energy(s,:)   = 0.d0
      elsewhere
        pressure(s,:) = exp( interp_log_h_to_p( log(enthalpy(s,:)) ) )
        energy(s,:)   = exp( interp_log_p_to_e( log(pressure(s,:)) ) )
      end where

      rho  (s,:) = rho  (s,:) * r_e_new_sq
      gama (s,:) = gama (s,:) * r_e_new_sq
      alpha(s,:) = alpha(s,:) * r_e_new_sq
      sphi (s,:) = sphi (s,:) * r_e_new
    end do

    ! --- Precomputation of derivatives & exponentials --------------------
    do s = 1, SDIV
      do m = 1, MDIV
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

    do s = 1, SDIV
      sgp = s_gp(s)
      s1  = sgp * (1.d0 - sgp)
      do m = 1, MDIV
        mum = mu(m)
        m1  = 1.d0 - mum*mum
        d2g_ss_cache(s,m) = s1*deriv_s(dg_s_cache,s,m) + (1.d0-2.d0*sgp)*dg_s_cache(s,m)
        d2g_mm_cache(s,m) = m1*deriv_m(dg_m_cache,s,m) - 2.d0*mum*dg_m_cache(s,m)
        e_gsm_cache(s,m)      = exp(  0.5d0 * gama(s,m) )
        e_rsm_cache(s,m)      = exp( -       rho(s,m) )
        e2alpha_r2_cache(s,m) = exp(  2.d0 * alpha(s,m) ) * r_e_new_sq
        Acoup_cache(s,m)      = exp( - sphi(s,m)**2 * B_coup / 4.d0 )
        Acoup4_cache(s,m)     = Acoup_cache(s,m)**4
      end do
    end do

    ! --- Build source terms ----------------------------------------------
    S_metric_rho   = 0.d0
    S_metric_sphi  = 0.d0
    S_metric_gama  = 0.d0
    S_metric_omega = 0.d0

    do m = 1, MDIV
      mum = mu(m)
      m1  = 1.d0 - mum*mum
      do s = 1, SDIV
        sgp  = s_gp(s)
        s1   = sgp*(1.d0 - sgp)
        s2   = (sgp/(1.d0 - sgp))**2
        gsm  = gama(s,m)
        rsm  = rho (s,m)
        wwsm = ww  (s,m)
        esm  = energy(s,m)  * Acoup4_cache(s,m)
        psm  = pressure(s,m)* Acoup4_cache(s,m)
        v2sm = velocity_sq(s,m)
        scal_p = sphi(s,m)
        ea = 16.d0 * pi * e2alpha_r2_cache(s,m)

        d_rho_s   = merge(0.d0, dr_s_cache(s,m), s == 1)
        d_rho_m   = merge(0.d0, dr_m_cache(s,m), s == 1)
        d_gama_s  = merge(0.d0, dg_s_cache(s,m), s == 1)
        d_gama_m  = merge(0.d0, dg_m_cache(s,m), s == 1)
        d_ww_s    = merge(0.d0, dww_s_cache(s,m), s == 1)
        d_ww_m    = merge(0.d0, dww_m_cache(s,m), s == 1)
        d_sphi_s  = merge(0.d0, ds_s_cache(s,m), s == 1)
        d_sphi_m  = merge(0.d0, ds_m_cache(s,m), s == 1)
        d_gama_ss = merge(0.d0, d2g_ss_cache(s,m), s == 1)
        d_gama_mm = merge(0.d0, d2g_mm_cache(s,m), s == 1)

        S_metric_rho(s,m) = e_gsm_cache(s,m)*( ea/2.d0 * (esm + psm) * s2 * (1.d0+v2sm)/(1.d0-v2sm) &
          + s2 * m1 * e_rsm_cache(s,m)**2 * ( (s1*d_ww_s)**2 + m1*d_ww_m**2 ) &
          + s1 * d_gama_s - mum * d_gama_m &
          + rsm / 2.d0 * ( ea*psm * s2 - s1*d_gama_s*(s1/2.d0*d_gama_s+1.d0) &
          - d_gama_m*(m1/2.d0*d_gama_m-mum)) )

        S_metric_gama(s,m) = e_gsm_cache(s,m) * ( ea*psm * s2 + gsm/2.d0 * ( ea*psm * s2 &
          - (s1*d_gama_s)**2/2.d0 - m1*d_gama_m**2/2.d0 ) )

        S_metric_omega(s,m) = e_gsm_cache(s,m) * e_rsm_cache(s,m) * &
            ( -ea*(Omega_c-wwsm)*(esm+psm)*s2/(1.d0-v2sm) &
          + wwsm * ( -0.5d0*ea*s2*(((1.d0+v2sm)*esm + 2.d0*v2sm*psm)/(1.d0-v2sm)) &
          - s1 *(2.d0*d_rho_s+0.5d0*d_gama_s) + mum*(2.d0*d_rho_m+0.5d0*d_gama_m) &
          + 0.25d0*s1**2*(4.d0*d_rho_s**2 - d_gama_s**2) + 0.25d0*m1*(4.d0*d_rho_m**2 - d_gama_m**2) &
          - m1*e_rsm_cache(s,m)**2*(sgp**4*d_ww_s**2 + s2*m1*d_ww_m**2) ) )

        S_metric_sphi(s,m) = -s1**2 * d_gama_s * d_sphi_s - m1 * d_gama_m * d_sphi_m &
          + scal_p * ( -2.d0 * pi * B_coup * (esm - 3.d0*psm) ) * s2 * e2alpha_r2_cache(s,m)
      end do
    end do

    ! --- Angular integration ------------------------------------------------
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
    end do

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
        sum_rho   = 0.d0
        sum_sphi  = 0.d0
        sum_gama  = 0.d0
        sum_omega = 0.d0
      end do
    end do

    ! --- Radial integration -------------------------------------------------
    n = 0
    do s = 1, SDIV
      do k = 1, SDIV
        Int_s(k) = f_rho(s,n+1,k) * D1_metric_rho(n+1,k)
      end do
      call d01gaf(s_gp, Int_s, SDIV, sum_rho, er2, ifail)

      do k = 1, SDIV
        Int_s(k) = f_rho(s,n+1,k) * D1_metric_sphi(n+1,k)
      end do
      call d01gaf(s_gp, Int_s, SDIV, sum_sphi, er2, ifail)

      D2_metric_rho  (s,n+1) = sum_rho
      D2_metric_sphi (s,n+1) = sum_sphi
      D2_metric_gama (s,n+1) = 0.d0
      D2_metric_omega(s,n+1) = 0.d0
      sum_rho = 0.d0
      sum_sphi= 0.d0
    end do

    do s = 1, SDIV
      do n = 1, LMAX
        do k = 1, SDIV
          Int_s(k) = f_rho(s,n+1,k) * D1_metric_rho(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_rho, er2, ifail)

        do k = 1, SDIV
          Int_s(k) = f_rho(s,n+1,k) * D1_metric_sphi(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_sphi, er2, ifail)

        do k = 1, SDIV
          Int_s(k) = f_gama(s,n+1,k) * D1_metric_gama(n+1,k)
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_gama, er2, ifail)

        do k = 1, SDIV
          Int_s(k) = merge( f_rho(s,n+1,k)*D1_metric_omega(n+1,k), &
                            f_gama(s,n+1,k)*D1_metric_omega(n+1,k), k < s )
        end do
        call d01gaf(s_gp, Int_s, SDIV, sum_omega, er2, ifail)

        D2_metric_rho  (s,n+1) = sum_rho
        D2_metric_sphi (s,n+1) = sum_sphi
        D2_metric_gama (s,n+1) = sum_gama
        D2_metric_omega(s,n+1) = sum_omega
        sum_rho   = 0.d0
        sum_sphi  = 0.d0
        sum_gama  = 0.d0
        sum_omega = 0.d0
      end do
    end do

    ! --- Summation of coefficients -----------------------------------------
    if (.not. allocated(target_rho)) then
      allocate(target_rho(SDIV,MDIV), target_gama(SDIV,MDIV), target_ww(SDIV,MDIV), target_sphi(SDIV,MDIV))
    end if

    do s = 1, SDIV
      do m = 1, MDIV
        gsm   = gama(s,m)
        rsm   = rho (s,m)
        sum_rho = -exp(-0.5d0*gsm) * P_2n(m,1) * D2_metric_rho(s,1)
        sum_sphi= -              P_2n(m,1) * D2_metric_sphi(s,1)
        sum_gama = 0.d0
        sum_omega= 0.d0
        temp1 = sin_theta(m)

        do n = 1, LMAX
          sum_rho  = sum_rho  - exp(-0.5d0*gsm) * P_2n(m,n+1) * D2_metric_rho(s,n+1)
          sum_sphi = sum_sphi -                  P_2n(m,n+1) * D2_metric_sphi(s,n+1)
          if (m == MDIV) then
            sum_gama  = sum_gama  - 2.d0 / pi * exp(-0.5d0*gsm) * D2_metric_gama(s,n+1)
            sum_omega = sum_omega + exp(rsm-0.5d0*gsm) * D2_metric_omega(s,n+1) / 2.d0
          else
            sum_gama  = sum_gama - 2.d0/pi*exp(-0.5d0*gsm)*(sin_2n_1_theta(m,n)/((2.d0*n-1)*temp1)) * D2_metric_gama(s,n+1)
            sum_omega = sum_omega - exp(rsm-0.5d0*gsm)*(P1_2n_1(m,n+1)/(2.d0*n*(2.d0*n-1)*temp1)) * D2_metric_omega(s,n+1)
          end if
        end do

        target_rho (s,m) = sum_rho
        target_gama(s,m) = sum_gama
        target_ww  (s,m) = sum_omega
        target_sphi(s,m) = sum_sphi
      end do
    end do

    call relaxation(target_rho, target_gama, target_ww, target_sphi, n_of_it)

    ! --- Divergence check & rigid rotation ---------------------------------
    if (abs(rho(2,1))>100.d0 .or. abs(gama(2,1))>300.d0 .or. abs(ww(2,1))>100.d0 &
        .or. abs(sphi(2,1))>10.d0) then
      write(*,"(i5,4es18.9)") n_of_it, rho(2,1), gama(2,1), ww(2,1), sphi(2,1)
      stop "spin_massless: divergence detected"
    end if

    if (r_ratio_const == 1.d0) then
      do s = 1, SDIV
        do m = 1, MDIV
          rho(s,m)  = rho(s,1)
          sphi(s,m) = sphi(s,1)
          gama(s,m) = gama(s,1)
          ww(s,m)   = 0.d0
        end do
      end do
    end if

    ! --- Fourth equation (alpha) -------------------------------------------
    alpha(:,:) = 0.d0
    if (r_ratio_const == 1.d0) then
      da_dm(:,:) = 0.d0
    else
      do s = 2, SDIV
        sgp = s_gp(s)
        s1  = sgp * (1.d0-sgp)
        do m = 1, MDIV
          da_dm(1,m) = 0.0
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

          temp1 = 2.d0 * sgp**2 * (sgp/(1.d0-sgp)) * m1 * d_ww_s * d_ww_m * (1.d0+s1*d_gama_s) &
            - ( (sgp**2 * d_ww_s)**2 - (sgp*d_ww_m/(1.d0-sgp))**2*m1 ) * (-mum + m1*d_gama_m)

          temp2 = 1.d0/( m1 * (1.d0+s1*d_gama_s)**2 + (-mum+m1*d_gama_m)**2 )

          temp3 = s1*d_gama_ss + (s1*d_gama_s)**2
          temp4 = d_gama_m*(-mum + m1*d_gama_m)
          temp5 = ( (s1*(d_rho_s+d_gama_s))**2 - m1*(d_rho_m+d_gama_m)**2 ) &
            * (-mum + m1*d_gama_m)
          temp6 = s1*m1*(  (d_rho_s+d_gama_s)*(d_rho_m+d_gama_m)/2.d0 &
            + d_gama_sm + d_gama_s*d_gama_m  ) * (1.d0 + s1*d_gama_s)
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
      alpha(s,:) = 0.d0
      alpha(s,2:MDIV) = dm * 0.5d0 * cumsum( da_dm(s,1:MDIV-1) + da_dm(s,2:MDIV) )
    end do

    do s = 1, SDIV
      do m = 1, MDIV
        alpha(s,m) = alpha(s,m) - alpha(s,MDIV) + ( gama(s,MDIV)-rho(s,MDIV) )/2.d0
        if(alpha(s,m).ge.300.d0) stop "spin_massless: alpha diverged"
        ww(s,m) = ww(s,m) / r_e_new
      end do
    end do

    n_of_it = n_of_it + 1
    if (n_of_it == 5000 .and. dif > 5.d-6 ) stop "spin_massless: did not converge"
  end do

  call cpu_time(finish_local)
  n_of_relaxation_steps = n_of_relaxation_steps + n_of_it
  !write(*,*) 'spin_massless: Relaxation steps =', n_of_it, ' time =', finish_local-start_local

  ! --- Final bookkeeping ---------------------------------------------------
  r_ratio = r_ratio_const
  omg(:,:) = Omega_c / r_e_new
  Omega_c  = Omega_c / r_e_new
  Omega_e  = Omega_c
  r_e      = r_e_new
  sphi_c   = sphi(1,1) * sqrt(B_coup)
  sphi_m   = maxval( sphi(:,1) * sqrt(B_coup) )
  rho_0    = n0_at_e( energy(1,1) ) * MB
  
  if (output) call output_helper
  
  if (allocated(target_rho)) then
  end if

contains

  subroutine anderson_acceleration(current_field, target_field, history_f, history_g, iter, m_hist)
    real(8), dimension(SDIV,MDIV), intent(inout) :: current_field
    real(8), dimension(SDIV,MDIV), intent(in)    :: target_field
    real(8), dimension(SDIV,MDIV,m_hist), intent(inout) :: history_f, history_g
    integer, intent(in) :: iter, m_hist
    real(8), dimension(m_hist) :: gamma
    real(8), dimension(m_hist, m_hist) :: F_mat
    real(8) :: residual(SDIV,MDIV)
    integer :: k, i, j, info
    integer, dimension(m_hist) :: ipiv
    integer, parameter :: conservative_steps = 5

    residual = current_field - target_field
    history_f(:,:,modulo(iter, m_hist)+1) = residual
    history_g(:,:,modulo(iter, m_hist)+1) = target_field
    k = min(iter, m_hist)

    if (iter < conservative_steps) then
      current_field = 0.9d0 * current_field + 0.1d0 * target_field
    else
      do i = 1, k
        do j = 1, k
          F_mat(i,j) = sum( (history_f(:,:,modulo(iter, m_hist)+1) - history_f(:,:,modulo(iter-i, m_hist)+1)) * &
                            (history_f(:,:,modulo(iter, m_hist)+1) - history_f(:,:,modulo(iter-j, m_hist)+1)) )
        end do
      end do

      gamma(1:k) = 1.d0
      call DGESV(k, 1, F_mat, m_hist, ipiv, gamma, m_hist, info)

      if (info == 0) then
        current_field = (1.d0 - sum(gamma(1:k))) * target_field
        do i=1,k
          current_field = current_field + gamma(i)*history_g(:,:,modulo(iter-i,m_hist)+1)
        end do
      else
        current_field = 0.5d0 * current_field + 0.5d0 * target_field
      end if
    end if
  end subroutine anderson_acceleration

  subroutine relaxation(target_rho, target_gama, target_ww, target_sphi, n_of_it)
    real(8), intent(in) :: target_rho(SDIV,MDIV), target_gama(SDIV,MDIV)
    real(8), intent(in) :: target_ww(SDIV,MDIV), target_sphi(SDIV,MDIV)
    integer, intent(in) :: n_of_it
    integer :: s, m

    integer, parameter :: m_hist = 5
    real(8), allocatable, save :: hist_f_rho(:,:,:), hist_g_rho(:,:,:)
    real(8), allocatable, save :: hist_f_gama(:,:,:), hist_g_gama(:,:,:)
    real(8), allocatable, save :: hist_f_ww(:,:,:), hist_g_ww(:,:,:)
    real(8), allocatable, save :: hist_f_sphi(:,:,:), hist_g_sphi(:,:,:)

    if (.not. allocated(hist_f_rho)) then
      allocate(hist_f_rho(SDIV,MDIV,m_hist), hist_g_rho(SDIV,MDIV,m_hist))
      allocate(hist_f_gama(SDIV,MDIV,m_hist), hist_g_gama(SDIV,MDIV,m_hist))
      allocate(hist_f_ww(SDIV,MDIV,m_hist), hist_g_ww(SDIV,MDIV,m_hist))
      allocate(hist_f_sphi(SDIV,MDIV,m_hist), hist_g_sphi(SDIV,MDIV,m_hist))
    end if

    call anderson_acceleration(rho, target_rho, hist_f_rho, hist_g_rho, n_of_it, m_hist)
    call anderson_acceleration(gama, target_gama, hist_f_gama, hist_g_gama, n_of_it, m_hist)
    call anderson_acceleration(ww, target_ww, hist_f_ww, hist_g_ww, n_of_it, m_hist)
    call anderson_acceleration(sphi, target_sphi, hist_f_sphi, hist_g_sphi, n_of_it, m_hist)

    where(sphi .ne. sphi) sphi = 0.d0
  end subroutine relaxation

  subroutine output_helper()
    call mass_radius
    write(fil1,"(f6.2)") ang_mom
    write(fil2,"(f16.5)") mass_0/MSUN
    write(fil3,"(es15.2)") B_coup
    write(fil4,"(es15.2)") 0.d0
    write(fil5,"(es15.3)") rho_0
    write(fil6,"(f15.3)") sphi_m

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
          velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), sphi(s,m)!*sqrt(B_coup)
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
          velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), sphi(s,m)!*sqrt(B_coup)
      end do
    end do
    close(99)
  end subroutine output_helper

end subroutine spin_massless

end module rotation_massless
