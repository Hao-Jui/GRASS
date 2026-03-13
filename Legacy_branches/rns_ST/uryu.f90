subroutine uryu
#include "option_macro.h"
  use para_mod
  implicit none
  integer :: m, s
  integer :: n, k, n_of_it
  real(8) :: r_p, s_p
  real(8) :: sum_rho,sum_gama,sum_omega, &  ! intermediate sum
             r_e_old, &                     ! equatorial radius in previus cycle
             dif                    ! difference | r_e_old - r_e |
  real(8) :: d_gama_s,d_gama_m, d_gama_ss,d_gama_mm,d_gama_sm, &
             d_rho_s,d_rho_m,d_ww_s,d_ww_m
  real(8) :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
  real(8) :: m1, s1, s2, ea, mum, sgp, s_1
  real(8) :: rsm, gsm, wwsm, esm, psm, v2sm, e_gsm, e_rsm
  real(8) :: grgr, term_in_Omega_h
  real(8) :: gama_pole_h,gama_max_h,gama_equator_h, &
              rho_pole_h, rho_max_h, rho_equator_h, &
             ww_equator_h, hh_max, omg_max_h, v2_max
  real(8) :: diff_Fmax, Fmax_old
  real(8) :: r_e_new,r_e_new_sq
  real(8) :: deriv_s,deriv_m,deriv_sm,deriv_ss,deriv_mm
  real(8) :: e_at_p,p_at_h,p_at_e,h_at_p,n0_at_e
  real(8) :: p, hh, rho_0, intF
  integer :: ifail, imax
  real(8) :: er2, Int_m(MDIV), Int_s(SDIV)
  real(8), dimension(SDIV) :: gama_mu_1,gama_mu_0,rho_mu_1
  real(8), dimension(SDIV) :: rho_mu_0 ,  ww_mu_0, ee_mu_0, omg_mu_0, v2_mu_0
  real(8), dimension(SDIV,MDIV) :: da_dm,dgds,dgdm
  real(8), dimension(SDIV,MDIV) :: S_metric_rho, S_metric_gama, S_metric_omega
  real(8), dimension(LMAX+1,SDIV) :: D1_metric_rho,D1_metric_gama,D1_metric_omega
  real(8), dimension(SDIV,LMAX+1) :: D2_metric_rho,D2_metric_gama,D2_metric_omega
  real(8) :: ax, fa
  character(8) :: fil1, fil2
  
  external diff_rotation_uryu, rotation_law_uryu

  dif = 1.d0; n_of_it =0
  r_e_new = r_e; r_e_new_sq = r_e_new**2
  
  do while( dif > accuracy .or. n_of_it <2 )
    
    do s = 1, SDIV
      do m = 1, MDIV
        rho  (s,m) = rho  (s,m) / r_e_new_sq ! hat
        gama (s,m) = gama (s,m) / r_e_new_sq ! hat
        alpha(s,m) = alpha(s,m) / r_e_new_sq ! hat
        ww   (s,m) = ww   (s,m) * r_e_new    ! hat
        omg  (s,m) = omg  (s,m) * r_e_new    ! hat
        F_j  (s,m) = F_j  (s,m) / r_e_new
      enddo
      ! along x-axis
      rho_mu_0 (s) = rho        (s,1) ! hat
      gama_mu_0(s) = gama       (s,1) ! hat
      ww_mu_0  (s) = ww         (s,1) ! hat
      ee_mu_0  (s) = energy     (s,1)
      v2_mu_0  (s) = velocity_sq(s,1)
      ! along z-axis
      rho_mu_1 (s) = rho (s,MDIV) ! hat
      gama_mu_1(s) = gama(s,MDIV) ! hat
    enddo
    
  !!! Compute new r_e; consider location at the maximum of energy density
    r_e_old = r_e_new
    r_p     = r_ratio * r_e_new
    s_p     = r_p / (r_p + r_e_new)
    imax    = maxloc( ee_mu_0, 1 )

    call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole_h)
    call interp(s_gp,  rho_mu_1, SDIV, s_p,  rho_pole_h)

    call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator_h)
    call interp(s_gp,  rho_mu_0, SDIV, s_e,  rho_equator_h)
    call interp(s_gp,   ww_mu_0, SDIV, s_e,   ww_equator_h)

    imax    = maxloc( ee_mu_0, 1 ); p = p_at_e(ee_mu_0(imax))
    if (imax > SDIV/2) stop "please check the profile of ee_mu_0; L76 in Uryu"
    rho_max_h = rho_mu_0 ( imax  ) ! hat
    gama_max_h= gama_mu_0( imax  ) ! hat
    hh_max    = h_at_p(p)
    omg_max_h = omg      ( imax,1)
    v2_max    = v2_mu_0  ( imax  )
    
    ! --- v1
      !grgr = 2.d0 * (hh_max-enthalpy_min) + log( 1.d0 - v2_max ) + 2.d0 * intF(omg_max_h, F_j(imax,1) )
      !r_e_new_sq = grgr / (gama_pole_h + rho_pole_h - gama_max_h - rho_max_h)
    ! --- v2
      r_e_new_sq = 2.d0 * (h_center - enthalpy_min) / (gama_pole_h + rho_pole_h - gama_mu_0(1) - rho_mu_0(1))
    !write(*,"(es15.6,A10,es15.6)") r_e_new, "--->", sqrt(r_e_new_sq)
    r_e_new = sqrt(r_e_new_sq)
  
    if ( r_e_new .ne. r_e_new .or. r_e_new/r_e_old > 2.d0 ) then
      write(*,"(10es15.6)") hh_max,h_center, v2_max, gama_max_h, gama_mu_0(1)
      write(*,"(A20,2es15.6)") "[r_e_new, r_e_old]",r_e_new,r_e_old 
      stop 'change in r_e is too dramatic; L94 in spin'
    endif

#if defined(restart)
    Fmax_h = 2.d-2
#else
    if ( omg (1,1)==0.d0 ) then
      grgr = gama_pole_h + rho_pole_h - gama_equator_h - rho_equator_h ! hat
      term_in_Omega_h = 1.d0 - exp( r_e_new_sq * grgr )
      if (term_in_Omega_h >= 0.d0) then
        Omega_e = ( ww_equator_h + exp(r_e_new_sq*rho_equator_h) * sqrt(term_in_Omega_h) ) ! hat
      else
        write(*,"(10es15.6)") rho_mu_0
        stop "L107 in uryu"
      endif
      Fmax_h = 2.d-3
    endif
#endif

    diff_Fmax = 1.d99; Fmax_old = Fmax_h; Fmax_h = Fmax_h/2.d0
    do while( abs(diff_Fmax) > 1.d-7)
        ax = Omega_e
        call zbrent_diff(ax, r_e_new,rho_equator_h,gama_equator_h,ww_equator_h, &
                        rho_pole_h,gama_pole_h,1.d-5, fa, diff_rotation_uryu) ! Fmax_h used here
        Omega_e = fa
        F_equator_h  = (Omega_e-ww_equator_h) / ( exp(2.d0*r_e_new_sq*rho_equator_h) - (Omega_e-ww_equator_h)**2 )
        if ( F_equator_h < 0.d0 ) stop "negative F_equator_h; L120 in uryu"
          
        Omega_c = Omega_e / lambda2
        omg_mu_0(1) = Omega_c
        do s = 2, SDIV*2/3
            rsm = rho(s,1) ! hat
            wwsm= ww (s,1) ! hat
            sgp = s_gp(s)
            ax  = omg_mu_0(s-1)
            call zbrent_rot( ax, r_e_new, rsm, wwsm, sgp, 0.d0, 1.d-5, omg_mu_0(s), rotation_law_uryu)
        enddo
        imax = maxloc( omg_mu_0, 1 )
        omg_max_h = omg_mu_0(imax)
        diff_Fmax = ( Omega_c * lambda1/omg_max_h - 1.d0 )
        Fmax_h = Fmax_h - diff_Fmax*2.d-2
        !write(*,"(3es15.6)") Fmax_h, diff_Fmax
    enddo 
    !write(*,"(A10,es15.6,A10,es15.6)") "Fmax_h", Fmax_old,"--->", Fmax_h
    !write(*,"(10es15.6)") Omega_c, Omega_e, F_equator_h

    !!! compute Omega profile
    Omg     (1,1:MDIV-1) = Omega_c
    Omg(1:3*SDIV/4,MDIV) = Omega_c
    do s = 2, SDIV*3/4
      do m = 1, MDIV-1
        rsm = rho(s,m) ! hat
        wwsm= ww (s,m) ! hat
        mum = mu(m)
        sgp = s_gp(s)
        ax  = Omg(s-1,m)
        call zbrent_rot( ax, r_e_new, rsm, wwsm, sgp, mum, 1.d-5, omg(s,m), rotation_law_uryu)
        F_j(s,m) = (omg(s,m)-wwsm) * sgp**2 * (1.d0-mum**2) &
                / (  (1.d0-sgp)**2 * exp(2.d0*r_e_new_sq*rsm) - (omg(s,m)-wwsm)**2 * sgp**2 * (1.d0-mum**2) )
      enddo
    enddo
    
!!! Compute velocity, energy density and pressure
    do s = 1, SDIV
      sgp = s_gp(s)
      do m = 1, MDIV
        rsm = rho(s,m) ! hat
        if ( r_ratio == 1.d0 ) then
          velocity_sq(s,m) = 0.d0
        else
          velocity_sq(s,m) = ( (omg(s,m)-ww(s,m)) * (sgp/(1.d0-sgp)) &
                * sin_theta(m) * exp(-rsm * r_e_new_sq) )**2
        endif

        if (velocity_sq(s,m) > 1.d0 .or. s > (SDIV/2+2) ) then
          velocity_sq(s,m) = 0.d0
        endif

        ! enthalpy_min is the assumed small value for the value at the pole
        enthalpy(s,m) = enthalpy_min + 5.d-1 * ( &
                r_e_new_sq * ( gama_pole_h + rho_pole_h - gama(s,m) - rsm ) &
                - log(1.d0-velocity_sq(s,m)) ) - intF( omg(s,m), F_j(s,m) )

        if (enthalpy(s,m) <= enthalpy_min .or. sgp > s_e) then
          enthalpy   (s,m) = enthalpy_min
          pressure   (s,m) = 0.d0
          energy     (s,m) = 0.d0
          velocity_sq(s,m) = 0.d0
        else
          pressure(s,m) = p_at_h(enthalpy(s,m))
          energy  (s,m) = e_at_p(pressure(s,m))
        endif
        ! Rescale back metric potentials (except omega)
        rho  (s,m) = rho  (s,m) * r_e_new_sq
        gama (s,m) = gama (s,m) * r_e_new_sq
        alpha(s,m) = alpha(s,m) * r_e_new_sq
      enddo
    enddo

#if defined(debug)
    open(78,file="./chech_uryu.dat")
      do s = 1, SDIV
        do m = 1, MDIV
          write(78,"(99es18.9)") s_gp(s),mu(m),F_j(s,m),enthalpy(s,m),omg(s,m),velocity_sq(s,m),ww(s,m) ! 11-12
        enddo
      enddo
    close(78)
    open(10,file="./Cont/Omega.dat")
    do s = 1, SDIV
      write(10,"(10es15.6)") s_gp(s)/(1.d0-s_gp(s)), omg_mu_0(s), F_j(s,1)
    enddo
    close(10)
#endif

    !!! Compute metric potentials
    S_metric_rho   = 0.d0
    S_metric_gama  = 0.d0
    S_metric_omega = 0.d0

    do s=1,SDIV
      do m=1,MDIV
        rsm   = rho(s,m)
        gsm   = gama(s,m)
        wwsm  = ww(s,m) ! hat
        esm   = energy(s,m)
        psm   = pressure(s,m)
        e_gsm = exp( gsm/2.d0 )
        e_rsm = exp(-rsm)
        v2sm  = velocity_sq(s,m)
        mum   = mu(m)
        m1    = 1.d0 - mum**2
        sgp   = s_gp(s)
        s_1   = 1.d0 - sgp
        s1    = sgp * s_1
        s2    = (sgp / s_1)**2
        ea    = 16.d0 * pi * exp(2.d0*alpha(s,m)) * r_e_new_sq

        if (s == 1) then
          ! at r=0; regularity on z-axis
          d_gama_s = 0.d0
          d_gama_m = 0.d0
          d_rho_s  = 0.d0
          d_rho_m  = 0.d0
          d_ww_s   = 0.d0
          d_ww_m   = 0.d0
        else
          d_rho_s  = deriv_s(rho ,s,m)
          d_rho_m  = deriv_m(rho ,s,m)
          d_gama_s = deriv_s(gama,s,m)
          d_gama_m = deriv_m(gama,s,m)
          d_ww_s   = deriv_s(ww  ,s,m)
          d_ww_m   = deriv_m(ww  ,s,m)
        endif

        S_metric_rho(s,m) = e_gsm*( 0.5*ea*(esm + psm)*s2 * (1.0+v2sm)/(1.0-v2sm) &
          + s2 * m1 * e_rsm**2 * ( (s1*d_ww_s)**2 + m1*d_ww_m**2 ) &
          + s1 * d_gama_s - mum * d_gama_m &
          + 0.5 * rsm * ( ea*psm*s2 - s1*d_gama_s*(0.5*s1*d_gama_s+1.0) &
          - d_gama_m*(0.5*m1*d_gama_m-mum)))

        S_metric_gama(s,m) = e_gsm*(ea*psm*s2 + 0.5*gsm*(ea*psm*s2 &
          - 0.5*(s1*d_gama_s)**2 - 0.5*m1*d_gama_m**2))

        S_metric_omega(s,m) = e_gsm*e_rsm*(-ea*(omg(s,m)-wwsm)*(esm+psm)*s2/(1.0-v2sm) &
          + wwsm*(-0.5*ea*s2*(((1.0+v2sm)*esm + 2.0*v2sm*psm)/(1.0-v2sm)) &
          - s1 *(2*d_rho_s+0.5*d_gama_s) &
          + mum*(2*d_rho_m+0.5*d_gama_m) &
          + 0.25*s1**2*(4*d_rho_s**2 - d_gama_s**2) &
          + 0.25*m1   *(4*d_rho_m**2 - d_gama_m**2) &
          - m1*e_rsm**2*(sgp**4*d_ww_s**2 + s2*m1*d_ww_m**2)))
#if defined(debug)
        if (S_metric_rho(s,m).ne.S_metric_rho(s,m) .or. S_metric_gama(s,m).ne.S_metric_gama(s,m) &
            .or. S_metric_omega(s,m).ne.S_metric_omega(s,m)) then 
            write(*,*) esm,psm,enthalpy(s,m),s_gp(s)/s_e,mu(m)
            write(*,*) n_of_it
            stop "L264 in uryu"
        endif
#endif
      enddo
    enddo


    !--- Angular Integration
    n = 0
    do k = 1, SDIV
      do m = 1, MDIV
        Int_m(m) = P_2n(m,n+1)*S_metric_rho(k,m)
      enddo
      call d01gaf(mu, Int_m, MDIV, sum_rho, er2, ifail)
      D1_metric_rho  (n+1,k) = sum_rho 
      D1_metric_gama (n+1,k) = 0.d0
      D1_metric_omega(n+1,k) = 0.d0
    enddo

    do n = 1, LMAX
      do k = 1, SDIV
        do m = 1, MDIV
          Int_m(m) = P_2n(m,n+1)*S_metric_rho(k,m)
        enddo
        call d01gaf(mu, Int_m, MDIV, sum_rho, er2, ifail)

        do m = 1, MDIV
          Int_m(m) = sin_2n_1_theta(m,n)*S_metric_gama(k,m)
        enddo
        call d01gaf(mu, Int_m, MDIV, sum_gama, er2, ifail)

        do m = 1, MDIV
          Int_m(m) = sin_theta(m)*P1_2n_1(m,n+1)*S_metric_omega(k,m)
        enddo
        call d01gaf(mu, Int_m, MDIV, sum_omega, er2, ifail)

        D1_metric_rho  (n+1,k) = sum_rho  
        D1_metric_gama (n+1,k) = sum_gama
        D1_metric_omega(n+1,k) = sum_omega
      enddo
    enddo

    !--- RADIAL INTEGRATION
    n = 0
    do s=1,SDIV
      do k = 1, SDIV
        Int_s(k) = f_rho(s,n+1,k)*D1_metric_rho(n+1,k)
      enddo
      call d01gaf(s_gp, Int_s, MDIV, sum_rho, er2, ifail)
      D2_metric_rho  (s,n+1) = sum_rho
      D2_metric_gama (s,n+1) = 0.0
      D2_metric_omega(s,n+1) = 0.0
    enddo

    do s=1,SDIV
      do n=1,LMAX
        do k = 1, SDIV
          Int_s(k) = f_rho(s,n+1,  k)*D1_metric_rho(n+1,  k)
        enddo
        call d01gaf(s_gp, Int_s, MDIV, sum_rho, er2, ifail)

        do k = 1, SDIV
          Int_s(k) = f_gama(s,n+1,  k)*D1_metric_gama(n+1,  k)
        enddo
        call d01gaf(s_gp, Int_s, MDIV, sum_gama, er2, ifail)

        do k = 1, SDIV
          Int_s(k) = merge( f_rho(s,n+1,k)*D1_metric_omega(n+1,k), &
              f_gama(s,n+1,k)*D1_metric_omega(n+1,k), k < s )
        enddo
        call d01gaf(s_gp, Int_s, MDIV, sum_omega, er2, ifail)
        D2_metric_rho(s,n+1)   = sum_rho
        D2_metric_gama(s,n+1)  = sum_gama
        D2_metric_omega(s,n+1) = sum_omega
      enddo
    enddo
    
    sum_rho   = 0.d0
    sum_gama  = 0.d0
    sum_omega = 0.d0
    ! SUMMATION of COEFFICIENTS
    do s = 1, SDIV
      do m = 1, MDIV
        gsm   = gama(s,m)
        rsm   = rho (s,m)
        wwsm  = ww  (s,m)
        e_gsm = exp(-0.5*gsm)
        e_rsm = exp(rsm)
        temp1 = sin_theta(m)

        sum_rho = sum_rho - e_gsm * P_2n(m,1) * D2_metric_rho(s,1)

        do n = 1, LMAX
          sum_rho = sum_rho - e_gsm * P_2n(m,n+1) * D2_metric_rho(s,n+1)
          if (m == MDIV) then
            sum_gama  = sum_gama  - 2.d0/pi*e_gsm     *D2_metric_gama(s,n+1)
            sum_omega = sum_omega + e_rsm*e_gsm*D2_metric_omega(s,n+1)/2.d0
          else
            sum_gama  = sum_gama - 2.d0/pi*e_gsm*(sin_2n_1_theta(m,n)/((2.0*n-1)*temp1)) &
              * D2_metric_gama(s,n+1)
            sum_omega = sum_omega - e_rsm*e_gsm*(P1_2n_1(m,n+1)/(2.0*n*(2.0*n-1)*temp1)) &
              * D2_metric_omega(s,n+1)
          endif
        enddo
        rho (s,m) = sum_rho!  rsm + cf*(sum_rho  - rsm)
        gama(s,m) = sum_gama! gsm + cf*(sum_gama - gsm)
        ww  (s,m) = sum_omega !wwsm + cf*(sum_omega- wwsm)
        sum_omega = 0.0
        sum_rho   = 0.0
        sum_gama  = 0.0
      enddo
    enddo
        
    ! check for divergence
    if (abs(rho(2,1))>100.0 .or. abs(gama(2,1))>300.0 .or. abs(ww(2,1))>100.0) then
      write(*,"(3es18.9)") rho(2,1), gama(2,1), ww(2,1)
      stop "Line 360 in spin"
    endif
          
    if (r_ratio == 1.d0) then
      do s=1,SDIV
        do m=1,MDIV
          rho (s,m) = rho(s,1)
          gama(s,m) = gama(s,1)
          ww  (s,m) = 0.d0
          omg (s,m) = 0.d0
        enddo
      enddo
    endif
    ! --- the fourth equation
    
    ! compute first order derivatives of gama

    do s=1,SDIV; do m=1,MDIV
      dgds(s,m) = deriv_s(gama,s,m)
      dgdm(s,m) = deriv_m(gama,s,m)
    enddo; enddo

    ! alpha
    alpha(:,:) = 0.d0
    if (r_ratio == 1.d0) then
      da_dm(:,:) = 0.0
    else
      do s = 2, SDIV
        do m = 1, MDIV
          da_dm(1,m) = 0.0
          sgp = s_gp(s)
          s1  = sgp * (1.d0-sgp)
          mum = mu(m)
          m1  = 1.d0 - mum**2
          d_gama_s  = dgds(s,m)
          d_gama_m  = dgdm(s,m)
          d_rho_s   = deriv_s(rho,s,m)
          d_rho_m   = deriv_m(rho,s,m)
          d_gama_sm = deriv_sm(gama,s,m)
                      
          d_ww_s    = deriv_s(ww,s,m)
          d_ww_m    = deriv_m(ww,s,m)
                      
          d_gama_ss = s1*deriv_s(dgds,s,m) + (1.d0-2.d0*sgp)*d_gama_s
          d_gama_mm = m1*deriv_m(dgdm,s,m) - 2.d0*mum*d_gama_m

          temp1 = 2.d0*sgp**2*(sgp/(1.d0-sgp))*m1*d_ww_s*d_ww_m*(1.d0+s1*d_gama_s) &
            - ( (sgp**2*d_ww_s)**2 - (sgp*d_ww_m/(1.d0-sgp))**2*m1 ) * (-mum + m1*d_gama_m)

          temp2 = 1.d0/( m1 * (1.d0+s1*d_gama_s)**2 + (-mum+m1*d_gama_m)**2 )

          temp3 = s1*d_gama_ss + (s1*d_gama_s)**2

          temp4 = d_gama_m*(-mum + m1*d_gama_m)

          temp5 = ( (s1*(d_rho_s+d_gama_s))**2 - m1*(d_rho_m+d_gama_m)**2 ) &
            * (-mum + m1*d_gama_m)

          temp6 = s1*m1*(  (d_rho_s+d_gama_s)*(d_rho_m+d_gama_m)/2.d0 &
            + d_gama_sm + d_gama_s*d_gama_m  ) * (1.d0 + s1*d_gama_s)

          temp7 = s1 * mum * d_gama_s * ( 1.d0 + s1 * d_gama_s )

          temp8 = m1 * exp( -2.d0 * rho(s,m) )
                
          da_dm(s,m) = - (d_rho_m+d_gama_m)/2.d0 &
            - temp2*( (temp3 - d_gama_mm - temp4)*(-mum+m1*d_gama_m)/2.d0 &
            + temp5/4.d0 - temp6  + temp7 + temp8*temp1/4.d0 )
        enddo
      enddo
    endif

    do s = 1, SDIV
      alpha(s,1) = 0.d0
      do m = 1, 2
        alpha(s,m+1) = alpha(s,m) + dm * ( da_dm(s,m+1) + da_dm(s,m) ) / 2.d0
      enddo
      do m = 4, MDIV
        call d01gaf( mu(1:m), da_dm(s,1:m), m, alpha(s,m), er2, ifail)
        !alpha(s,m+1) = alpha(s,m) + dm * ( da_dm(s,m+1) + da_dm(s,m) ) / 2.d0
      enddo
    enddo

    do s = 1, SDIV
      do m = 1, MDIV
        alpha(s,m) = alpha(s,m) - alpha(s,MDIV) + ( gama(s,MDIV)-rho(s,MDIV) )/2.d0
        if(alpha(s,m).ge.300.0) then
          write(*,*) n_of_it
          stop "L440, alpha fails"
        endif
        ww(s,m) =  ww(s,m) / r_e_new
        omg(s,m)= omg(s,m) / r_e_new
        F_j(s,m)= F_j(s,m) * r_e_new
      enddo
    enddo
    do s = 2, SDIV
      sgp = s_gp(s)
      do m = 2, MDIV
        mum = mu(m)
        F_j(s,m) = (omg(s,m)-ww(s,m))*sgp**2*(1.d0-mum**2)*exp(-2.d0*r_e_new_sq*rsm) &
                / (  (1.d0-sgp**2) - (omg(s,m)-ww(s,m))**2 * sgp**2 * (1.d0-mum**2) * exp(-2.d0*r_e_new_sq*rsm) )
      enddo
    enddo

    dif = abs(r_e_old-r_e_new)/r_e_new
    n_of_it = n_of_it + 1
    !write(*,"(i5,es15.6)") n_of_it, dif
    if (n_of_it == 100) stop "Cannot converge; L404 in uryu"
  enddo
  ! --- End of iteration 
  
  Omega_c = Omega_c / r_e_new
  Omega_e = Omega_e / r_e_new
  r_e = r_e_new

  if (output) then
    write(fil1,"(f4.2)") ang_mom
    write(fil2,"(f4.2)") mass_0/MSUN
    open(98,file="./Cont/J"//trim(adjustl(fil1))//"_Mb"//trim(adjustl(fil2))// &
                  "_Uryu.dat")
#if defined(matlab)
#else
        write(98,"(2i5,99es27.17)") SDIV, MDIV, r_e*sqrt(KAPPA)/1.d5, &
        energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)) , Omega_c* (C/sqrt(kappa)) 
#endif
    do s = 1, SDIV
      do m = 1, MDIV
        ! r, \theta, \apha, \gamma, \rho, \omega, \phi, \varepsilon, \rho_0, p
        if (enthalpy(s,m) > enthalpy_min) then 
          rho_0 = n0_at_e( energy(s,m) ) * MB
        else 
          rho_0 = 0.d0
          !omg(s,m) = 0.d0
        endif
        write(98,"(99es27.17)") s_gp(s),mu(m),alpha(s,m),gama(s,m),rho(s,m),ww(s,m)* (C/sqrt(kappa)), & ! 1-6
        pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), enthalpy(s,m), rho_0, & ! 7-10
        velocity_sq(s,m), omg(s,m)* (C/sqrt(kappa)) ! 11-12
      enddo
    enddo
    close(98)
  endif

end subroutine uryu

subroutine diff_rotation_uryu(x, fx, re, rho_e, g_e, w_e, rho_p, g_p)
! Compute Omega_e
  use para_mod, only: lambda1, lambda2, uyru_p, uyru_q, Fmax_h
  implicit none
  real(8), intent(in) :: x, re, rho_e, g_e, w_e, rho_p, g_p
  real(8), intent(out):: fx
  real(8) :: F_e, ocre, RHS, AA_h, BB_h, aa, bb

  ocre = x / lambda2
  F_e  = (x - w_e) / ( exp(2.d0*re**2*rho_e) - (x-w_e)**2 )

  aa = AA_h(F_e, Fmax_h)
  bb = BB_h(F_e, Fmax_h)

  RHS = F_e * x - aa * ocre / 4.d0 * &
        (2.d0 * aa / bb * atan(F_e**2/aa**2) &
        - sqrt(2.d0) * ( atan(1.d0 - F_e*sqrt(2.d0)/aa) - atan(1.d0 + F_e*sqrt(2.d0)/aa) ) &
        + sqrt(2.d0) * ATANH( aa * F_e * sqrt(2.d0) / (F_e**2 + aa**2) ) )

  fx = re**2 * (g_e + rho_e - g_p - rho_p) + log( 1.d0 - ((x-w_e)*exp(-re**2*rho_e))**2 ) &
      + 2.d0 * RHS
  !write(*,"(A10,19es15.6)") "brent:", x, fx, ((x-w_e)*exp(-re**2*rho_e))**2
  !if (fx.ne.fx) stop "L536"

end subroutine diff_rotation_uryu

subroutine rotation_law_uryu(x, fx, re, rho_p, ww_p, sgp, mugp)
  use para_mod, only: uyru_p,uyru_q,Omega_c,F_equator_h, Fmax_h
  implicit none
  real(8), intent(in) :: x, re, rho_p, ww_p, sgp, mugp
  real(8), intent(out):: fx
  real(8) :: Fj, AA_h, BB_h, aa, bb

  aa = AA_h(F_equator_h, Fmax_h)
  bb = BB_h(F_equator_h, Fmax_h)

  Fj = (x-ww_p) * sgp**2 * (1.d0-mugp**2) / ( (1.d0-sgp)**2 * exp(2.d0*re**2*rho_p) - (x-ww_p)**2 * sgp**2 * (1.d0-mugp**2) )
  fx = x / Omega_c * ( 1.d0 + (Fj / aa)**(uyru_p+uyru_q) ) - ( 1.d0 + (Fj / bb)**uyru_p )
  
  !write(*,"(10es15.6)") x,fx,tmp2,BB,AA

end subroutine rotation_law_uryu

real(8) function intF(x, F_at_x)
  use para_mod, only: F_equator_h, Fmax_h, Omega_c
  implicit none
  real(8), intent(in) :: x, F_at_x
  real(8) :: AA_h, BB_h, aa, bb

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

real(8) function AA_h(F_e, F_m)
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

real(8) function BB_h(F_e, F_m)
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