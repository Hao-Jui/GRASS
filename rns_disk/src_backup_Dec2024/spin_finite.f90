subroutine spin
#include "option_macro.h"
  use para_mod
  implicit none
  integer :: m, s
  integer :: n, k, n_of_it
  real(8) :: r_p, s_p
  real(8) :: sum_rho,sum_gama,sum_omega, &  ! intermediate sum
             r_e_old, &                     ! equatorial radius in previus cycle
             dif                            ! difference | r_e_old - r_e |
  real(8) :: d_gama_s,d_gama_m, d_gama_ss,d_gama_mm,d_gama_sm, &
             d_rho_s,d_rho_m,d_ww_s,d_ww_m
  real(8) :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
  real(8) :: m1, s1, s2, ea, mum, sgp, s_1
  real(8) :: rsm, gsm, wwsm, esm, psm, v2sm, e_gsm, e_rsm
  real(8) :: grgr, term_in_Omega_h
  real(8) :: gama_pole_h,gama_center_h,gama_equator_h, &
             rho_pole_h,rho_center_h,rho_equator_h, ww_equator_h
  real(8) :: r_e_new,r_e_new_sq
  real(8) :: deriv_s,deriv_m,deriv_sm!,deriv_ss,deriv_mm
  real(8) :: e_at_p,p_at_h,n0_at_e
  real(8) :: rho_0
  integer :: ifail
  real(8) :: er2, Int_m(MDIV), Int_s(SDIV)
  real(8), dimension(SDIV) :: gama_mu_1,gama_mu_0,rho_mu_1,rho_mu_0,ww_mu_0
  real(8), dimension(SDIV,MDIV) :: da_dm,dgds,dgdm
  real(8), dimension(SDIV,MDIV) :: S_metric_rho, S_metric_gama, S_metric_omega
  real(8), dimension(LMAX+1,SDIV) :: D1_metric_rho, D1_metric_gama, D1_metric_omega
  real(8), dimension(SDIV,LMAX+1) :: D2_metric_rho, D2_metric_gama, D2_metric_omega
  character(8) :: fil1, fil2


  dif = 1.d0; n_of_it =0

  r_e_new = r_e; r_e_new_sq = r_e_new**2
  
  do while( dif > accuracy .or. n_of_it <2 )
      do s = 1, SDIV
        do m = 1, MDIV
          rho  (s,m) = rho  (s,m) / r_e_new_sq ! hat
          gama (s,m) = gama (s,m) / r_e_new_sq ! hat
          alpha(s,m) = alpha(s,m) / r_e_new_sq ! hat
          ww   (s,m) = ww   (s,m) * r_e_new    ! hat
        enddo
        ! along x-axis
        rho_mu_0 (s) = rho (s,1) ! hat
        gama_mu_0(s) = gama(s,1) ! hat
        ww_mu_0  (s) = ww  (s,1) ! hat
        ! along z-axis
        rho_mu_1 (s) = rho (s,MDIV) ! hat
        gama_mu_1(s) = gama(s,MDIV) ! hat
      enddo
    
      ! Compute new r_e
      r_e_old = r_e_new ! only to compute dif
      r_p     = r_ratio * r_e_new
      s_p     = r_p / (r_p + r_e_new)

      call interp(s_gp, gama_mu_1, SDIV, s_p, gama_pole_h)
      call interp(s_gp, gama_mu_0, SDIV, s_e, gama_equator_h)
      gama_center_h = gama(1,1) ! hat

      call interp(s_gp,  rho_mu_1, SDIV, s_p, rho_pole_h)
      call interp(s_gp,  rho_mu_0, SDIV, s_e, rho_equator_h)
      rho_center_h = rho(1,1) ! hat
    
      grgr = gama_pole_h + rho_pole_h - gama_center_h - rho_center_h ! hat
      r_e_new_sq = 2.d0 * ( h_center - enthalpy_min ) / grgr
      r_e_new = sqrt( r_e_new_sq )
      if (r_e_new .ne. r_e_new .or. r_e_new/r_e_old > 2) then
        write(*,*) " "
        stop 'change in r_e is too dramatic; L72 in spin'
      endif

      ! Compute angular velocity Omega, only for rigid rotation
      if(r_ratio == 1.0) then
        Omega_c = 0.d0
        ww_equator_h = 0.d0
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

      ! Compute velocity, energy density and pressure
      do s = 1, SDIV
        sgp = s_gp(s)
        do m = 1, MDIV
          rsm = rho(s,m) ! hat
          if (r_ratio == 1.d0) then
            velocity_sq(s,m) = 0.d0
          else
            velocity_sq(s,m) = ( (Omega_c-ww(s,m))*(sgp/(1.d0-sgp)) * sin_theta(m) * exp(-rsm*r_e_new_sq) )**2
          endif

          if (velocity_sq(s,m) > 1.d0) velocity_sq(s,m) = 0.d0

          ! enthalpy_min is the assumed small value for the value at the pole
          enthalpy(s,m) = enthalpy_min + 5.d-1 * ( &
                r_e_new_sq * ( gama_pole_h + rho_pole_h - gama(s,m) - rsm ) - log(1.d0-velocity_sq(s,m))  )

          if (enthalpy(s,m) <= enthalpy_min .or. sgp > s_e) then
            enthalpy(s,m) = enthalpy_min
            pressure(s,m) = 0.d0
            energy(s,m)   = 0.d0
          else
            pressure(s,m) = p_at_h(enthalpy(s,m))
            energy(s,m)   = e_at_p(pressure(s,m))
          endif
          ! Rescale back metric potentials (except omega)
          rho  (s,m) = rho  (s,m) * r_e_new_sq
          gama (s,m) = gama (s,m) * r_e_new_sq
          alpha(s,m) = alpha(s,m) * r_e_new_sq
        enddo
      enddo

      ! Compute metric potentials
      S_metric_rho   = 0.d0
      S_metric_gama  = 0.d0
      S_metric_omega = 0.d0

      do s = 1, SDIV
        do m = 1, MDIV
          rsm   = rho     (s,m)
          gsm   = gama    (s,m)
          wwsm  = ww      (s,m) ! hat
          esm   = energy  (s,m)
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

          S_metric_rho(s,m) = e_gsm*( ea/2.d0 * (esm + psm) * s2 * (1.d0+v2sm) / (1.d0-v2sm) &
              + s2 * m1 * e_rsm**2 * ( (s1*d_ww_s)**2 + m1*d_ww_m**2 ) &
              + s1 * d_gama_s - mum * d_gama_m &
              + rsm / 2.d0 * ( ea*psm*s2 - s1*d_gama_s*(s1/2.d0*d_gama_s+1.d0) &
              - d_gama_m*(m1/2.d0*d_gama_m-mum)))

          S_metric_gama(s,m) = e_gsm*(ea*psm*s2 + gsm/2.d0*(ea*psm*s2 &
              - (s1*d_gama_s)**2/2.d0 - m1*d_gama_m**2/2.d0))

          S_metric_omega(s,m) = e_gsm*e_rsm*(-ea*(Omega_c-wwsm)*(esm+psm)*s2/(1.d0-v2sm) &
              + wwsm*(-0.5*ea*s2*(((1.d0+v2sm)*esm + 2.d0*v2sm*psm)/(1.0-v2sm)) &
              - s1 *(2*d_rho_s+0.5*d_gama_s) &
              + mum*(2*d_rho_m+0.5*d_gama_m) &
              + 0.25*s1**2*(4.d0*d_rho_s**2 - d_gama_s**2) &
              + 0.25*m1   *(4.d0*d_rho_m**2 - d_gama_m**2) &
              - m1*e_rsm**2*(sgp**4*d_ww_s**2 + s2*m1*d_ww_m**2)))
          if (S_metric_rho(s,m).ne.S_metric_rho(s,m) .or. S_metric_gama(s,m).ne.S_metric_gama(s,m) &
            .or. S_metric_omega(s,m).ne.S_metric_omega(s,m)) then 
            write(*,*) esm,psm,enthalpy(s,m),s_gp(s)/s_e,mu(m)
            stop "L184 in spin"
          endif
        enddo
      enddo

      sum_rho = 0.d0
      !--- Angular Integration
      n = 0
      do k = 1, SDIV
        do m = 1, MDIV-2, 2
          sum_rho = sum_rho + (DM/3.d0)* &
            (          P_2n(m  ,n+1) * S_metric_rho(k,m  ) &
              + 4.d0 * P_2n(m+1,n+1) * S_metric_rho(k,m+1) &
              +        P_2n(m+2,n+1) * S_metric_rho(k,m+2) )
        enddo
        D1_metric_rho  (n+1,k) = sum_rho 
        D1_metric_gama (n+1,k) = 0.d0
        D1_metric_omega(n+1,k) = 0.d0
      enddo

      sum_rho = 0.d0
      sum_gama = 0.d0
      sum_omega = 0.d0
      do n = 1, LMAX
        do k = 1, SDIV
          do m = 1, MDIV-2, 2
            sum_rho = sum_rho + (DM/3.d0) * &
                    (          P_2n(m  ,n+1) * S_metric_rho(k,  m) &
                      + 4.d0 * P_2n(m+1,n+1) * S_metric_rho(k,m+1) &
                      +        P_2n(m+2,n+1) * S_metric_rho(k,m+2) )
                       
            sum_gama = sum_gama + (DM/3.d0) * &
                    (          sin_2n_1_theta(m  ,n) * S_metric_gama(k,m  ) &
                      + 4.d0 * sin_2n_1_theta(m+1,n) * S_metric_gama(k,m+1) &
                      +        sin_2n_1_theta(m+2,n) * S_metric_gama(k,m+2) )
  
            sum_omega = sum_omega + (DM/3.d0) * &
                    (          sin_theta(m  ) * P1_2n_1(m,  n+1) * S_metric_omega(k,m  ) &
                      + 4.d0 * sin_theta(m+1) * P1_2n_1(m+1,n+1) * S_metric_omega(k,m+1) &
                      +        sin_theta(m+2) * P1_2n_1(m+2,n+1) * S_metric_omega(k,m+2) )
          enddo
          
          D1_metric_rho  (n+1,k) = sum_rho  
          D1_metric_gama (n+1,k) = sum_gama
          D1_metric_omega(n+1,k) = sum_omega
          sum_rho   = 0.d0
          sum_gama  = 0.d0
          sum_omega = 0.d0
        enddo
      enddo

      !--- RADIAL INTEGRATION
      n = 0
      do s = 1, SDIV
        do k = 1, SDIV-2, 2
          sum_rho = sum_rho + (DS/3.d0) * &
                  (          f_rho(s,n+1,k  ) * D1_metric_rho(n+1,k  ) &
                    + 4.d0 * f_rho(s,n+1,k+1) * D1_metric_rho(n+1,k+1) &
                    +        f_rho(s,n+1,k+2) * D1_metric_rho(n+1,k+2) )
        enddo
        
        !write(*,"(es15.6)") sum_rho,D1_metric_rho(n+1,k),D1_metric_rho(n+1,k+1),D1_metric_rho(n+1,k+2); stop
        D2_metric_rho  (s,n+1) = sum_rho
        D2_metric_gama (s,n+1) = 0.d0
        D2_metric_omega(s,n+1) = 0.d0
        sum_rho = 0.d0
      enddo

      do s = 1, SDIV
        do n = 1, LMAX
          do k = 1, SDIV-2, 2
            sum_rho = sum_rho + (DS/3.d0) * &
                    (          f_rho(s,n+1,k  ) * D1_metric_rho(n+1,k  ) &
                      + 4.d0 * f_rho(s,n+1,k+1) * D1_metric_rho(n+1,k+1) &
                      +        f_rho(s,n+1,k+2) * D1_metric_rho(n+1,k+2) )
            sum_gama = sum_gama + (DS/3.d0) * &
                    (          f_gama(s,n+1,k  ) * D1_metric_gama(n+1,k  ) &
                      + 4.d0 * f_gama(s,n+1,k+1) * D1_metric_gama(n+1,k+1) &
                      +        f_gama(s,n+1,k+2) * D1_metric_gama(n+1,k+2) )
            if ( k < s .and. k+2 <= s ) then
              sum_omega = sum_omega + (DS/3.d0) * &
                    (          f_rho(s,n+1,k  ) * D1_metric_omega(n+1,k  ) &
                      + 4.d0 * f_rho(s,n+1,k+1) * D1_metric_omega(n+1,k+1) &
                      +        f_rho(s,n+1,k+2) * D1_metric_omega(n+1,k+2) )
            else
              if (k>=s) then
                sum_omega = sum_omega + (DS/3.d0) * &
                    (          f_gama(s,n+1,k  ) * D1_metric_omega(n+1,k  ) &
                      + 4.d0 * f_gama(s,n+1,k+1) * D1_metric_omega(n+1,k+1) &
                      +        f_gama(s,n+1,k+2) * D1_metric_omega(n+1,k+2) ) 
              else
                sum_omega = sum_omega + (DS/3.d0) * &
                    (         f_rho (s,n+1,k  ) * D1_metric_omega(n+1,k  ) &
                      + 4.d0* f_rho (s,n+1,k+1) * D1_metric_omega(n+1,k+1) &
                      +       f_gama(s,n+1,k+2) * D1_metric_omega(n+1,k+2) )
              endif
            endif
          enddo
          D2_metric_rho(s,n+1)   = sum_rho
          D2_metric_gama(s,n+1)  = sum_gama
          D2_metric_omega(s,n+1) = sum_omega
          sum_rho   = 0.d0
          sum_gama  = 0.d0
          sum_omega = 0.d0
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
          e_gsm = exp(-gsm/2.d0)
          e_rsm = exp(rsm)
          temp1 = sin_theta(m)

          sum_rho = sum_rho - e_gsm * P_2n(m,1) * D2_metric_rho(s,1)

          do n = 1, LMAX
            sum_rho = sum_rho - e_gsm * P_2n(m,n+1) * D2_metric_rho(s,n+1)
            if (m == MDIV) then
              sum_gama  = sum_gama  - 2.d0 / pi * e_gsm * D2_metric_gama(s,n+1)
              sum_omega = sum_omega + e_rsm*e_gsm*D2_metric_omega(s,n+1)/2.d0
            else
              sum_gama  = sum_gama - 2.d0/pi*e_gsm*(sin_2n_1_theta(m,n)/((2.d0*n-1)*temp1)) &
                * D2_metric_gama(s,n+1)
              sum_omega = sum_omega - e_rsm*e_gsm*(P1_2n_1(m,n+1)/(2.d0*n*(2.d0*n-1)*temp1)) &
                * D2_metric_omega(s,n+1)
            endif
          enddo
          rho (s,m) = sum_rho    !  rsm + cf*(sum_rho  - rsm)
          gama(s,m) = sum_gama   !  gsm + cf*(sum_gama - gsm)
          ww  (s,m) = sum_omega  ! wwsm + cf*(sum_omega-wwsm)
          sum_omega = 0.d0
          sum_rho   = 0.d0
          sum_gama  = 0.d0
        enddo
      enddo

          
      ! check for divergence
      if (abs(rho(2,1))>100.d0 .or. abs(gama(2,1))>300.d0 .or. abs(ww(2,1))>100.d0) then
        write(*,"(3es18.9)") rho(2,1), gama(2,1), ww(2,1); stop "Line 300 in spin"
      endif
            
      if (r_ratio == 1.d0) then
        do s = 1, SDIV
          do m = 1, MDIV
            rho(s,m) = rho(s,1)
            gama(s,m)= gama(s,1)
            ww(s,m)  = 0.d0
          enddo
        enddo
      endif
        
    ! --- the fourth equation
    
    ! compute first order derivatives of gama

    do s = 1, SDIV; do m = 1, MDIV
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

          temp8 = m1 * exp( -2.d0 * rho(s,m) )
                
          da_dm(s,m) = - (d_rho_m+d_gama_m)/2.d0 &
            - temp2*( (temp3 - d_gama_mm - temp4)*(-mum+m1*d_gama_m)/2.d0 &
            + temp5/4.d0 - temp6  + temp7 + temp8*temp1/4.d0 )
        enddo
      enddo
    endif

    do s = 1, SDIV
      alpha(s,1) = 0.0
      do m = 1, MDIV-1
        alpha(s,m+1) = alpha(s,m) + dm * ( da_dm(s,m+1) + da_dm(s,m) ) / 2.d0
      enddo
    enddo

    do s = 1, SDIV
      do m = 1, MDIV
        alpha(s,m) = alpha(s,m) - alpha(s,MDIV) + ( gama(s,MDIV)-rho(s,MDIV) )/2.d0
        if(alpha(s,m).ge.300.0) then
          write(*,*) s,m
          stop "L392, alpha fails"
        endif
        ww(s,m) = ww(s,m) / r_e_new
      enddo
    enddo

    if (smax == 1.0) then
      alpha(SDIV,:) = 0.0
    endif
    dif = abs(r_e_old-r_e_new)/r_e_new
    n_of_it = n_of_it + 1
    !write(*,"(3es15.6)") r_e_old,r_e_new,dif
    if (n_of_it == 100) stop "Cannot converge; L404 in spin"
  enddo
! --- End of iteration 

  ! compute omega
  omg(:,:) = Omega_c / r_e_new
  Omega_c  = Omega_c / r_e_new
  r_e      = r_e_new

  if (output) then
    write(fil1,"(f4.2)") ang_mom
    write(fil2,"(f4.2)") mass_0/MSUN
    open(98,file="./Cont/J"//trim(adjustl(fil1))//"_Mb"//trim(adjustl(fil2))//".dat")
    !open(98,file="./Cont/check2D.dat")
#if defined(matlab)
#else
    write(98,"(2i5,99es27.17)") SDIV, MDIV, r_e*sqrt(KAPPA)/1.d5, energy(1,1)/(C*C*KSCALE), r_ratio
#endif
    do s = 1, SDIV
      do m = 1, MDIV
        ! r, \theta, \apha, \gamma, \rho, \omega, \phi, \varepsilon, \rho_0, p
        if (enthalpy(s,m) > enthalpy_min) then 
          rho_0 = n0_at_e( energy(s,m) ) * MB
        else 
          rho_0 = 0.d0
        endif
        write(98,"(99es27.17)") s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), ww(s,m) * (C/sqrt(kappa)), & ! 1-6
        pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), enthalpy(s,m), rho_0, & ! 7-10
        velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)) ! 11-12
      enddo
    enddo
    close(98)
  endif


end subroutine spin


