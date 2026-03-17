subroutine set_disk(r_eq)
  use toolkit_mod, only: interp
  use para_mod, only: SDIV, MDIV, res, s_pwr, &
                      s_gp, mu, s_inner, j_disk, &
                      rho, gama, ww, omg, &
                      enthalpy, enthalpy_min, pressure, energy, velocity_sq
  implicit none
  real(8), intent(in) :: r_eq
  real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0, ww_mu_0
  real(8) :: rho_in, gama_in, ww_in, w_0, ww2, hh, m1
  real(8) :: tmp, u_phi
  real(8) :: r_h, r_in
  integer :: s, m
  logical :: set = .false.
  real(8) :: p_at_h, e_at_h

  rho_mu_0 (:) = rho (:,1)
  gama_mu_0(:) = gama(:,1)
  ww_mu_0  (:) = ww  (:,1)

  call interp(s_gp,  rho_mu_0, SDIV, s_inner,  rho_in)
  call interp(s_gp, gama_mu_0, SDIV, s_inner, gama_in)
  call interp(s_gp,   ww_mu_0, SDIV, s_inner,   ww_in)

  ! > equating u_phi at s_inner with j-disk
  r_in = ( s_inner/(1.d0-s_inner) )**s_pwr
  ww2  = (r_eq / j_disk)**2 * exp(rho_in+gama_in) + exp(2.d0*rho_in) / r_in**2
  w_0  = ww_in - sqrt( ww2 ) ! hat; minus sign may be changed to plus

  do s = res+1, SDIV
    if ( s_gp(s) > s_inner ) then
      r_h = ( s_gp(s)/(1.d0-s_gp(s)) )**s_pwr
      do m = 1, MDIV
          m1 = 1.d0 - mu(m)**2
          ! u_phi here is actually u_phi / r_e
          tmp = (ww(s,m)-w_0)**2 / exp(rho(s,m)+gama(s,m)) - exp(rho(s,m)-gama(s,m)) / m1 / r_h**2
          hh  = tmp / ( (ww_in-w_0)**2 / exp(rho_in+gama_in) - exp(rho_in-gama_in) / r_in**2 )
          enthalpy (s,m) = log( sqrt(hh) )

          if ( (enthalpy(s,m) < 0.d0) .or. (enthalpy(s,m).ne.enthalpy(s,m)) ) then
            enthalpy(s,m) = enthalpy_min
          endif

          if ( enthalpy(s,m) <= enthalpy_min ) then
            pressure(s,m) = 0.d0
            energy  (s,m) = 0.d0
            velocity_sq(s,m) = 0.d0
            omg(s,m) = 0.d0
          else
            u_phi = 1.d0 / sqrt( tmp )
            tmp   = exp(gama(s,m)-rho(s,m)) * r_h**2 * m1
            velocity_sq(s,m) = u_phi**2 / ( tmp + u_phi**2 )
            
            if ( velocity_sq(s,m) > 1.d0 .or. velocity_sq(s,m) < 0.d0 ) velocity_sq(s,m) = 0.d0

            omg(s,m) = ww(s,m) + exp(rho(s,m)) * sqrt( velocity_sq(s,m) / m1 ) / r_h
            pressure(s,m) = p_at_h( enthalpy(s,m) )
            energy  (s,m) = e_at_h( enthalpy(s,m) )
            !if ( m==1 ) write(*,"(10es15.6)") u_phi, velocity_sq(s,1), omg(s,1)
          endif
      enddo
      if ( enthalpy(s,1) > enthalpy_min ) then 
          set = .true.
          !write(*,"(10es15.6)") u_phi, velocity_sq(s,1), ww(s,1)
      endif
    endif
  enddo

  !if ( set .and. enthalpy(SDIV,1) == enthalpy_min ) write(*,*) "Disk bestowed!"


end subroutine set_disk

