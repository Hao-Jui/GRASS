! ********************************************* !
!
!         solving the stellar structure         !
!
! ********************************************* !
subroutine sphere

  use para_mod
    implicit none
    integer s,m
    real(8) r_is_s, r_is_final, r_final, m_final, &
           lambda_s, nu_s, e_s, gama_eq, rho_eq
    real(8), dimension(SDIV) :: gama_mu_0, rho_mu_0
    real(8), dimension(RDIV) :: r_is_gp, lambda_gp, nu_gp, e_d_gp

    write(*,*) " "
    write(*,*) "Configurating spherical guess ..."
    write(*,"(4A15)") "|      runs","r_is (km)","r (km)","m (M_o)"
!    write(*,"(A15,1P3E15.2)") "|  Real TOV",11.4484/sqrt(1-2*2.0572/11.4484), &
!        11.4484,2.0572/1.47664


    do s = 1,3
      call TOV(s, r_is_gp, lambda_gp, nu_gp, e_d_gp, r_is_final, r_final, m_final)
      write(*,"(A5,i10,3f15.5)") "|", &
          s, r_is_final*sqrt(KAPPA)/1.d5, r_final*sqrt(KAPPA)/1.d5, &
          m_final*sqrt(KAPPA)*C*C/G/MSUN
    enddo

  ! map the static NS to 2D data
    do s = 1, SDIV
      r_is_s = r_is_final * (s_gp(s)/(1-s_gp(s)))
      if (r_is_s <= r_is_final) then
        call interp(r_is_gp, lambda_gp, RDIV, r_is_s, lambda_s)
        call interp(r_is_gp,     nu_gp, RDIV, r_is_s,     nu_s)
        call interp(r_is_gp,    e_d_gp, RDIV, r_is_s,      e_s)
      else
        e_s      = 1.d-3*(C*C*KSCALE)
        lambda_s = 2.d0*log(1.d0+m_final/(2.0*r_is_s))
        nu_s     = log((1.d0-m_final/(2.d0*r_is_s))/(1.d0+m_final/(2*r_is_s)))
      endif
        sphi (s,:) = ( 1.d0 - exp(nu_s) ) /1.d1 * exp(-sqrt(mphi_r)*r_is_s) 
        rho  (s,:) = nu_s-lambda_s
        gama (s,:) = lambda_s+nu_s
        alpha(s,:) = (lambda_s-nu_s)/2.d0
        energy(s,:)= e_s
        rho_mu_0 (s) = nu_s-lambda_s
        gama_mu_0(s) = lambda_s+nu_s
    enddo
    ww(:,:) = 0.d0
    omg(:,:)= 0.d0
    
  ! since s=0.5 is not necessary at the s-grid, 
  ! we need to find gama and rho at the equator
    call interp(s_gp, gama_mu_0, SDIV, s_e, gama_eq)
    call interp(s_gp,  rho_mu_0, SDIV, s_e,  rho_eq)
    
  ! r_e is roughly r_is_final
    r_e = r_final * exp(0.5*(rho_eq-gama_eq))


  end subroutine sphere

! *****************************
!              TOV
! *****************************
subroutine TOV(i_check, r_is_gp, lambda_gp, nu_gp, e_d_gp, &
    r_is_final, r_final, m_final)

    use para_mod
    integer n_nearest, i_check, i
    real(8),intent(inout):: r_is_final
    real(8),intent (out) :: r_final,m_final
    real(8),dimension(RDIV) :: r_is_gp, lambda_gp, nu_gp, e_d_gp
    real(8),dimension(RDIV) :: r_gp, m_gp
    real(8) r, r_is, r_is_est, r_is_check, dr_is_save, &
            e_d, p, h, m, nu_s, hh, rho_0, &
            a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4, &
            k_rescale
    real(8) dm_dr_is,dp_dr_is,dr_dr_is,h_at_p,p_at_e,e_at_p,n0_at_e

    ! use estimate of r to set the step size h
    if (i_check == 1) then
      r_is_est = 1.5d6/sqrt(kappa) ! 15 km
      h = r_is_est/100
    else
      r_is_est= r_is_final
      h = r_is_est/10000
      dr_is_save = r_is_final/rdiv
      r_is_check = dr_is_save
    endif

    r_is = 0.0                            ! initial isotropic radius
    r = 0.0                               ! initial radius
    m = 0.0                               ! initial mass
    p = p_center                          ! initial pressure (code unit)

    r_is_gp(1) = 0.0
    r_gp(1) = 0.0
    m_gp(1) = 0.0
    lambda_gp(1) = 0.0
    e_d_gp(1) = e_center

    i = 2
    !write(*,"(3es15.6)") r,m,p    
    do while(p >= p_surface)
      e_d = e_at_p(p)

      if (i_check == 3 .and. r_is > r_is_check .and. i < rdiv) then
        r_is_gp(i) = r_is
        r_gp(i) = r
        m_gp(i) = m
        e_d_gp(i) = e_d
        i = i+1
        r_is_check = r_is_check + dr_is_save
      endif
      r_is_final = r_is
      r_final = r
      m_final = m

      a1 = dr_dr_is(r_is, r, m)
      b1 = dm_dr_is(r_is, r, m, p)
      c1 = dp_dr_is(r_is, r, m, p)

      a2 = dr_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0)
      b2 = dm_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0)
      c2 = dp_dr_is(r_is+h/2.0, r+h*a1/2.0, m+h*b1/2.0, p+h*c1/2.0)

      a3 = dr_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0)
      b3 = dm_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0)
      c3 = dp_dr_is(r_is+h/2.0, r+h*a2/2.0, m+h*b2/2.0, p+h*c2/2.0)

      a4 = dr_dr_is(r_is+h, r+h*a3, m+h*b3)
      b4 = dm_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3)
      c4 = dp_dr_is(r_is+h, r+h*a3, m+h*b3, p+h*c3)

      r = r + (h/6.0)*(a1+2*a2+2*a3+a4)
      m = m + (h/6.0)*(b1+2*b2+2*b3+b4)
      p = p + (h/6.0)*(c1+2*c2+2*c3+c4)

      r_is = r_is+h
      !write(*,"(3es15.6)") r,m,p
    enddo
    e_d_gp (rdiv) = 0.d0
    r_is_gp(rdiv) = r_is_final
    r_gp   (rdiv) = r_final
    m_gp   (rdiv) = m_final

    ! Rescale r_is and compute lambda, nu
    
    if (i_check == 3) then
      k_rescale = 0.5*(r_final/r_is_final)* &
        (1.0-m_final/r_final + sqrt(1.0-2.0*m_final/r_final) )
      !write(*,"(es18.9)") k_rescale
      r_is_final = r_is_final * k_rescale
      ! lapse = e^{2 nu_s}
      nu_s = log( (1.0-m_final/(2.0*r_is_final))/ &
        (1.0+m_final/(2.0*r_is_final)) )
          
      open(988,file="./Cont/checkTOV.dat")
      do i = 1, rdiv
        r_is_gp(i) = r_is_gp(i)*k_rescale
        ! psi^4 = e^{2 lambda_gp}
        if (i == 1) then
          lambda_gp(1) = log(1.0/k_rescale)
        else
          lambda_gp(i) = log(r_gp(i)/r_is_gp(i))
        endif
        
        if(e_d_gp(i) < e_surface) then
          hh = 0.d0
        else
          p = p_at_e(e_d_gp(i))
          hh = h_at_p(p)
          rho_0 = n0_at_e(e_d_gp(i))
        endif
        nu_gp(i) = nu_s - hh
      
        if ( e_d_gp(i)/(C*C*KSCALE) > 1.d16 ) stop " bug, L177 Sphere"
        
        write(988,"(99es18.9)") r_is_gp(i)*sqrt(KAPPA)/1.d5, &
              r_gp(i)*sqrt(KAPPA)/1.d5, m_gp(i)*sqrt(KAPPA)/1.d5, e_d_gp(i)/(C*C*KSCALE), &
              p/KSCALE, rho_0*MB, nu_gp(i), lambda_gp(i)
      enddo
      close(988)
      nu_gp(rdiv)=nu_s
    endif

end subroutine TOV

real(8) function dm_dr_is(r_is,r,m,p)

	use para_mod, only : p_surface,tov_rmin,e_center,pi
	implicit none
	real(8) r_is,r,m,p, e_d,e_at_p

	if(p < p_surface) then
		e_d = 0.d0
	else
    	e_d = e_at_p(p)
	endif
	if(r_is < tov_rmin) then
    	dm_dr_is=4*pi * e_center * r**2 * (1+4*pi*e_center*r**2/3)
	else
    	dm_dr_is=4*pi * e_d * r**3 * sqrt(1-2*m/r) / r_is
	endif

end function dm_dr_is

real(8) function dp_dr_is(r_is,r,m,p)

	use para_mod, only : p_surface,tov_rmin,e_center,pi
	implicit none
	real(8) r_is,r,m,p, e_d,e_at_p

	if(p<p_surface) then
		e_d = 0.d0
	else
    	e_d = e_at_p(p)
	endif
	if(r_is < tov_rmin) then
    	dp_dr_is=-4*pi*(e_center+p)*(e_center+3*p)*r*(1+4*e_center*r**2/3)/3
	else
    	dp_dr_is=-(e_d+p)*(m+4*pi*r**3*p)/(r*r_is*sqrt(1.0-2.0*m/r))
	endif

end function dp_dr_is

real(8) function dr_dr_is(r_is,r,m)

	use para_mod, only : tov_rmin
	implicit none
	real(8) r_is,r,m

	if(r_is<tov_rmin) then
    	dr_dr_is=1.d0
	else
    	dr_dr_is=(r/r_is)*sqrt(1-2*m/r)
	endif

end function dr_dr_is


