subroutine fix_J(J_goal)

	use para_mod
	implicit none
	real(8), intent(in) :: J_goal
	real(8) :: p_at_e, h_at_p, n0_at_e, n0_at_h, e_at_h
	real(8) :: mass, mass_0, chi
	real(8) :: Ju, Jd
	real(8) :: rup, rdo
	integer :: it
	real(8) :: er,rho0,ee


	e_center = 9.d14
	e_center = e_center * C * C * KSCALE

	p_center = p_at_e(e_center)
	h_center = h_at_p(p_center)
	call sphere

	rup  = 1.d0
	rdo  = 9.d-1
	output   = .false.

	write(*,*) " "
	write(*,*) "Initial guess:"

	r_ratio = rup; call spin; call mass_radius(mass,mass_0,chi); Ju = ang_mom
	write(*,*) " "
	r_ratio = rdo; call spin; call mass_radius(mass,mass_0,chi); Jd = ang_mom

	if ( (Ju-J_goal)*(Jd-J_goal) > 0.d0 ) stop "rp not bracketed"

	it = 1; er = 1.d99
	write(*,*) " "
	write(*,*) "Iterating:"
	do while (it < 10000 .and. abs(er) > 1.d-6)
		write(*,*) " ===================================="
		write(*,"(A10,i5)") "iter  :",it
		r_ratio = (rup + rdo) / 2.d0
		call spin
		call mass_radius(mass,mass_0,chi)
		er = ang_mom-J_goal
		if ( er*(Ju-J_goal) > 0 ) then 
			rup = r_ratio
		else
			rdo = r_ratio
		endif
		if (rup == rdo) stop "line 41"

		rho0 = n0_at_h(h_center)
		ee   = e_at_h (h_center)
		write(*,"(A10,es18.9,A10)") "rho_c   :", rho0*MB,"g/cm^3"
		write(*,"(A10,es21.11,A10)") "  e_c   :", ee/(C * C * KSCALE),"g/cm^3"
		write(*,"(A10,es21.11)")  "r_ratio :", r_ratio 
		write(*,"(A10,es18.9)")  "er  :", er 
		write(*,*) " ===================================="
		write(*,*) " "
		it = it + 1
	enddo

	output = .true.
	call spin



end subroutine fix_J
