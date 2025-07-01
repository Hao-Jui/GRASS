subroutine J_seq

  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  real(8) :: er,ee,rho0
  integer :: it, cnt

  r_ratio  = 1.d0
  output   = .false.

  e_center = 6.d14
  e_center = e_center * C * C * KSCALE
  ee       = e_center
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center) 
  call sphere
  
  open(23,file="./Cont/J1.8_seq_MPA1.dat")
  
  write(*,*) " "
  do cnt = 1, 1000
    if ( ee/(C * C * KSCALE) < 2.d15 ) then
      h_center = h_center * 1.1d0
    else
      h_center = h_center * 1.005d0
    endif
    it = 1
    er = 1.d99
    do while ( er > 1.d-5 )

      call spin
      call mass_radius

      rho0 = n0_at_h(h_center)
      ee   = e_at_h (h_center)

      call update_J_seq(r_ratio, er)
      it = it + 1
      if ( it == 100 ) exit
    enddo

    write(*,"(2i5,14es18.9)") cnt, it, ee/(C * C * KSCALE), mass/MSUN
    write(23,"(14es18.9)") ee/(C * C * KSCALE), r_ratio, mass_0/MSUN, mass/MSUN, &
                r_e, omega_c/2.d0/pi* (C/sqrt(kappa)), ang_mom
    if (ee/(C * C * KSCALE) > 3.d15) exit
  enddo
  close(23)
  write(*,*) " "
  write(*,*) "Completed!"
  
end subroutine J_seq

subroutine update_J_seq(rep, er)
  use para_mod
  implicit none
  real(8), intent(inout) :: rep
  real(8), intent(out) :: er
  real(8) :: new

  new = min( 1.d0, rep + (ang_mom-J_goal)*3.d-2 )
  !write(*,"(A10,es15.6,A10,es15.6)") "rp/re:",rep,"-->", new
  rep = new
  
  er = abs(ang_mom-J_goal)

end subroutine update_J_seq

