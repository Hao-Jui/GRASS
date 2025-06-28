subroutine shoot_v2
#include "option_macro.h"
  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it
  real(8) :: er, rho0, ee

  r_ratio  = 9.340690294E-01
  output   = .false.

#if defined(restart)
  call restart_read
  e_center = 2.54d15
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
  output = .true.
#else
  e_center = 1.616234914E+15
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center) 
  call sphere
#endif
  write(*,*) " "
  
  it = 1
  er = 1.d99
  do while ( er > 1.d-4 .and. .not.output )

    call spin
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    if ( mod(it,10) == 0 ) then
    write(*,*) " ===================================="
    write(*,"(A10,i6)")         " iter :", it/10
    write(*,"(A10,es18.9,A10)") "rho_c :", rho0*MB,"g/cm^3"
    write(*,"(A10,es18.9,A10,es18.9)") "e_c   :", ee/(C * C * KSCALE),"g/cm^3", ee/(C * C * KSCALE)*rho_uni
    write(*,"(A10,es18.9)")     "rp/re :", r_ratio
    write(*,"(A10,es18.9,A4)")  "OMG_c :", omega_c/2.d0/pi* (C/sqrt(kappa)),"Hz"
    write(*,"(A10,es18.9,A4)")  "Mass  :", Mass/MSUN,"M_o"
    write(*,"(A10,es18.9,A4)")  "M_0   :", Mass_0/MSUN,"M_o"
    write(*,"(A10,es18.9)")     "J     :", ang_mom
    write(*,"(A10,es18.9)")     "chi   :", chi
    write(*,"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
    write(*,"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
    write(*,"(A10,3es18.9)")    "er    :", er
    write(*,*) " ===================================="
    write(*,*) " "
    endif
  
    call update_v2(h_center, r_ratio, er)
    it = it + 1
    if (it==1000) stop "Iteration may never converge."
  enddo
  
  !write(*,*) "no 2D output "
  !stop " "
  output = .true.
  call spin
  call mass_radius
  rho0 = n0_at_h(h_center)
  ee   = e_at_h (h_center)
  write(*,*) " "
  write(*,"(10es18.9)") ee/(C * C * KSCALE), r_ratio, r_e

  open(221,file="./Cont/properties.dat")
  write(*,*) " ===================================="
  write(*,*) "              Converged              "
  do it = 1, 2
    write(6+215*(it-1),"(A10,es18.9,A10)") "rho_c :", rho0*MB,"g/cm^3"
    write(6+215*(it-1),"(A10,es18.9,A10,es18.9)") "e_c   :", ee/(C * C * KSCALE),"g/cm^3", ee/(C * C * KSCALE)*rho_uni
    write(6+215*(it-1),"(A10,es18.9)")     "rp/re :", r_ratio
    write(6+215*(it-1),"(A10,es18.9,A4)")  "OMG_c :", omega_c/2.d0/pi* (C/sqrt(kappa)),"Hz"
    write(6+215*(it-1),"(A10,es18.9,A4)")  "Mass  :", Mass/MSUN,"M_o"
    write(6+215*(it-1),"(A10,es18.9,A4)")  "M_0   :", Mass_0/MSUN,"M_o"
    write(6+215*(it-1),"(A10,es18.9)")     "J     :", ang_mom
    write(6+215*(it-1),"(A10,es18.9)")     "chi   :", chi
    write(6+215*(it-1),"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
    write(6+215*(it-1),"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
  enddo
  write(*,*) " ===================================="
  close(221)
  write(*,*) " "
  write(*,*) "Completed!"

end subroutine shoot_v2


subroutine update_v2(hc, rep, er)
  use para_mod
  implicit none
  real(8), intent(inout) :: hc, rep
  real(8), intent(out) :: er
  real(8) :: new

  er = abs(mass_0/MSUN-Mb_goal) + abs(ang_mom-J_goal)
  if ( er < 1.d-5 ) return

  if (abs(mass_0/MSUN - Mb_goal) > abs(ang_mom - J_goal) ) then
    new = hc + (Mb_goal - mass_0/MSUN) * 1.d-1
    write(*,"(A10,es15.6,A10,es15.6)") "hc:",hc,"-->",new
    hc = new
  else
    new = min( 1.d0, rep + (ang_mom-J_goal)*1.d-1 )
    write(*,"(A10,es15.6,A10,es15.6)") "rp/re:",rep,"-->", new
    rep = new
  endif

end subroutine update_v2
