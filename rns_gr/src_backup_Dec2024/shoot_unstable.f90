subroutine shoot_unstable
#include "option_macro.h"
  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it
  real(8) :: er,rho0,ee,e_c,e_c_uni

  r_ratio  = 8.5d-1
  output   = .false.

  e_c      = 2.65d15
  e_center = e_c
  
  if ( e_c > 2.d15 ) then
#if defined(restart)
    call restart_read
    e_center = e_center * C * C * KSCALE
#else
    e_center = e_center * C * C * KSCALE / 2.d0
#endif
    e_c_uni  = e_center
    p_center = p_at_e(e_center)
    h_center = h_at_p(p_center)
  else
    e_center = e_center * C * C * KSCALE
    p_center = p_at_e(e_center)
    h_center = h_at_p(p_center)
  endif  
  call sphere
  
  write(*,*) " "
  it = 1
  er = 1.d99
  do while (e_center < e_c* C * C * KSCALE .or. er > 1.d-5)

    if ( e_c > 2.d15 ) then
      e_center = e_c_uni * min( 2.d0, 1.005d0**(it-1) )
      p_center = p_at_e(e_center)
      h_center = h_at_p(p_center)
    endif

    if ( mod(it,20) == 0 ) output = .true.

    call spin
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)
    
    call update_unstable(r_ratio, er, J_goal )

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
    it = it + 1
    output = .false.
  enddo
  
  !write(*,*) "Completed!"
  !write(*,*) "no 2D output "
  !stop " "
  output = .true.
  call spin
  call mass_radius

  open (221,file="./Cont/properties.dat")
    write(221,"(A10,i6)")         " iter :", it
    write(221,"(A10,es18.9,A10)") "rho_c :", rho0*MB,"g/cm^3"
    write(221,"(A10,es18.9,A10,es18.9)") "e_c   :", ee/(C * C * KSCALE),"g/cm^3", ee/(C * C * KSCALE)*rho_uni
    write(221,"(A10,es18.9)")     "rp/re :", r_ratio
    write(221,"(A10,es18.9,A4)")  "OMG_c :", omega_c/2.d0/pi* (C/sqrt(kappa)),"Hz"
    write(221,"(A10,es18.9,A4)")  "Mass  :", Mass/MSUN,"M_o"
    write(221,"(A10,es18.9,A4)")  "M_0   :", Mass_0/MSUN,"M_o"
    write(221,"(A10,es18.9)")     "J     :", ang_mom
    write(221,"(A10,es18.9)")     "chi   :", chi
    write(221,"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
    write(221,"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
  close(221)

end subroutine shoot_unstable


subroutine update_unstable(rep, er, goal_B)
  use para_mod
  implicit none
  real(8), intent(inout) :: rep
  real(8), intent(out) :: er
  real(8), intent(in) :: goal_B
  real(8) :: new


  new = min( 1.d0, rep + (ang_mom-goal_B)*3.d-2 )
  write(*,"(A10,es15.6,A10,es15.6)") "rp/re:",rep,"-->", new
  rep = new

  er = abs(ang_mom - goal_B)

end subroutine update_unstable
