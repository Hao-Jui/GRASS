subroutine shoot_v2
#include "option_macro.h"
  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it
  real(8) :: er, rho0, ee

  r_ratio  = 7.632149492E-01
  output   = .false.

#if defined(restart)
  call refine_read
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#else
  e_center = 7.208248462E+14
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center) 
  call sphere
#endif
  write(*,*) " "
  
  it = 1
  er = 1.d99
  do while ( er > accuracy .and. .not.output )

    call uryu
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
    write(*,"(A10,es18.9)")     "T/W   :", T_kin/abs(Mass - Mass_0 - Mass_p + T_kin)
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
  call uryu
  call mass_radius
  rho0 = n0_at_h(h_center)
  ee   = e_at_h (h_center)
  write(*,*) " "

  open(221,file="./Cont/properties.dat")
  write(*,*) " ===================================="
  write(*,*) "              Converged              "
  do it = 1, 2
    write(6+215*(it-1),"(A10,es18.9,A10)") "rho_c :", rho0*MB,"g/cm^3"
    write(6+215*(it-1),"(A10,es18.9,A10,es18.9)") "e_c   :", ee/(C * C * KSCALE),"g/cm^3", ee/(C * C * KSCALE)*rho_uni
    write(6+215*(it-1),"(A10,es18.9)")     "rp/re :", r_ratio
    write(6+215*(it-1),"(A10,es18.9,A4)")  "OMG_c :", omega_c/2.d0/pi* (C/sqrt(kappa)),"Hz"
    write(6+215*(it-1),"(A10,es18.9,A4)")  "Oc/Oe :", omega_c/Omega_e
    write(6+215*(it-1),"(A10,es18.9,A4)")  "Mass  :", Mass/MSUN,"M_o"
    write(6+215*(it-1),"(A10,es18.9,A4)")  "M_0   :", Mass_0/MSUN,"M_o"
    write(6+215*(it-1),"(A10,es18.9)")     "J     :", ang_mom
    write(6+215*(it-1),"(A10,es18.9)")     "chi   :", chi
    write(6+215*(it-1),"(A10,es18.9)")     "T/W   :", T_kin/abs(Mass_p - Mass + T_kin)
    write(6+215*(it-1),"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
    write(6+215*(it-1),"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
  enddo
  write(*,*) " "
  write(*,*) "In code unit:"
  write(*,"(A10,es18.9)") "h_c", h_center, "r_e", r_e, "Fmax", Fmax_h, "Omega_e", Omega_e * r_e
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
  real(8) :: new, deviA, deviB

  deviA = mass_0/MSUN - Mb_goal
  deviB = ang_mom - J_goal


  er = abs(deviA) + abs(deviB)
  if ( er < accuracy ) return

  if ( er > 1.d-1 ) then
    if ( abs(deviA) > abs(deviB) ) then
      new = hc - deviA * 3.d-2
      write(*,"(A10,es15.6,A10,es15.6)") "hc:",hc,"-->",new
      hc = new
    else
      new = min( 1.d0, rep + deviB*3.d-2 )
      write(*,"(A10,es15.6,A10,es15.6)") "rp/re:",rep,"-->", new
      rep = new
    endif
  else
    if (abs(deviA) > abs(deviB) ) then
      new = hc - deviA * 1.d-1
      write(*,"(A10,es15.6,A10,es15.6)") "hc:",hc,"-->",new
      hc = new
    else
      new = min( 1.d0, rep + deviB*1.d-1 )
      write(*,"(A10,es15.6,A10,es15.6)") "rp/re:",rep,"-->", new
      rep = new
    endif
  endif

end subroutine update_v2
