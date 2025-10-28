subroutine MRcurve
#include "option_macro.h"
  use para_mod
  use rotation_dispatch, only: call_rotation_solver
  implicit none
  real(8), external :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it, unit
  real(8) :: rho0, ee, sphi_max
  real(8) :: B_ini = 16.d0
  character(16) :: fil1, fil2

  call initialize_starting_model(p_at_e, h_at_p)

  mphi_r  = (mphi_goal/l_uni)**2 * KAPPA / 1.d10
  B_coup  = B_goal
  rho0    = n0_at_h(h_center)
  
  it = 1
  output = .true.
  do while ( rho0*MB > 1.d14 )
    call call_rotation_solver()
    output = .false.
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)
    sphi_max = maxval(sphi(:,1)) * sqrt(B_coup) 
    
    if ( mod(it,10) == 0 ) call check_point()

    write(fil1,"(f10.1)") mphi_goal
    write(fil2,"(f12.0)") B_goal
    open(newunit=unit,file="/Users/horay/Data4Projects/crazy/Seq/"//trim(adjustl(eos_file))// &
            "/1d_mphi"//trim(adjustl(fil1))//"_B"//trim(adjustl(fil2))//"dat",access='append')
    write(unit,"(99es18.9e3)") ee/(C * C * KSCALE), rho0*MB, h_center, &
        Mass/MSUN, r_circ/1e5, Mass_0/MSUN, &
        sphi_c, sphi_max, T_kin/abs(Mass_p - Mass + T_kin)
    close(unit)
    it = it + 1
    h_center = h_center - 5.d-3
  enddo 

  write(*,*) " "
  write(*,*) "Completed!"

contains
  subroutine check_point
    write(*,*) " ===================================="
    write(*,"(A10,i6)")         " iter :", it
    write(*,"(A10,2es18.9)")     " sphi :", sphi_c, sphi_max
    write(*,"(A10,es18.9,A10)") "rho_c :", rho0*MB,"g/cm^3"
    write(*,"(A10,es18.9,A10,es18.9)") "e_c   :", ee/(C * C * KSCALE),"g/cm^3", ee/(C * C * KSCALE)*rho_uni
    write(*,"(A10,es18.9)")     "rp/re :", r_ratio
    write(*,"(A10,es18.9,A4)")  "OMG_c :", omega_c/2.d0/pi* (C/sqrt(kappa)),"Hz"
    write(*,"(A10,es18.9,A4)")  "Mass  :", Mass/MSUN,"M_o"
    write(*,"(A10,es18.9,A4)")  "M_0   :", Mass_0/MSUN,"M_o"
    write(*,"(A10,es18.9)")     "J     :", ang_mom
    write(*,"(A10,es18.9)")     "chi   :", chi
    write(*,"(A10,es18.9)")     "T/W   :", T_kin/abs(Mass_p - Mass + T_kin)
    write(*,"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
    write(*,"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
    write(*,"(A10,es18.9,A4)")  "r_e/M :", r_circ/1e5 / (Mass/MSUN*1.4769994423016508d0)
    write(*,"(A10,es18.9,A4)")  "Ome_K :", Omega_K
    write(*,"(A10,3es18.9)")    "T_kin :", T_kin/MSUN
    write(*,*) " ===================================="
    write(*,*) " "
  end subroutine check_point

end subroutine MRcurve 

