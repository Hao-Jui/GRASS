subroutine shoot_v2
#include "option_macro.h"
  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it, i
  real(8) :: er, rho0, ee

  r_ratio  = 8.418424129E-01
  output   = .false.

#if defined(uryu_law)
  ! with Uryu subroutine, non-rotating NS cannot be constructed
  if ( r_ratio == 1.d0 ) r_ratio = 7.d-1
#endif

#if defined(restart)
  call restart_read
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#elif defined (refine)
  call refine_read
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#else
  e_center = 5.337208333E+14
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
  call sphere
#endif

  write(*,*) " "
  
  it = 1
  er = 1.d99
  do while ( er > accuracy .and. .not.output )
#if defined(uniform)
    call spin
#elif defined (constJ)
    call const_j
#elif defined (fitting)
    call poly
#elif defined (uryu_law)
    call uryu
#endif

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
      write(*,"(A10,es18.9)")     "T/W   :", T_kin/abs(Mass_p - Mass + T_kin)
      write(*,"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
      write(*,"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
      write(*,"(A10,es18.9,A4)")  "r_e/M :", r_circ/1e5 / (Mass/MSUN*1.4769994423016508d0)
      write(*,"(A10,es18.9,A4)")  "Ome_K :", Omega_K
      write(*,"(A10,3es18.9)")    "er    :", er,T_kin/MSUN
      write(*,*) " ===================================="
      write(*,*) " "
    endif
    
    call update_v2(h_center, r_ratio, er)

    it = it + 1
    if (it==1000) stop "Iteration may never converge."
  enddo

  output = .true.
#if defined(uniform)
    call spin
#elif defined (constJ)
    call const_j
#elif defined (fitting)
    call poly
#elif defined (uryu_law)
    call uryu
#endif

  rho0 = n0_at_h(h_center)
  ee   = e_at_h (h_center)
  write(*,*) " "

  !open(74,file="./model_log/model_log.dat",access='append')
  !write(74,"(10es18.9)") ee/(C * C * KSCALE), omega_c/2.d0/pi* (C/sqrt(kappa)), Mass_0/MSUN, &
  !    ang_mom, chi, T_kin/abs(Mass_p - Mass + T_kin)
  !close(74)

  open(221,file="./Cont/properties.dat")
  write(*,*) " ===================================="
  write(*,*) "              Converged              "
  do i = 1, 2
    write(6+215*(i-1),"(A18,es18.9,A8,es18.9)") &
                                          "   Central rho =", rho0*MB,"g/cm^3", rho0*MB*rho_uni
    write(6+215*(i-1),"(A18,es18.9,A10)") "Central energy =", ee / (C * C * KSCALE),"g/cm^3"
    write(6+215*(i-1),"(A18,es18.9)")     "   Axial ratio =", r_ratio
    write(6+215*(i-1),"(A18,f18.9,A4)")   " Central Omega =", omega_c / 2.d0 / pi * (C/sqrt(kappa)),"Hz"
    write(6+215*(i-1),"(A18,f18.9,A4)")   " Equator Omega =", Omega_e / 2.d0 / pi * (C/sqrt(kappa)),"Hz"
    write(6+215*(i-1),"(A18,f18.9,A4)")   "      ADM Mass =", Mass / MSUN,"M_o"
    write(6+215*(i-1),"(A18,f18.9,A16,es18.9,A2)") &
                                          " Baryonic Mass =", Mass_0 / MSUN,"M_o ( binding =", (Mass - Mass_0) / MSUN, " )"
    write(6+215*(i-1),"(A18,f18.9,A10,f7.4,A2)") &
                                          " Ang. Momentum =", ang_mom, " ( chi =", chi, ")"
    write(6+215*(i-1),"(A18,f18.9)")      "        M2/M^3 =", M2
    write(6+215*(i-1),"(A18,f18.9)")      "        S3/M^4 =", S3
    write(6+215*(i-1),"(A18,f18.9)")      "        M4/M^5 =", M4
    write(6+215*(i-1),"(A18,f18.9)")      "           T/W =", T_kin/abs(Mass_p - Mass + T_kin)
    write(6+215*(i-1),"(A18,f18.9,A4)")   "       Coord R =", r_e*sqrt(KAPPA)/1.d5,"km"
    write(6+215*(i-1),"(A18,f18.9,A4)")   "       Areal R =", r_circ/1e5,"km"
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

  deviA = mass/MSUN/M_goal - 1.d0
  !deviB = Omega_c / 2.d0 / pi * (C/sqrt(kappa))/omc_goal - 1.d0
  deviB = chi - chi_goal

  er    = abs(deviA) + abs(deviB)

  if ( abs(deviA) > abs(deviB) ) then
    new = hc - deviA * 1.d-1
    write(*,"(A10,es15.6,A6,es15.6,A8,es15.6)") "hc:", hc, "-->", new, " er:", er
    hc = new
  else
    new = min( 1.d0, rep + deviB*1.d-1 )
    write(*,"(A10,es15.6,A6,es15.6,A8,es15.6)") "rp/re:", rep, "-->", new, " er:", er
    rep = new
  endif

end subroutine update_v2
