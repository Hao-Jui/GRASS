subroutine shoot_v2
#include "option_macro.h"
  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it
  real(8) :: er, rho0, ee, mphi0, mphi1
  real(8) :: B_ini = 16.d0

  output   = .false.

#if defined(restart)
  call restart_read
  !e_center = 1.75d15
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#elif defined(refine)
  call refine_read
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#else
  r_ratio  = 1.d0
  e_center = 1.89d15
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
  call sphere
  
  mphi_r  = (0.01d0/l_uni)**2 * KAPPA / 1.d10
  mphi_ep = mphi_r
  B_coup = B_ini
  call spin

  !do while ( sqrt(mphi_r*1.d10/KAPPA)*l_uni < mphi_ )
  !  mphi_r  = mphi_r * 1.5d0; mphi_ep = mphi_r
  !  call spin
  !  if ( sphi_c < 0.0 1d0 ) then 
  !      B_coup = B_coup * 1.5d0
  !  else
  !      B_coup = B_coup * 1.2d0
  !  endif
  !  if ( sphi_c < 0.001d0 ) then
  !      mphi_r = mphi_r / 1.2d0
  !  elseif( sphi_c > 0.5d0 ) then
  !      mphi_r = mphi_r * 4.d0
  !  endif
  !  write(*,"(4es18.9)") B_coup, sqrt(mphi_r*1.d10/KAPPA)*l_uni, sphi_c
  !enddo

  write(*,*) " "
  write(*,*) "Starts with sphi_c:", sphi_c
#endif
  it = 1
  er = 1.d99


  mphi_r  = (mphi_/l_uni)**2 * KAPPA / 1.d10
  B_coup  = B_
  rho0    = n0_at_h(h_center)
  
!  do while ( rho0*MB > 1.d14 )
  mphi_r  = mphi_r !* 1.2d0
  mphi_ep = mphi_r !(mphi_/l_uni)**2 * (C**4/G/1.d12) * KSCALE

  do while ( er > accuracy .and. .not.output )

    if (mphi_==0.d0) then 
      call spin_massless
    else
      call spin
    endif
    call mass_radius

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)

    if ( mod(it,10) == 0 ) then
    write(*,*) " ===================================="
    write(*,"(A10,i6)")         " iter :", it/10
    write(*,"(A10,es18.9)")     " sphi :", sphi_c
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
    if ( it == 1000 ) stop "Iteration may never converge."
  enddo
  
  output = .true.

  if ( mphi_ == 0.d0 ) then 
    call spin_massless
  else
    call spin
  endif

#if 0
  if ( sphi_c < 0.5d0 ) then 
      B_coup = B_coup * 1.5d0
  else
      B_coup = B_coup * 1.2d0
  endif
  if ( sphi_c < 0.1d0 ) then
      mphi_r = mphi_r / 1.2d0
  elseif( sphi_c > .8d0 ) then
      mphi_r = mphi_r * 1.5d0
  endif
#endif

  call mass_radius
  rho0 = n0_at_h(h_center)
  ee   = e_at_h (h_center)
  write(*,*) " "

  open(221,file="./Cont/properties_unstable.dat")
  write(*,*) " ============================================"
  write(*,*) "                  Converged                  "
  write(*,*) " "
  do it = 1, 2
    write(6+215*(it-1),"(A18,f18.9)")      "    Coupling B =", B_coup
    write(6+215*(it-1),"(A18,f18.9)")      "   Scalar mass =", sqrt(mphi_r*1.d10/KAPPA)*l_uni
    write(6+215*(it-1),"(A18,es18.9,A4,es18.9)") "   Central rho =", rho0*MB,"g/cm^3", rho0*MB*rho_uni
    write(6+215*(it-1),"(A18,es18.9,A10)") "Central energy =", ee/(C * C * KSCALE),"g/cm^3"
    write(6+215*(it-1),"(A18,f18.9)")      "     varphi(0) =", sphi_c 
    write(6+215*(it-1),"(A18,es18.9)")     "   Axial ratio =", r_ratio
    write(6+215*(it-1),"(A18,f18.9,A4)")   " Central Omega =", omega_c/2.d0/pi* (C/sqrt(kappa)),"Hz"
    write(6+215*(it-1),"(A18,f18.9,A4)")   " Equator Omega =", Omega_e,"Hz"
    write(6+215*(it-1),"(A18,f18.9,A4)")   "      ADM Mass =", Mass/MSUN,"M_o"
    write(6+215*(it-1),"(A18,f18.9,A4)")   " Baryonic Mass =", Mass_0/MSUN,"M_o"
    write(6+215*(it-1),"(A18,f18.9)")      " Ang. Momentum =", ang_mom
    write(6+215*(it-1),"(A18,f18.9)")      "           chi =", chi
    write(6+215*(it-1),"(A18,f18.9)")      "           T/W =", T_kin/abs(Mass_p - Mass + T_kin)
    write(6+215*(it-1),"(A18,f18.9,A4)")   "       Coord R =", r_e*sqrt(KAPPA)/1.d5,"km"
    write(6+215*(it-1),"(A18,f18.9,A4)")   "       Areal R =", r_circ/1e5,"km"
  enddo
  write(*,*) " "
  write(*,*) "In code unit:"
  write(*,"(A10,es18.9)") "h_c", h_center, "r_e", r_e, "Fmax", Fmax_h, "Omega_e", Omega_e * r_e
  write(*,*) " ============================================"
  close(221)

      open(74,file="/Users/horay/Data4Projects/crazy/Seq/"//trim(adjustl(eos_file))// &
              "/1d_mphi10_B30000.dat",access='append')
      write(74,"(10es18.9)") ee/(C * C * KSCALE), rho0*MB, h_center, Mass/MSUN, r_circ/1e5, &
          Mass_0/MSUN, sphi_c, T_kin/abs(Mass_p - Mass + T_kin)
      close(74)

!  h_center = h_center - 1.d-3
!  enddo 
  
  write(*,*) " "
  write(*,*) "Completed!"

end subroutine shoot_v2


subroutine update_v2(hc, rep, er)
  use para_mod
  implicit none
  real(8), intent(inout) :: hc, rep
  real(8), intent(out) :: er
  real(8) :: new, deviA, deviB

  deviA = mass/MSUN - M_goal
  deviB = ang_mom - J_goal


  er = abs(deviA) + abs(deviB)

  if (abs(deviA) > abs(deviB) ) then
    new = hc - deviA * 1.d-1
    write(*,"(A10,es15.6,A10,es15.6,A10,es15.6)") "hc:",hc,"-->",new," er:",er
    hc = new
  else
    new = min( 1.d0, rep + sign( sqrt(abs(deviB)), deviB)*3.d-3 )
    write(*,"(A10,es15.6,A10,es15.6,A10,es15.6)") "| rp/re:",rep,"-->", new," er:",er
    rep = new
  endif

end subroutine update_v2
