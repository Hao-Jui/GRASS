subroutine shoot_v2
#include "option_macro.h"
  use para_mod
  use shoot_solver_mod
  use shoot_newton_helpers, only: evaluate_solution, build_jacobian
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it, i, i_ratio
  real(8) :: er, rho0, ee, min_Vrr
  real(8) :: F(2), x(2), delta_x(2), rhs(2), h_new, r_new
  type(newton_state) :: solver_state
  logical :: step_ok
  external :: p_at_e, h_at_p, n0_at_h, e_at_h

  ! > NOTE: with Uryu subroutine, non-rotating NS cannot be constructed
  r_ratio = 6.868048263E-01

#if defined(restart)
  call restart_read
  r_ratio  = 0.5
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#elif defined (refine)
  call refine_read
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
#else
  e_center = .580240E+15
  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
  call sphere
#endif

  write(*,*) " "

  do i_ratio = 1, 40
    it = 1
    er = 1.d99
    call reset_newton_state(solver_state)
    output   = .false.
    do
        call evaluate_solution(h_center, r_ratio, F, rho0, ee, er)
        call to_solver_coords(h_center, r_ratio, x)

        if (solver_state%has_jacobian) then
            call broyden_update(solver_state, x, F)
        else
            call build_jacobian(solver_state, x, F, h_center, r_ratio, rho0, ee, er)
            call to_solver_coords(h_center, r_ratio, x)
        endif

        call commit_state(solver_state, x, F)

        if ( mod(it,10) == 0 ) call print_iter_status(it, rho0, ee, er)
        write(*,"(A10,es15.6)") "er :", er
        if ( er <= accuracy .or. output ) exit

        rhs = -F
        step_ok = solve_linear(solver_state%J, rhs, delta_x)
        if (.not. step_ok) then
            call reset_newton_state(solver_state)
            cycle
        endif

        call clamp_step(delta_x)
        x = x + delta_x
        call from_solver_coords(x, h_new, r_new)
        h_center = h_new
        r_ratio  = r_new

        solver_state%has_jacobian = .true.

        it = it + 1
        if (it==1000) stop "Iteration may never converge."
    end do ! > Newton solve

    !output = .true.; call call_rotation_solver

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)
    write(*,*) " "
    call print_converged_block(rho0, ee)

    !open(74,file="./model_log/model_log.dat",access='append')
    !write(74,"(10es18.9)") ee/(C * C * KSCALE), omega_c/2.d0/pi* (C/sqrt(kappa)), Mass_0/MSUN, &
    !    ang_mom, chi, T_kin/abs(Mass_p - Mass + T_kin)
    !close(74)

    ! > Locate ISCO indices (±)
    min_Vrr = 1.d10
    i_isco_m = 1
    do i = res, res*5/3
        if ( min_Vrr > abs(V_rr_m(i)) ) then
            i_isco_m = i
            min_Vrr = abs(V_rr_m(i))
        endif
    enddo
    
    min_Vrr = 1.d10
    i_isco_p = 1
    do i = res, res*5/3
        if ( min_Vrr > abs(V_rr_p(i)) ) then
            i_isco_p = i
            min_Vrr = abs(V_rr_p(i))
        endif
    enddo

    open(771,file="./Cont/Kep_"//trim(adjustl(eos_file))//".log",POSITION='APPEND')
    write(771,"(99es18.9)") omega_c / 2.d0 / pi * (C/sqrt(kappa)), &
                            chi, &
                            s_gp(i_isco_m) / (1.d0-s_gp(i_isco_m)), &
                            s_gp(i_isco_p) / (1.d0-s_gp(i_isco_p)), &
                            r_e*sqrt(KAPPA) / 1.d5, &
                            Mass / MSUN, &
                            (C/sqrt(kappa)) * v_minus(i_isco_m) / r_e, &
                            (C/sqrt(kappa)) * v_plus(i_isco_p) / r_e, &
                            !omega_c / 2.d0 / pi * (C/sqrt(kappa)) * r_e*sqrt(KAPPA) / 1.d5, &
                            ( Omega_e * (C/sqrt(kappa)) ) / Omega_K, Mb_goal
    close(771)
    Mb_goal = Mb_goal + 0.05d0
    !r_ratio = merge( 1.d0 - 0.48d0 * sin( pi / 2.d0 * dble(i_ratio) / 50.d0 )**3, r_ratio - 2.d-4*i_ratio, i_ratio > 5)
  enddo

end subroutine shoot_v2

subroutine print_iter_status(it, rho0, ee, er)
  use para_mod
  implicit none
  integer, intent(in) :: it
  real(8), intent(in) :: rho0, ee, er

  write(*,*) " ===================================="
  write(*,"(A10,I6)")           " iter :", it/10
  write(*,"(A10,ES18.9,A10)")    "rho_c :", rho0*MB,"g/cm^3"
  write(*,"(A10,ES18.9,A10,ES18.9)") "e_c   :", ee/(C*C*KSCALE),"g/cm^3", ee/(C*C*KSCALE)*rho_uni
  write(*,"(A10,ES18.9)")       "rp/re :", r_ratio
  write(*,"(A10,ES18.9,A4)")    "OMG_c :", omega_c/(2.d0*pi)*(C/sqrt(kappa)),"Hz"
  write(*,"(A10,ES18.9,A4)")    "Mass  :", Mass/MSUN,"M_o"
  write(*,"(A10,ES18.9,A4)")    "M_0   :", Mass_0/MSUN,"M_o"
  write(*,"(A10,ES18.9)")       "J     :", ang_mom
  write(*,"(A10,ES18.9)")       "chi   :", chi
  write(*,"(A10,ES18.9)")       "T/W   :", T_kin/abs(Mass_p - Mass + T_kin)
  write(*,"(A10,ES18.9,A4)")    "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
  write(*,"(A10,ES18.9,A4)")    "R_cir :", r_circ/1.d5,"km"
  write(*,"(A10,ES18.9,A4)")    "r_e/M :", (r_circ/1.d5)/(Mass/MSUN*1.4769994423016508d0)
  write(*,"(A10,2ES18.9)")      "Om_K/e:", Omega_K/(2.d0*pi), Omega_e/(2.d0*pi)*(C/sqrt(kappa))
  write(*,"(A10,3ES18.9)")      "er    :", er
  write(*,*) " ===================================="
  write(*,*) " "
end subroutine print_iter_status

subroutine print_converged_block(rho0, ee)
#include "option_macro.h"
  use para_mod
  implicit none
  real(8), intent(in) :: rho0, ee
  integer :: i

  open(221, file="./Cont/properties.dat")
  write(*,*) " ===================================="
  write(*,*) "              Converged              "
  do i = 1, 2
     write(6+215*(i-1),"(A18,ES18.9,A8,ES18.9)")  &
          "   Central rho =", rho0*MB,"g/cm^3", rho0*MB*rho_uni
     write(6+215*(i-1),"(A18,ES18.9,A10)")        &
          "Central energy =", ee/(C*C*KSCALE), "g/cm^3"
     write(6+215*(i-1),"(A18,ES18.9)")            "   Axial ratio =", r_ratio
     write(6+215*(i-1),"(A18,F18.9,A4)")          " Central Omega =", Omega_c/(2.d0*pi)*(C/sqrt(kappa)), "Hz"
     write(6+215*(i-1),"(A18,F18.9,A4)")          " Equator Omega =", Omega_e/(2.d0*pi)*(C/sqrt(kappa)), "Hz"
     write(6+215*(i-1),"(A18,F18.9,A4)")          "   Kepler Omega =", Omega_K/(2.d0*pi), "Hz"
     write(6+215*(i-1),"(A18,F18.9,A4)")          "      ADM Mass =", Mass/MSUN, "M_o"
     write(6+215*(i-1),"(A18,F18.9,A16,ES18.9,A2)") &
          " Baryonic Mass =", Mass_0/MSUN, "M_o ( binding =", (Mass - Mass_0)/MSUN, " )"
     write(6+215*(i-1),"(A18,F18.9,A10,F7.4,A2)") &
          " Ang. Momentum =", ang_mom, " ( chi =", chi, ")"
#ifndef Fishbone
     write(6+215*(i-1),"(A18,F18.9)")             "        M2/M^3 =", M2
     write(6+215*(i-1),"(A18,F18.9)")             "        S3/M^4 =", S3
     write(6+215*(i-1),"(A18,F18.9)")             "        M4/M^5 =", M4
#endif
     write(6+215*(i-1),"(A18,F18.9)")             "           T/W =", T_kin/abs(Mass_p - Mass + T_kin)
     write(6+215*(i-1),"(A18,F18.9,A4)")          "       Coord R =", r_e*sqrt(KAPPA)/1.d5,"km"
     write(6+215*(i-1),"(A18,F18.9,A4)")          "       Areal R =", r_circ/1.d5,"km"
  end do

  write(*,*) " "
  write(*,*) "In code unit:"
  write(*,"(A6,ES18.9,2X,A4,ES18.9,2X,A5,ES18.9,2X,A8,ES18.9)")  &
       "h_c", h_center, "r_e", r_e, "Fmax", Fmax_h, "Omega_e", Omega_e*r_e
  write(*,*) " ===================================="
  close(221)
  write(*,*) " "
  write(*,*) "Completed!"
end subroutine print_converged_block

