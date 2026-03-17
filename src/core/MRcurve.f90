subroutine MRcurve
  use para_mod, only: wp, h_center, mass, mass_0, mass_p, r_ratio, r_circ, &
                      sphi, sphi_c, B_coup, Omega_e, KAPPA, &
                      MB, MSUN, C, KSCALE, rho_uni, n_sat, &
                      I_inertia, Love2, M2, M4, S3, chi, T_kin, &
                      mphi_goal, B_goal, eos_file, sound_speed
  use rotation_uniform,  only: rotation_solver
  use starting_model_mod, only: initialize_starting_model
  implicit none
  real(wp), external :: p_at_e, h_at_p, n0_at_h, e_at_h
  integer :: it, unit
  real(wp) :: rho0, ee, pp, sphi_max, mass_prev
  real(wp) :: Q_bar, T_over_W, traceT
  character(16) :: fil1, fil2
  logical :: ascending = .false.

  call initialize_starting_model(p_at_e, h_at_p, n0_at_h, e_at_h)
  rho0 = n0_at_h(h_center); mass_prev = mass
  r_ratio  = 1.e0_wp

  it = 0
  write(*,"(A5,6A15)") "Iter", "mass density", "ADM mass", "MoI", "Love", "Q bar", "varphi"
  do while (merge(mass_prev <= mass .and. mass > 2.e0_wp, rho0*MB > 5.e13_wp, ascending))
    mass_prev = mass
    if (ascending) then
      h_center = h_center * 1.01e0_wp
    else
      h_center = h_center * 1.01e0_wp
    end if
    call rotation_solver()
    call solution_properties()

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)
    pp   = p_at_e (ee)
    sphi_max = maxval(sphi(:,1)) * sqrt(B_coup) 

    call output_seq()
    
    if ( mod(it,10) == 0 ) write(*,"(i5,6es15.6)") it, rho0*MB, Mass/MSUN, I_inertia, Love2, Q_bar, sphi_max
    it = it + 1
  enddo 

contains
  subroutine output_seq()
    traceT = 3.e0_wp*pp*1.80171810e-39_wp/KSCALE - ee*rho_uni/(C * C * KSCALE)
    Q_bar = merge(-1.d0, M2/chi**2, chi==0.d0)
    T_over_W=merge(-1.d0, T_kin/abs(Mass_p - Mass + T_kin), chi==0.d0)
    write(fil1,"(f10.1)") mphi_goal
    write(fil2,"(es12.1)") B_goal
    if (r_ratio /= 1.e0_wp) then
      open(newunit=unit,file="/Users/horay/Data4Projects/crazy/Seq/"//trim(adjustl(eos_file))// &
            "/mphi"//trim(adjustl(fil1))//"_B"//trim(adjustl(fil2))//"_spin.dat",access='append')
    else
      open(newunit=unit,file="/Users/horay/Data4Projects/crazy/Seq/"//trim(adjustl(eos_file))// &
            "/mphi"//trim(adjustl(fil1))//"_B"//trim(adjustl(fil2))//".dat",access='append')
    end if
    write(unit,"(99es18.9e3)") ee/(C * C * KSCALE), rho0*MB/n_sat, h_center, sound_speed(1), & ! 1-4
            Mass/MSUN, Mass_0/MSUN, r_circ/1e5, & ! 5-7
            sphi_c, sphi_max, & ! 8, 9
            I_inertia, Love2, Q_bar, & ! 10-12
            T_over_W, M2, S3, M4, chi, Omega_e* (C/sqrt(kappa)), traceT ! 13-19
    close(unit)
  end subroutine output_seq
end subroutine MRcurve 

