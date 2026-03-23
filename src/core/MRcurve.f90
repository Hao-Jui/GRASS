module MRcurve_mod
  implicit none
contains

subroutine MRcurve
  use analysis_mod, only: solution_properties
  use eos_mod, only: p_at_e, n0_at_h, e_at_h
  use para_mod, only: wp, h_center, mass, mass_0, mass_p, r_ratio, r_circ, &
                      sphi, sphi_c, B_coup, Omega_e, KAPPA, &
                      MB, MSUN, C, KSCALE, rho_uni, n_sat, &
                      I_inertia, Love2, M2, M4, S3, chi, T_kin, &
                      mphi_goal, B_goal, eos_file, sound_speed
  use rotation_uniform,  only: rotation_solver
  use starting_model_mod, only: initialize_starting_model
  implicit none
  integer :: it, unit
  real(wp) :: rho0, ee, pp, sphi_max, mass_prev
  real(wp) :: Q_bar, T_over_W, traceT
  character(16) :: fil1, fil2
  logical :: ascending = .false.

  call initialize_starting_model()
  rho0 = n0_at_h(h_center); mass_prev = mass
  r_ratio  = 1.e0_wp

  it = 0
  write(*,"(A5,6A15)") "Iter", "mass density", "ADM mass", "MoI", "Love", "Q bar", "varphi"
  do while (rho0*MB > 5.e13_wp .and. mass/MSUN > 0.1_wp)
    mass_prev = mass
    if (ascending) then
      h_center = h_center * 1.005_wp
    else
      h_center = h_center / 1.005_wp
    end if
    call rotation_solver()
    call solution_properties()

    rho0 = n0_at_h(h_center)
    ee   = e_at_h (h_center)
    pp   = p_at_e (ee)
    sphi_max = maxval(sphi(:,1)) * sqrt(B_coup) 

    call output_seq()
    
    if ( mod(it,10) == 0 ) write(*,"(i5,6es15.6)") it, rho0*MB, Mass/MSUN, I_inertia, Love2, Q_bar, sphi_max
    if (mass_prev > mass .and. mass/MSUN > 2.0_wp) exit
    it = it + 1
  enddo 

contains
  subroutine output_seq()
    character(len=1024) :: filename
    character(len=32)   :: suffix
    integer :: ios
    traceT  = 3.0_wp*pp*1.80171810e-39_wp/KSCALE - ee*rho_uni/(C * C * KSCALE)
    Q_bar   = merge(-1.0_wp, M2/chi**2, abs(chi) < epsilon(chi))
    T_over_W= merge(-1.0_wp, T_kin/abs(Mass_p - Mass + T_kin), abs(chi) < epsilon(chi))

    write(fil1,"(f10.1)") mphi_goal
    write(fil2,"(es12.1)") B_goal

    if (abs(r_ratio - 1.0_wp) >= epsilon(r_ratio)) then
        suffix = "_spin.dat"
    else
        suffix = ".dat"
    end if

    write(filename, '(A, A, A, A, A, A, A)') &
      "/Users/horay/Data4Projects/crazy/PT/", trim(eos_file), &
      "/mphi", trim(adjustl(fil1)), "_B", trim(adjustl(fil2)), trim(suffix)

    open(newunit=unit, file=trim(filename), access='append', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,I0,A)') 'ERROR: Cannot open output file. IOSTAT = ', ios, trim(filename)
    end if
    write(unit,"(99es18.9e3)") &
            ee/(C * C * KSCALE), rho0*MB/n_sat, & ! 1-2
            h_center, sound_speed(1), &           ! 3-4
            Mass/MSUN, Mass_0/MSUN, r_circ/1e5_wp, & ! 5-7
            sphi_c, sphi_max, &                   ! 8, 9
            I_inertia, Love2, Q_bar, T_over_W, &  ! 10-13
            M2, S3, M4, chi, Omega_e*(C/sqrt(kappa)), traceT ! 14-19
    close(unit)
  end subroutine output_seq
end subroutine MRcurve 

end module MRcurve_mod

