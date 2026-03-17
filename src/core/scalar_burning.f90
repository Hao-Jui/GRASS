module scalar_burning_mod
  use precision_mod, only: wp
  use para_mod, only: active_theory, THEORY_GR, &
                      B_coup, B_burn_init, &
                      mphi_r, mphi_burn_seed, l_uni, KAPPA, &
                      scalar_burn_max_iter, sphi_c, sphi_m, output
  use rotation_uniform,  only: rotation_solver
  implicit none
contains
  subroutine perform_scalar_burn(target_mphi)
    real(wp), intent(in) :: target_mphi
    real(wp) :: current_mphi
    integer :: burn_iter
    character(100) :: string

    if (active_theory == THEORY_GR) return

    B_coup = B_burn_init
    mphi_r = (mphi_burn_seed / l_uni)**2 * KAPPA / 1.e10_wp
    print *, " "
    write(*,'("Solving for B =",es15.6," and mphi =",es15.6)') B_burn_init, mphi_burn_seed
    call rotation_solver()
    current_mphi = sqrt(mphi_r*1.e10_wp/KAPPA) * l_uni
    print *, " "
    print *, "scalar burn stage:  (B, mphi, sphi_c, sphi_m)"

    burn_iter = 0
    do while (current_mphi < target_mphi .and. burn_iter < scalar_burn_max_iter)
      call rotation_solver()

      B_coup = merge(B_coup * 1.5e0_wp, B_coup * 1.1e0_wp, current_mphi > 30.e0_wp)

      if (sphi_m < 0.4e0_wp) then
        mphi_r = mphi_r * 1.1e0_wp
      else
        mphi_r = mphi_r * 1.5e0_wp
      end if

      if (sphi_m < 0.5e0_wp) B_coup = B_coup * 1.3e0_wp

      current_mphi = sqrt(mphi_r*1.e10_wp/KAPPA) * l_uni
      if (mod(burn_iter, 10) == 0) then
        write(*,"(i3,A,2(es12.3),A,2(1X,ES18.9))") burn_iter+1, ") ", &
              B_coup, current_mphi, "  |", sphi_c, sphi_m
      end if
      burn_iter = burn_iter + 1
    end do
    output = .true.; call rotation_solver(); output = .false.
    
    write(*,"(A)") " ", " Solution saved for restart after burning.", " "
    write(string,"(f12.4)") sphi_m
    if (burn_iter >= scalar_burn_max_iter .and. current_mphi < target_mphi) then
      write(unit=*, fmt=*) "scalar burn stage reached iteration limit before hitting target mass."
    end if
    write(*,"(A)") " ", merge("*** initial guess is non-scalarized", &
                      "*** Starts with sphi max : "//trim(adjustl(string)), &
                      sphi_m < 1.e-5_wp), " "
  end subroutine perform_scalar_burn
end module scalar_burning_mod
