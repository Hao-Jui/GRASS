module scalar_burning_mod
  use precision_mod, only: wp
  use para_mod, only: active_theory, THEORY_GR, &
                      B_coup, B_burn_init, n_of_relaxation_steps, &
                      mphi_r, mphi_burn_seed, l_uni, KAPPA, &
                      scalar_burn_max_iter, sphi_c, sphi_m, output
  use rotation_solver_mod,  only: rotation_solver
  implicit none
contains
  subroutine perform_scalar_burn(target_mphi)
    real(wp), intent(in) :: target_mphi
    real(wp) :: current_mphi, t0, t1
    integer :: burn_iter

    if (active_theory == THEORY_GR) return

    B_coup = B_burn_init
    mphi_r = (mphi_burn_seed / l_uni)**2 * KAPPA / 1.e10_wp
    print *, " "
    write(*,'("Solving for B =",es15.6," and mphi =",es15.6)') B_burn_init, mphi_burn_seed
    call cpu_time(t0); call rotation_solver(); call cpu_time(t1); write(*,"(A, f10.4)") "Elapsed time [s]:", t1-t0
    print *, " "
    print *, "scalar burn stage:  (B, mphi, sphi_c, sphi_m)"

    current_mphi = sqrt(mphi_r*1.e10_wp/KAPPA) * l_uni

    burn_iter = 0; n_of_relaxation_steps = 0; call cpu_time(t0)
    do while (current_mphi < target_mphi .and. burn_iter < scalar_burn_max_iter)
      call rotation_solver()

      B_coup = merge(B_coup * 1.5_wp, B_coup * 1.2_wp, current_mphi > 5.0_wp)

      if (sphi_m < 0.4_wp) then
        mphi_r = mphi_r * 1.1_wp
      else
        mphi_r = mphi_r * 1.4_wp
      end if

      if (sphi_m < 0.4_wp) B_coup = B_coup * 1.7_wp

      current_mphi = sqrt(mphi_r*1.e10_wp/KAPPA) * l_uni
      if (mod(burn_iter, 10) == 0) then
        write(*,"(i3,A,2(es12.3),A,2(1X,ES18.9))") burn_iter+1, ") ", &
              B_coup, current_mphi, "  |", sphi_c, sphi_m
      end if
      burn_iter = burn_iter + 1
    end do
    output = .true.; call rotation_solver(); output = .false.; call cpu_time(t1)

    write(*,'("Total iterations: ",i0,"  Elapsed time [s]:", f10.4)') n_of_relaxation_steps, t1-t0
    write(*,"(A)") " ", " Solution saved for restart after burning.", " "

    if (burn_iter >= scalar_burn_max_iter .and. current_mphi < target_mphi) then
      write(unit=*, fmt=*) "scalar burn stage reached iteration limit before hitting target mass."
    end if
    write(*,"(A)") " "

    if (sphi_m < 1.e-5_wp) then
      write(*,"(A)") "*** initial guess is non-scalarized"
    else
      write(*,"(A, f10.6, 2es15.6)") "*** Starts with sphi {max|B|mphi} : ", &
              sphi_m, B_coup, current_mphi
    end if

    write(*,"(A)") " "
  end subroutine perform_scalar_burn
end module scalar_burning_mod
