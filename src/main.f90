program grass
  use eos_mod, only: loadEos
  use toolkit_mod, only: debug_mod_bessel
  use para_mod, only: wp, initialize_theory, &
                      has_scalar, mphi_goal, scalarton, B_goal, &
                      solver_type, A_diff, lambda1, lambda2, &
                      s_gp, SDIV, MDIV, s_pwr, e_surface, C, KSCALE
  use constrain_mod, only: hamiltonian
  implicit none
  real(wp) :: hamL2

  call initialize_theory()
  
  print *, " "
  print *, " "
  if (has_scalar) then
    write(*,"(A32,es12.6,A11,es12.6,A2)") &
    "Theory: Scalar-Tensor  ( mphi = ", mphi_goal*scalarton, "eV,  B = ", B_goal, " )"
  else
    write(*,*) "Theory: General Relativity"
  end if
  write(*,*) "Relaxation scheme: Anderson"

  call loadEos

  write(*,'(1X,A)') repeat('-', 36)
  select case (trim(adjustl(solver_type)))
  case ("const_j")
    write(*,'(1X,A,F5.3,A)') "Rotation Law: Constant-J   (A^-1 = ", 1.d0 / A_diff, ")"
  case ("uryu")
    write(*,'(1X,A,F5.3,A,F5.3,A)') "Rotation Law: Uryu   (lambda1 = ", lambda1, "  lambda2 = ", lambda2, ")"
  case default
    write(*,'(1X,A)') "Rotation Law: Uniform"
  end select
  write(*,'(1X,A)') repeat('-', 36)

  call make_grid
  call GridTrig

  write(*,fmt="(A18,es18.9)",Advance='NO') " Outer boundary:", (s_gp(SDIV-1) / ( 1.d0 - s_gp(SDIV-1) ))**s_pwr
  write(*,"(A18,i0,A,i0)") "Resolution: ", SDIV, " x ", MDIV
  write(*,"(A18,es18.9)") "Surface eps:", e_surface / (C*C*KSCALE)
  
  !call debug_mod_bessel; stop "debug_mod_bessel output written"
  
  !call MRcurve
  call shoot_v2

  call hamiltonian(hamL2)
  write(*,"(A18,es12.4)") "Ham L2:", hamL2
  write(*,*) " "
  write(*,*) "Completed!"

end program grass
