program grass
  use toolkit_mod, only: debug_mod_bessel
  use para_mod
  implicit none
  character(len=64) :: theory_arg

  call initialize_theory()

  print *, " "
  print *, " "
  if (active_theory == THEORY_ST) then
    write(*,"(A32,es12.6,A11,es12.6,A2)") &
    "Theory: Scalar-Tensor  ( mphi = ", mphi_goal*scalarton, "eV,  B = ", B_goal, " )"
  else
    write(*,*) "Theory: General Relativity"
  end if
  call loadEos

  write(*,'(1X,A)') repeat('-', 36)
  select case (trim(adjustl(solver_type)))
  case ("const_j")
    write(*,'(1X,A,F8.3,A)') "Rotation Law: Constant-J (A^-1 = ", 1.d0 / A_diff, ")"
  case default
    write(*,'(1X,A)') "Rotation Law: Uniform"
  end select
  write(*,'(1X,A)') repeat('-', 36)

  call make_grid
  call GridTrig

  write(*,fmt="(A18,es18.9)",Advance='NO') " Outer boundary:", (s_gp(SDIV-1) / ( 1.d0 - s_gp(SDIV-1) ))**s_pwr
  write(*,"(A18,i5,A4,i5)") "Resolution:", SDIV, "x", MDIV
  write(*,"(A18,es18.9)") "Surface eps:", e_surface / (C*C*KSCALE)
  
  !call debug_mod_bessel; stop "debug_mod_bessel output written"
  
  !call MRcurve
  call shoot_v2

end program grass
