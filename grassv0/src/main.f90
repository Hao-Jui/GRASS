program rns
#include "option_macro.h"
  use toolkit_mod, only: debug_mod_bessel
  use para_mod
  implicit none
  character(len=64) :: theory_arg
  integer :: parse_status

  call get_command_argument(1, theory_arg)
  call initialize_theory_from_string(trim(theory_arg), parse_status)
  if (parse_status /= 0) then
    write(*,*) "Unknown theory option: ", trim(theory_arg)
    write(*,*) "Use 'gr' or 'st'."
    stop 1
  endif
  print *, " "
  if (has_scalar) then
    write(*,"(A32,es12.6,A11,f10.2,A2)") &
    "Theory: Scalar-Tensor  ( mphi = ", mphi_goal*1.33d-10, "eV,  B = ", B_goal, " )"
  else
    write(*,*) "Theory: General Relativity"
  end if
  call loadEos
  print *, " "

  write(*,'(1X,A)') repeat('-', 36)
  write(*,'(1X,A)') "Rotation Law: Uniform"
  write(*,'(1X,A)') repeat('-', 36)

  call make_grid
  call GridTrig

  write(*,fmt="(A18,es18.9)",Advance='NO') " Outer boundary:", (s_gp(SDIV-1) / ( 1.d0 - s_gp(SDIV-1) ))**s_pwr
  write(*,"(A18,i5,A4,i5)") "Resolution:", SDIV, "x", MDIV
  write(*,"(A18,es18.9)") "Surface eps:", e_surface / (C*C*KSCALE)
  
  !call debug_mod_bessel; stop "debug_mod_bessel output written"
  
  call shoot_v2

end program rns
