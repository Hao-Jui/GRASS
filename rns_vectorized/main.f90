program rns
#include "option_macro.h"
  use para_mod
  implicit none
  
  call loadEos
  print *, " "
#if defined(uniform)
  write(*,*) "Rotation law: Uniform"
#elif defined (constJ)
  write(*,"(A32,2f9.3)") "Rotation law: Const J with A^-1 =", 1.d0/A_diff
#elif defined (fitting)
  write(*,"(A32,2f9.4)") "Rotation law: Poly with parA and B =", parA, parB
#elif defined (uryu_law)
  write(*,"(A32,2f9.3)") "Rotation law: Uryu with lambda =", lambda1, lambda2
#endif
  call make_grid
  call GridTrig

  !j_disk = j_disk * (G*Msun/C) / (sqrt(KAPPA)*C)
  write(*,fmt="(A18,es18.9)",Advance='NO') " Outer boundary:", (s_gp(SDIV-1) / ( 1.d0 - s_gp(SDIV-1) ))**s_pwr
  write(*,"(A18,i5,A4,i5)") "Resolution:", SDIV, "x", MDIV
  write(*,"(A18,es18.9)") "Surface eps:", e_surface / (C*C*KSCALE)
#if defined(Fishbone)
  write(*,"(A18,es18.9)") " Fishbone j:", j_disk
#endif
  call shoot_v2
  !call shoot_unstable
  !call J_seq

end program rns

subroutine call_rotation_solver
  use para_mod
  implicit none

  select case (trim(adjustl(solver_type)))
  case ("uniform")
     call spin
  case ("constJ")
     call const_j
  case ("uryu_law")
     call uryu
  case default
     write(*,*) "ERROR: Unknown solver type:", solver_type
     stop
  end select
end subroutine call_rotation_solver


