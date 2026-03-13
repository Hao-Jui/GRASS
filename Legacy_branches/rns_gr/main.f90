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
  write(*,"(A13,i5,A4,i5)") "Resolution:", SDIV, "x", MDIV

  call shoot_v2
  !call shoot_unstable
  !call J_seq

end program rns


