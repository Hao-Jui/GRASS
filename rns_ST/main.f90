program rns
#include "option_macro.h"
  use para_mod
  implicit none
  

  call loadEos
  print *, " "
  
  call make_grid
  call GridTrig
  write(*,"(A13,i5,A4,i5)") "Resolution:", SDIV, "x", MDIV
  write(*,"(A13,es15.6)") "Scalar mass:",mphi_
  write(*,"(A13,es15.6)") "Coupling B:",B_

  call shoot_v2
  !call shoot_unstable
  !call J_seq

end program rns


