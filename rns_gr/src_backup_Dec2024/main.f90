program rns
#include "option_macro.h"
  use para_mod
  implicit none
  
  eos_file = "./eos/MPA1.dat"
  call loadEos
  print *, " "
  
  call make_grid
  call GridTrig
  write(*,"(A13,i5,A4,i5)") "Resolution:", SDIV, "x", MDIV

  call shoot_v2
  !call shoot_unstable
  !call J_seq

end program rns
