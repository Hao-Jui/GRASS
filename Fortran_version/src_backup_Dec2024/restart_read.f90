subroutine restart_read
  use para_mod
  implicit none
  integer :: res_r, res_t
  real(8) :: dum
  integer :: s, m
  
  open(89, file="./Cont/check2D.dat")
  read(89,"(2i5,99es27.17)") res_r, res_t, r_e, e_center, r_ratio
  r_e = r_e / (sqrt(KAPPA)/1.d5)
  if (res_r .ne. SDIV) stop "Difference in the resolution."
  
  do s = 1, SDIV
    do m = 1, MDIV
      read(89,"(99es27.17)") dum,dum,alpha(s,m),gama(s,m),rho(s,m),ww(s,m), & ! 1-6
        dum, dum, enthalpy(s,m), dum, & ! 7-10
        velocity_sq(s,m), omg(s,m) ! 11-12
      ww(s,m) = ww(s,m) / (C/sqrt(kappa))
      omg(s,m) = omg(s,m) / (C/sqrt(kappa))
    enddo
  enddo
  close(89)
  write(*,*) "Restart OK!"


end subroutine restart_read
