subroutine restart_read
  use para_mod
  implicit none
  integer :: res_r, res_t
  real(8) :: dum
  integer :: s, m
  
  open(89, file="./Res/res.dat")
  read(89,"(2i5,99es27.17)") res_r, res_t, r_e, e_center, r_ratio
  r_e = r_e / (sqrt(KAPPA)/1.d5)
  if (res_r .ne. SDIV) stop "Difference in the resolution."
  
  do s = 1, SDIV
    do m = 1, MDIV
      read(89,"(99es27.17)") dum,dum,alpha(s,m),gama(s,m),rho(s,m),ww(s,m), & ! 1-6
        dum, dum, enthalpy(s,m), dum, & ! 7-10
        velocity_sq(s,m), omg(s,m), sphi(s,m) ! 11-12
      ww(s,m) = ww(s,m) / (C/sqrt(kappa))
      omg(s,m) = omg(s,m) / (C/sqrt(kappa))
      sphi(s,m)= sphi(s,m) / sqrt(B_)
    enddo
  enddo
  close(89)
  write(*,*) " "
  write(*,*) "Restart OK!"


end subroutine restart_read

subroutine refine_read
  use para_mod
  implicit none
  integer :: res_r, res_t
  real(8) :: dum
  integer :: s, m, k
  
  open(89, file="./Res/res.dat" )
  read(89,"(2i5,99es27.17)") res_r, res_t, r_e, e_center, r_ratio, Omega_e, Omega_c
  r_e = r_e / (sqrt(KAPPA)/1.d5)
  Omega_e = Omega_e / (C/sqrt(kappa)) * r_e
  Omega_c = Omega_c / (C/sqrt(kappa)) * r_e
  if ( res_r*2 .ne. SDIV+1 ) then 
    write(*,*) res_r*2, SDIV+1
    stop "Difference in the resolution."
  endif
  
  do s = 1, SDIV/2+1
    do m = 1, MDIV/2+1
      if ( s.ne.SDIV/2+1 .and. m.ne.MDIV/2+1 ) then
        read(89,"(99es27.17)") dum,dum,alpha(s*2,m*2),gama(s*2,m*2), & ! 1-4
          rho(s*2,m*2),ww(s*2,m*2), & ! 5-6
          dum, dum, enthalpy(s*2,m*2), dum, & ! 7-10
          velocity_sq(s*2,m*2), omg(s*2,m*2), sphi(s*2,m*2) ! 11-12

          do k = 1, 3
            alpha(s*2-mod(k,2),m*2-k/2) = alpha(s*2,m*2)
            gama (s*2-mod(k,2),m*2-k/2) = gama (s*2,m*2)
            rho  (s*2-mod(k,2),m*2-k/2) = rho  (s*2,m*2)
            ww   (s*2-mod(k,2),m*2-k/2) = ww   (s*2,m*2)
            omg  (s*2-mod(k,2),m*2-k/2) = omg  (s*2,m*2)
            sphi (s*2-mod(k,2),m*2-k/2) = sphi (s*2,m*2)
            enthalpy(s*2-1,m*2-1) = enthalpy(s*2,m*2)
          enddo
      elseif ( s == SDIV/2+1 .and. m.ne.MDIV/2+1 ) then
        read(89,"(99es27.17)") dum,dum,alpha(s*2-1,m*2),gama(s*2-1,m*2), & ! 1-4
          rho(s*2-1,m*2),ww(s*2-1,m*2), & ! 5-6
          dum, dum, enthalpy(s*2-1,m*2), dum, & ! 7-10
          velocity_sq(s*2-1,m*2), omg(s*2-1,m*2), sphi(s*2-1,m*2) ! 11-12

          alpha(s*2-1,m*2-1) = alpha(s*2-1,m*2)
          gama (s*2-1,m*2-1) = gama (s*2-1,m*2)
          rho  (s*2-1,m*2-1) = rho  (s*2-1,m*2)
          ww   (s*2-1,m*2-1) = ww   (s*2-1,m*2)
          omg  (s*2-1,m*2-1) = omg  (s*2-1,m*2)
          sphi (s*2-1,m*2-1) = sphi (s*2-1,m*2)
      elseif ( m == MDIV/2+1 .and. s.ne.SDIV/2+1 ) then
        read(89,"(99es27.17)") dum,dum,alpha(s*2,m*2-1),gama(s*2,m*2-1), & ! 1-4
          rho(s*2,m*2-1),ww(s*2,m*2-1), & ! 5-6
          dum, dum, enthalpy(s*2,m*2-1), dum, & ! 7-10
          velocity_sq(s*2,m*2-1), omg(s*2,m*2-1), sphi(s*2,m*2-1) ! 11-12

          alpha(s*2-1,m*2-1) = alpha(s*2,m*2-1)
          gama (s*2-1,m*2-1) = gama (s*2,m*2-1)
          rho  (s*2-1,m*2-1) = rho  (s*2,m*2-1)
          ww   (s*2-1,m*2-1) = ww   (s*2,m*2-1)
          omg  (s*2-1,m*2-1) = omg  (s*2,m*2-1)
          sphi (s*2-1,m*2-1) = sphi (s*2,m*2-1)
      else
        read(89,"(99es27.17)") dum,dum,alpha(s*2-1,m*2-1),gama(s*2-1,m*2-1), & ! 1-4
          rho(s*2-1,m*2-1),ww(s*2-1,m*2-1), & ! 5-6
          dum, dum, enthalpy(s*2-1,m*2-1), dum, & ! 7-10
          velocity_sq(s*2-1,m*2-1), omg(s*2-1,m*2-1), sphi(s*2-1,m*2-1) ! 11-12
      endif
    enddo
  enddo
  close(89)
  ww  = ww  / (C/sqrt(kappa))
  omg = omg / (C/sqrt(kappa))
  sphi= sphi / sqrt(B_)
  write(*,*) " "
  write(*,*) "Refine-read OK!"

end subroutine refine_read
