!-*-f90-*-
module parameters

  real*8 Gamma1,Gamma2,Gamma3,log_p1
  real*8 eps,prs,hh,cspead,rho,log_rho,dd(6)
  real*8 rho_match,pp,K1,K2,K3,step
  real*8 K_practical,G_practical,err,aa,aaa(7)

  character(len=128):: eos_name = "APR1"

  real*8, parameter :: clite = 2.99792458d10
  real*8, parameter :: dens_division1 = 14.7
  real*8, parameter :: dens_division2 = 15.0

  real*8, parameter :: crust_Gamma1 = 1.58425
  real*8, parameter :: crust_Gamma2 = 1.28733
  real*8, parameter :: crust_Gamma3 = 0.62223
  real*8, parameter :: crust_Gamma4 = 1.35692
  real*8, parameter :: crust_K1 = 6.80110d-9
  real*8, parameter :: crust_K2 = 1.06186d-6
  real*8, parameter :: crust_K3 = 5.32697d1
  real*8, parameter :: crust_K4 = 3.99874d-8
  real*8, parameter :: crust_devision1 = log10(2.44034d7)
  real*8, parameter :: crust_devision2 = log10(3.78358d11)
  real*8, parameter :: crust_devision3 = log10(2.62780d12)
  
  real(8), parameter:: MB = 1.6749286d-24 ! neutron mass
  real(8), parameter:: rho_uni = 1.d0/1.61930347d-18


end module parameters
program piecewisepoly

  use parameters
  implicit none
  integer ll
  

  call setup_table

  ! matching low density EOS with the high density one
  pp = log10(crust_K4) + crust_Gamma4*crust_devision3
  rho_match = (log_p1-2*log10(clite)-pp-Gamma1*dens_division1+crust_Gamma4*crust_devision3)/(crust_Gamma4-Gamma1)
  write(*,"(A12,2es18.9)") "rho_match=",10**rho_match, 10**rho_match/rho_uni
  ! fix K's
  K1 = crust_K4*(1.d1**rho_match)**(crust_Gamma4-Gamma1)
  K2 = K1*(1.d1**dens_division1)**(Gamma1-Gamma2)
  K3 = K2*(1.d1**dens_division2)**(Gamma2-Gamma3)

  dd(1) = 2.44034d7
  dd(2) = 3.78358d11
  dd(3) = 2.62780d12
  dd(4) = 1.d1**rho_match
  dd(5) = 1.d1**dens_division1
  dd(6) = 1.d1**dens_division2

  aaa(1) = 0
  aaa(2) = ((1+aaa(1))*dd(1) + crust_K1*dd(1)**crust_Gamma1/(crust_Gamma1-1))/dd(1) - 1 &
    - crust_K2*dd(1)**(crust_Gamma2-1)/(crust_Gamma2-1)
  aaa(3) = ((1+aaa(2))*dd(2) + crust_K2*dd(2)**crust_Gamma2/(crust_Gamma2-1))/dd(2) - 1 &
    - crust_K3*dd(2)**(crust_Gamma3-1)/(crust_Gamma3-1)
  aaa(4) = ((1+aaa(3))*dd(3) + crust_K3*dd(3)**crust_Gamma3/(crust_Gamma3-1))/dd(3) - 1 &
    - crust_K4*dd(3)**(crust_Gamma4-1)/(crust_Gamma4-1)
  aaa(5) = ((1+aaa(4))*dd(4) + crust_K4*dd(4)**crust_Gamma4/(crust_Gamma4-1))/dd(4) - 1 &
    - K1*dd(4)**(Gamma1-1)/(Gamma1-1)
  aaa(6) = ((1+aaa(5))*dd(5) + K1*dd(5)**Gamma1/(Gamma1-1))/dd(5) - 1 &
    - K2*dd(5)**(Gamma2-1)/(Gamma2-1)
  aaa(7) = ((1+aaa(6))*dd(6) + K2*dd(6)**Gamma2/(Gamma2-1))/dd(6) - 1 &
    - K3*dd(6)**(Gamma3-1)/(Gamma3-1)


  write(*,*) 'rhonuc:'
  write(*,"(10es18.9e3)") dd
  write(*,*) 'cnuc:'
  write(*,"(10es18.9e3)") aaa
  write(*,*) 'kappa:'
  write(*,"(10es18.9e3)") crust_K1,crust_K2,crust_K3,crust_K4,K1,k2,k3
  write(*,*) 'for SACRA'
  write(*,"(20es18.9e3)") crust_Gamma1,crust_Gamma2,crust_Gamma3,crust_Gamma4, &
          Gamma1,Gamma2,Gamma3,dd*1.61930347d-18,crust_K1*1.61930347d-18**(1-crust_Gamma1)
  !stop "Only want to check parameters; L81"
  open(665,file='../'//trim(adjustl(eos_name))//'.dat',status="unknown")
  write(665,*) 5931
  
  log_rho = 1
  step = (log10(7.d15)-log_rho)/5931
  
  do ll=1,5931
    rho = 1.d1**(log_rho)
    if (log_rho.le.crust_devision1) then
      aa = aaa(1)
      G_practical = crust_Gamma1
      K_practical = crust_K1
    elseif (log_rho.le.crust_devision2) then
      aa = aaa(2)
      G_practical = crust_Gamma2
      K_practical = crust_K2
    elseif (log_rho.le.crust_devision3) then
      aa = aaa(3)
      G_practical = crust_Gamma3
      K_practical = crust_K3
    elseif (log_rho.le.rho_match) then
      aa = aaa(4)
      G_practical = crust_Gamma4
      K_practical = crust_K4
    elseif (log_rho.le.dens_division1) then
      aa = aaa(5)
      G_practical = Gamma1
      K_practical = K1
    elseif (log_rho.le.dens_division2) then
      aa = aaa(6)
      G_practical = Gamma2
      K_practical = K2
    else
      aa = aaa(7)
      G_practical = Gamma3
      K_practical = K3
    endif
    prs = K_practical*rho**G_practical
    eps = rho*(1+aa) + prs/(G_practical-1)
    cspead = sqrt(G_practical*prs/(eps+prs))
    hh  = (1+aa)+ K_practical*rho**(G_practical-1)*G_practical/(G_practical-1)
    !write(666,"(1P5E18.11)") eps,prs*clite**2,cspead,rho
    !write(665,"(1P5E18.11)") eps,prs*clite**2,cspead,rho
    
    ! for rns
    write(665,"(5es18.11)") eps, prs*clite**2, hh, rho/MB, cspead**2
    log_rho = log_rho + step
  enddo
  !close(666)
  close(665)

end program piecewisepoly

subroutine setup_table

  use parameters
  implicit none

  if(eos_name.eq.'PAL6') then
    Gamma1 = 2.227
    Gamma2 = 2.189
    Gamma3 = 2.159
    log_p1 = 34.380
  elseif(eos_name.eq.'SLy') then
    Gamma1 = 3.005
    Gamma2 = 2.988
    Gamma3 = 2.851
    log_p1 = 34.384
  elseif(eos_name.eq.'APR1') then
    Gamma1 = 2.442
    Gamma2 = 3.256
    Gamma3 = 2.908
    log_p1 = 33.943
  elseif(eos_name.eq.'APR2') then
    Gamma1 = 2.643
    Gamma2 = 3.014
    Gamma3 = 2.945
    log_p1 = 34.126
  elseif(eos_name.eq.'APR3') then
    Gamma1 = 3.166
    Gamma2 = 3.573
    Gamma3 = 3.281
    log_p1 = 34.392
  elseif(eos_name.eq.'APR4') then
    Gamma1 = 2.830
    Gamma2 = 3.445
    Gamma3 = 3.348
    log_p1 = 34.269
  elseif(eos_name.eq.'FPS') then
    Gamma1 = 2.985
    Gamma2 = 2.863
    Gamma3 = 2.600
    log_p1 = 34.283
  elseif(eos_name.eq.'WFF1') then
    Gamma1 = 2.519
    Gamma2 = 3.791
    Gamma3 = 3.660
    log_p1 = 34.031
  elseif(eos_name.eq.'WFF2') then
    Gamma1 = 2.888
    Gamma2 = 3.475
    Gamma3 = 3.517
    log_p1 = 34.233
  elseif(eos_name.eq.'WFF3') then
    Gamma1 = 3.329
    Gamma2 = 2.952
    Gamma3 = 2.589
    log_p1 = 34.283
  elseif(eos_name.eq.'BBB2') then
    Gamma1 = 3.418
    Gamma2 = 2.835
    Gamma3 = 2.832
    log_p1 = 34.331
  elseif(eos_name.eq.'ENG') then
    Gamma1 = 3.514
    Gamma2 = 3.130
    Gamma3 = 3.168
    log_p1 = 34.437
  elseif(eos_name.eq.'MPA1') then
    Gamma1 = 3.446
    Gamma2 = 3.572
    Gamma3 = 2.887
    log_p1 = 34.495
  elseif(eos_name.eq.'MS1') then ! 34.858 3.224 3.033 1.325
    Gamma1 = 3.224
    Gamma2 = 3.033
    Gamma3 = 1.325
    log_p1 = 34.858
  elseif(eos_name.eq.'MS2') then
    Gamma1 = 2.447
    Gamma2 = 2.184
    Gamma3 = 1.855
    log_p1 = 34.605
  elseif(eos_name.eq.'MS1b') then ! 34.855 3.456 3.011 1.425
    Gamma1 = 3.456
    Gamma2 = 3.011
    Gamma3 = 1.425
    log_p1 = 34.855
  elseif(eos_name.eq.'PS') then
    Gamma1 = 2.216
    Gamma2 = 1.640
    Gamma3 = 2.365
    log_p1 = 34.671
  elseif(eos_name.eq.'BGN1H1') then
    Gamma1 = 3.258
    Gamma2 = 1.472
    Gamma3 = 2.464
    log_p1 = 34.623
  elseif(eos_name.eq.'GNH3') then
    Gamma1 = 2.664
    Gamma2 = 2.194
    Gamma3 = 2.304
    log_p1 = 34.648
  elseif(eos_name.eq.'H1') then
    Gamma1 = 2.595
    Gamma2 = 1.845
    Gamma3 = 1.897
    log_p1 = 34.564
  elseif(eos_name.eq.'H2') then
    Gamma1 = 2.775
    Gamma2 = 1.855
    Gamma3 = 1.858
    log_p1 = 34.617
  elseif(eos_name.eq.'H4') then
    Gamma1 = 2.909
    Gamma2 = 2.246
    Gamma3 = 2.144
    log_p1 = 34.669
  elseif(eos_name.eq.'H3') then
    Gamma1 = 2.787
    Gamma2 = 1.951
    Gamma3 = 1.901
    log_p1 = 34.646
  elseif(eos_name.eq.'H5') then
    Gamma1 = 2.793
    Gamma2 = 1.974
    Gamma3 = 1.915
    log_p1 = 34.609
  elseif(eos_name.eq.'H7') then
    Gamma1 = 2.621
    Gamma2 = 2.048
    Gamma3 = 2.006
    log_p1 = 34.559
  elseif(eos_name.eq.'PCL2') then
    Gamma1 = 2.554
    Gamma2 = 1.880
    Gamma3 = 1.977
    log_p1 = 34.507
  elseif(eos_name.eq.'ALF1') then
    Gamma1 = 2.013
    Gamma2 = 3.389
    Gamma3 = 2.033
    log_p1 = 34.055
  elseif(eos_name.eq.'ALF2') then
    Gamma1 = 4.070
    Gamma2 = 2.411
    Gamma3 = 1.890
    log_p1 = 34.616
  elseif(eos_name.eq.'ALF3') then
    Gamma1 = 2.883
    Gamma2 = 2.653
    Gamma3 = 1.952
    log_p1 = 34.283
  elseif(eos_name.eq.'ALF4') then
    Gamma1 = 3.009
    Gamma2 = 3.438
    Gamma3 = 1.803
    log_p1 = 34.314
  elseif(eos_name.eq.'HB') then
    Gamma1 = 3.000
    Gamma2 = 3.000
    Gamma3 = 3.000
    log_p1 = 34.40364
  elseif(eos_name.eq.'H') then
    Gamma1 = 3.000
    Gamma2 = 3.000
    Gamma3 = 3.000
    log_p1 = 34.50364
  elseif(eos_name.eq.'B') then
    Gamma1 = 3.000
    Gamma2 = 3.000
    Gamma3 = 3.000
    log_p1 = 34.30364
  elseif(eos_name.eq.'15H') then
    Gamma1 = 3.000
    Gamma2 = 3.000
    Gamma3 = 3.000
    log_p1 = 34.70364
  elseif(eos_name.eq.'125H') then
    Gamma1 = 3.000
    Gamma2 = 3.000
    Gamma3 = 3.000
    log_p1 = 34.60364
  else
    stop 'no data for this eos!'
  endif


end subroutine setup_table
