module parameters
  implicit none
  
  integer, parameter :: piece = 6
  integer, parameter :: len = 5931
  real(8), parameter :: clite = 2.99792458d10
  character(32) :: eos_name = "APR4" 
  real(8) :: rho_(piece), gam(0:piece), kappa(0:piece), aa(0:piece)
end module parameters

program gen_table
  use parameters
  implicit none
  integer :: ll, jj
  real(8) :: log_rho, step
  real(8) :: eps, pres, vs2, rho_p
  
  call setup_table
  log_rho = 1.d0
  step = ( log10(7.d15) - log_rho ) / dble(5931)

  open( 665, file=trim(adjustl(eos_name))//'.dat' )
  write(665,*) 5931

  do ll = 1, len
      rho_p = 1.d1**(log_rho)
      do jj = 1, piece 
          if ( log_rho <= rho_(piece) ) exit
      enddo
      pres = kappa(jj-1) * rho_p**gam(jj-1)
      eps  = (1.d0 + aa(jj-1)) * rho_p + pres / ( gam(jj-1) - 1.d0 ) 
      vs2  = gam(jj-1) * pres / ( eps + pres )
      !write(*,"(20es18.9)") rho_p,eps,pres; stop
      write(665,"(5es18.11)") eps, pres*clite**2, 0.d0, gam(jj-1), gam(jj-1)
      log_rho = log_rho + step
  enddo
  close(665)
end program

subroutine setup_table
  use parameters
  implicit none
  integer :: ll
  if ( eos_name == 'BSk21' ) then
      kappa(0) = 1.d1**(-12.4958d0) 
      gam(0) = 1.6357d0; aa  (0) = 0.d0 
      gam(1) = 1.3107d0; rho_(1) = 6.9433d0 
      gam(2) = 0.7452d0; rho_(2) = 11.3651d0
      gam(3) = 1.2571d0; rho_(3) = 12.3329d0
      gam(4) = 3.4841d0; rho_(4) = 14.1610d0
      gam(5) = 3.1032d0; rho_(5) = 14.6921d0
      gam(6) = 2.8012d0; rho_(6) = 14.9021d0
  elseif ( eos_name == 'APR4' ) then
      kappa(0) = 6.80110d-9
      gam(0) = 1.58425d0; aa  (0) = 0.d0 
      gam(1) = 1.28733d0; rho_(1) = log10(2.44034d7)
      gam(2) = 0.62223d0; rho_(2) = log10(3.78358d11)
      gam(3) = 1.35692d0; rho_(3) = log10(2.62780d12)
      gam(4) = 2.830d0; rho_(4) = log10(1.511879208d14)
      gam(5) = 3.445d0; rho_(5) = 14.7d0
      gam(6) = 3.348d0; rho_(6) = 15.d0
  endif

  do ll = 1, piece
      kappa(ll) = kappa(ll-1) * ( 1.d1**rho_(ll) )**( gam(ll-1) - gam(ll) )
      aa   (ll) = aa   (ll-1) &
        + kappa(ll-1) / ( gam(ll-1) - 1.d0 ) * (1.d1**rho_(ll))**( gam(ll-1) - 1.d0 ) &
        - kappa(ll)   / ( gam(ll)   - 1.d0 ) * (1.d1**rho_(ll))**( gam(ll)   - 1.d0 )
  enddo

  write(*,*) 'rhonuc:'
  write(*,"(10es18.9e3)") 1.d1**rho_
  write(*,*) 'cnuc:'
  write(*,"(10es18.9e3)") aa
  write(*,*) 'kappa:'
  write(*,"(10es18.9e3)") kappa
  write(*,*) " "

end subroutine
