subroutine make_grid
  use para_mod
  implicit none
  integer m,s

  do s = 1, SDIV
    s_gp(s) = (dble(s)-1.d0) * DS
  enddo
  do m = 1, MDIV
    mu(m) = (dble(m)-1.d0) * DM
  enddo
  mu(MDIV) = 1.d0
  !open(72,file="./trig/s_axis.dat")
  !write(72,"(es18.9)") s_gp
  !open(73,file="./trig/mu_axis.dat")
  !write(73,"(es18.9)") mu
  !close(72); close(73)
end subroutine make_grid

subroutine GridTrig
  use toolkit_mod, only: legendre, plgndr
  use para_mod
  implicit none
  real(8) :: f2n(LMAX+1,SDIV), theta(MDIV)
  integer :: i,j,k,n,m
  real(8) :: sk, sj, sk1, sj1
  
  P_2n(:,:) = 0.d0
  P1_2n_1(:,:) = 0.d0
  sin_2n_1_theta(:,:) = 0.d0
  f_rho (:,:,:) = 0.d0
  f_gama(:,:,:) = 0.d0
  
  f2n(:,:) = 0.d0
  theta(:) = 0.d0

  !write(*,*) "compute_f2n"
  do n = 0, LMAX
    do i = 2, SDIV
      f2n(n+1,i) = ( (1.d0-s_gp(i)) / s_gp(i) )**(2*s_pwr*n)
    enddo
  enddo

  !write(*,*) "compute_f_rho_gamma"
  do j = 2, SDIV ! s
    do n = 1, LMAX
      do k = 2, SDIV ! s'
        sk = s_gp(k)
        sj = s_gp(j)
        sk1 = 1.d0 - sk
        sj1 = 1.d0 - sj
        if ( k < j ) then
          f_rho (j,n+1,k) = f2n(n+1,j) * (sj1 / sj)**s_pwr * sk**(s_pwr-1) / ( f2n(n+1,k) * sk1**(s_pwr+1) ) * dble(s_pwr)
          f_gama(j,n+1,k) = f2n(n+1,j) / ( f2n(n+1,k) * sk * sk1 ) * dble(s_pwr)
        else
          f_rho (j,n+1,k) = f2n(n+1,k) / ( f2n(n+1,j) * sk * sk1 ) * dble(s_pwr)
          f_gama(j,n+1,k) = f2n(n+1,k) * sj1**(2*s_pwr) * sk**(2*s_pwr-1) &
                            / ( sj**(2*s_pwr) * f2n(n+1,j) * sk1**(2*s_pwr+1) ) * dble(s_pwr)
        endif
      enddo
    enddo
  enddo
    
  ! > specify the coefficients for s = 0
  j = 1
  n = 0
  do k = 2, SDIV
      sk = s_gp(k)
      f_rho(j,n+1,k) = dble(s_pwr) / ( sk * (1.d0-sk) )
  enddo

  n = 1
  do k = 2, SDIV
      sk = s_gp(k)
      sk1 = 1.d0 - sk
      f_rho (j,n+1,k) = 0.d0
      f_gama(j,n+1,k) = dble(s_pwr) / ( sk * sk1 )
  enddo
    
  do n = 2, LMAX
    do k = 1, SDIV
      f_rho (j,n+1,k) = 0.d0
      f_gama(j,n+1,k) = 0.d0
    enddo
  enddo

  ! > specify the coefficients for s' = 0
  k = 1
  n = 0
  do j = 1, SDIV
      f_rho(j,n+1,k) = 0.d0
  enddo
  do j = 1, SDIV
    do n = 1, LMAX
      f_rho(j,n+1,k) = 0.d0
      f_gama(j,n+1,k)= 0.d0
    enddo
  enddo

  ! > specify the coefficients when n = 0 for f_rho
  n = 0
  do j = 2, SDIV
    do k = 2, SDIV
      sk = s_gp(k)
      sj = s_gp(j)
      sk1 = 1.d0 - sk
      sj1 = 1.d0 - sj

      if ( k < j ) then
        f_rho(j,n+1,k) = ( sj1 / sj )**s_pwr * sk**(s_pwr-1) / sk1**(s_pwr+1) * dble(s_pwr) 
      else
        f_rho(j,n+1,k) = dble(s_pwr) / ( sk * sk1 )
      endif
    enddo
  enddo
  

  ! > compute base functions
  n=0
  !write(fil1,"(i2)") n*2; open(100+n,file="./trig/P_"//trim(adjustl(fil1))//".dat")
  do i = 1, MDIV
    P_2n(i,n+1) = legendre(2*n,mu(i))
    !write(100+n,"(2es18.9)") mu(i), P_2n(i,n+1)
  enddo
  close(100+n)

  do m = 1, MDIV
    sin_theta(m) = sqrt(1.d0-mu(m)**2)
    theta(m) = asin(sin_theta(m))
  enddo

  do n = 1, LMAX
    !write(fil1,"(i2)") n*2
    !open(100+n,file="./trig/P_"//trim(adjustl(fil1))//".dat")
    !open(200+n,file="./trig/P1_"//trim(adjustl(fil1))//"_1.dat")
    do i = 1, MDIV
      P_2n          (i,n+1) = legendre(2*n,mu(i))
      P1_2n_1       (i,n+1) = plgndr(2*n-1,1,mu(i))
      sin_2n_1_theta(i,n  ) = sin(  (2.d0*n-1.d0)*theta(i) )
      !write(100+n,"(2es18.9)") mu(i), P_2n(i,n+1)
      !write(200+n,"(2es18.9)") mu(i), sin_2n_1_theta(i,n)
    enddo
  enddo
  close(100+n)
  close(200+n)

end subroutine GridTrig
