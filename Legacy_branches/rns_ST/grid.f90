subroutine GridTrig

  use para_mod
  implicit none
  real(8) :: f2n(LMAX+1,SDIV), theta(MDIV)
  integer :: i,j,k,n,m
  real(8) :: sk, sj, sk1, sj1, r_gp
  real(8) :: legendre, plgndr
  real(8) :: besselk, besseli
  character(4) :: fil1, fil2
  
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
      f2n(n+1,i) = ( (1.d0-s_gp(i)) / s_gp(i) )**(2*n)
    enddo
  enddo

  !write(*,*) "compute_f_rho_gamma"
  if(SMAX .ne. 1.0) then
    do j = 2, SDIV
      do n = 1, LMAX
        do k = 2, SDIV
          sk = s_gp(k)
          sj = s_gp(j)
          sk1 = 1.d0 - sk
          sj1 = 1.d0 - sj
          if ( k < j ) then
            f_rho (j,n+1,k) = f2n(n+1,j) * sj1 / ( sj * f2n(n+1,k) * sk1**2 )
            f_gama(j,n+1,k) = f2n(n+1,j)/(f2n(n+1,k)*sk*sk1)
          else
            f_rho (j,n+1,k) = f2n(n+1,k) / ( f2n(n+1,j) * sk * sk1 )
            f_gama(j,n+1,k) = f2n(n+1,k)*sj1**2*sk/(sj**2*f2n(n+1,j)*sk1**3)
          endif
        enddo
      enddo
    enddo
    

    j = 1
    n = 0
    do k = 2, SDIV
      sk = s_gp(k)
      f_rho(j,n+1,k) = 1.d0 / ( sk*(1.d0-sk) )
    enddo

    n = 1
    do k = 2, SDIV
      sk = s_gp(k)
      sk1 = 1.d0 - sk
      f_rho (j,n+1,k) = 0.d0
      f_gama(j,n+1,k) = 1.d0/(sk*sk1)
    enddo
    
    do n=2,LMAX; do k=1,SDIV
      f_rho (j,n+1,k) = 0.d0
      f_gama(j,n+1,k) = 0.d0
    enddo; enddo

    k = 1
    n = 0
    do j = 1, SDIV
      f_rho(j,n+1,k) = 0.d0
    enddo
    do j = 1, SDIV; do n = 1, LMAX
      f_rho(j,n+1,k) = 0.d0
      f_gama(j,n+1,k)= 0.d0
    enddo; enddo

    n = 0
    do j = 2, SDIV; do k = 2, SDIV
      sk = s_gp(k)
      sj = s_gp(j)
      sk1 = 1.d0-sk
      sj1 = 1.d0-sj

      if ( k < j ) then
        f_rho(j,n+1,k) = sj1/(sj*sk1**2)
      else
        f_rho(j,n+1,k) = 1.d0/(sk*sk1)
      endif
    enddo; enddo
  else
    stop " SMAX is better to be less than 1; L110 grid.f90"
  endif

  n = 0
  do i = 1, MDIV
    P_2n(i,n+1) = legendre(2*n,mu(i))
  enddo

  do m = 1, MDIV
    sin_theta(m) = sqrt(1.d0-mu(m)**2)
    theta(m) = asin(sin_theta(m))
  enddo

  do n = 1, LMAX
    do i = 1, MDIV
      P_2n          (i,n+1) = legendre(2*n,mu(i))
      P1_2n_1       (i,n+1) = plgndr(2*n-1,1,mu(i))
      sin_2n_1_theta(i,n  ) = sin(  (2.d0*n-1.d0)*theta(i) )
    enddo
#if 0
    write(fil1,"(i2)") n*2
    open(100+n,file="./trig/mod_bessel_"//trim(adjustl(fil1))//".dat")
    do k = 1, SDIV
      sk = s_gp(k)
      r_gp= sk / (1.d0-sk)
      write(100+n,"(3es18.9)") r_gp, besselk(n*2,r_gp), besseli(n*2,r_gp)
    enddo
    close(100+n)
#endif
  enddo

end subroutine GridTrig

subroutine make_grid
  use para_mod
  implicit none
  integer :: m,s

  do s = 1, SDIV
    s_gp(s) = (s-1.0) * DS
  enddo
  do m = 1, MDIV
    mu(m) = (m-1.0) * DM
  enddo
  mu(MDIV) = 1.d0

end subroutine make_grid