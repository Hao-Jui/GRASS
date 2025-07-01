subroutine Muller
  use para_mod
  implicit none
  real(8) :: p_at_e, h_at_p, n0_at_e, n0_at_h, e_at_h
  real(8) :: rho0,ee
  real(8) :: aaa, bbb
  integer :: pass = 0
  complex(8) :: ccx
  complex(8) :: cxin(3)
  complex(8) :: xa,xb,xc,qq,s(3),fs(3),s4(2),dd
  integer :: Jmax = 2500
  integer :: i
  real(8) :: arg,theta1
  complex(8) :: cu = dcmplx(0.d0, 1.d0)
  complex(8) :: alfa

  arg(alfa) = atan2( imag(alfa), real(alfa))

  r_ratio  = 1.d0
  e_center = 9.d14

  e_center = e_center * C * C * KSCALE
  p_center = p_at_e(e_center)
  h_center = h_at_p(p_center)
  call sphere

  write(*,*) " "
  output = .false.

  cxin(1) = dcmplx(h_center     ,r_ratio)
  cxin(2) = dcmplx(h_center*1.1 ,   0.95)
  cxin(3) = dcmplx(h_center*0.95,    0.7)

  do i = 1, 3
    s(i) = cxin(i)
    call fun_hc  ( real(s(i)), aaa)
    call fun_rerp( imag(s(i)), bbb)
    fs(i) = dcmplx(aaa, bbb)
    if (abs(fs(i)) < accuracy) then
      ccx=s(i)
      pass = 1
      exit
    endif
    if (fs(i).ne.fs(i)) then
      stop "func is not defined; L46 in Muller"
    endif
  enddo

  if (pass == 1) then 
  else
    do i = 1, Jmax
      qq = (s(3)-s(2)) / (s(2)-s(1))
      if (qq.ne.qq) then
        write(*,*) "Candidates:"
        write(*,*) s(1),s(2),s(3)
        write(*,*) "Difference:",abs(fs(1)),abs(fs(2)),abs(fs(3))
        write(*,*) "muller: no new generation",i
        stop
      endif
      xa = qq*fs(3) - qq*(1.+qq)*fs(2) + qq*qq*fs(1)
      xb = (2.*qq+1.)*fs(3) - (1.+qq)**2*fs(2) + qq*qq*fs(1)
      xc = (1.+qq)*fs(3)
      dd = xb*xb - 4.*xa*xc
      theta1 = arg(dd)
      !write(*,"(1P10E25.16)") theta1
      if (theta1 > pi/2.d0 .or. theta1 < -pi/2.d0) theta1 = theta1 + pi
    
      s4(1) = xb + sqrt(abs(dd))*exp(cu*theta1/2.d0)
      s4(2) = xb - sqrt(abs(dd))*exp(cu*theta1/2.d0)
      !write(*,"(1P10E15.6)") qq,a,b,c,dd
      !write(*,"(1P10E15.6)") theta1,s4(1),s4(2)
      !stop
      s(1) = s(2)
      s(2) = s(3)
      if (abs(s4(1)) > abs(s4(2))) then
        s(3) = s(2) - (s(2)-s(1)) * (2.d0*xc/s4(1))
      else
        s(3) = s(2) - (s(2)-s(1)) * (2.d0*xc/s4(2))
      endif

      fs(1)= fs(2)
      fs(2)= fs(3)
      call fun_hc  ( real(s(3)), aaa)
      call fun_rerp( imag(s(3)), bbb)
      fs(3) = dcmplx(aaa, bbb)
      
      write(*,"(A15,2es18.9,A10,es18.9)") "[hc, rep]",s(3),"er:",abs(fs(3))

      if (abs(fs(3)) < accuracy) then
        ccx = s(3)
        exit
      endif
      if (i==Jmax) stop "TOO many Iterations in Muller!"
    enddo
  endif
  h_center = real(ccx)
  r_ratio  = imag(ccx)

  output = .true.
  call const_j
  call mass_radius

  rho0 = n0_at_h(h_center)
  ee   = e_at_h (h_center)

  write(*,*) " ===================================="
  write(*,"(A10,i6)")         " iter :", i
  write(*,"(A10,es18.9,A10)") "rho_c :", rho0*MB,"g/cm^3"
  write(*,"(A10,es18.9,A10,es18.9)") "e_c   :", ee/(C * C * KSCALE),"g/cm^3", ee/(C * C * KSCALE)*rho_uni
  write(*,"(A10,es18.9)")     "rp/re :", r_ratio
  write(*,"(A10,es18.9,A4)")  "OMG_c :", omega_c/2.d0/pi,"Hz"
  write(*,"(A10,es18.9,A4)")  "Mass  :", Mass/MSUN,"M_o"
  write(*,"(A10,es18.9,A4)")  "M_0   :", Mass_0/MSUN,"M_o"
  write(*,"(A10,es18.9)")     "J     :", ang_mom
  write(*,"(A10,es18.9)")     "chi   :", chi
  write(*,"(A10,es18.9,A4)")  "R_is  :", r_e*sqrt(KAPPA)/1.d5,"km"
  write(*,"(A10,es18.9,A4)")  "R_cir :", r_circ/1e5,"km"
  write(*,"(A10,3es18.9)")    "er    :", fs(3)
  write(*,*) " ===================================="

end subroutine Muller