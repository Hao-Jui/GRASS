  subroutine loadEos

    use para_mod
    implicit none
    integer :: i

    open(16,file="./eos/"//trim(adjustl(eos_file))//".dat",status="unknown",form='formatted')
    read(16,*) num_tab
    allocate(log_e(num_tab),log_p(num_tab),log_h(num_tab),log_n0(num_tab))
    do i = 1, num_tab
      read(16,*) log_e(i),log_p(i),log_h(i),log_n0(i)
      if (i==1) then
        enthalpy_min =   log( log_h(i) )
      endif
      log_e(i) = log( log_e(i)*C*C*KSCALE)
      log_p(i) = log( log_p(i)*KSCALE)
      log_h(i) = log( log( log_h(i) ) )
      log_n0(i)= log( log_n0(i))
    enddo
    close(16)
    write(*,*) ' '
    write(*,*) ' '
    write(*,*) 'EOS data is in with log(h_min) =',enthalpy_min
    
  end subroutine loadEos


  real(8) function e_at_p(pp)
    use para_mod, only : log_p, log_e, num_tab
    implicit none
    real(8) pp,pwr
    call interp(log_p, log_e, num_tab, log(pp), pwr)
    e_at_p = exp(pwr)
  end function e_at_p

  real(8) function p_at_e(ee)
    use para_mod, only : log_p, log_e, num_tab
    implicit none
    real(8) ee,pwr
    call interp(log_e, log_p, num_tab, log(ee), pwr)
    p_at_e = exp(pwr)
  end function p_at_e

  real(8) function p_at_h(hh)
    use para_mod, only : log_p, log_h, num_tab
    implicit none
    real(8) hh,pwr
    call interp(log_h, log_p, num_tab, log(hh), pwr)
    p_at_h = exp(pwr)
  end function p_at_h

  real(8) function h_at_p(pp)
    use para_mod, only : log_p, log_h, num_tab
    implicit none
    real(8) pp,pwr
    call interp(log_p, log_h, num_tab, log(pp), pwr)
    h_at_p = exp(pwr)
  end function h_at_p

  real(8) function n0_at_e(ee)
    use para_mod, only : log_n0, log_e, num_tab
    implicit none
    real(8) ee,pwr
    call interp(log_e, log_n0, num_tab, log(ee), pwr)
    n0_at_e = exp(pwr)
  end function n0_at_e

  real(8) function n0_at_h(hh)
    use para_mod, only : log_n0, log_h, num_tab
    implicit none
    real(8) hh,pwr
    call interp(log_h, log_n0, num_tab, log(hh), pwr)
    n0_at_h = exp(pwr)
  end function n0_at_h

  real(8) function e_at_h(hh)
    use para_mod, only : log_e, log_h, num_tab
    implicit none
    real(8) hh,pwr
    call interp(log_h, log_e, num_tab, log(hh), pwr)
    e_at_h = exp(pwr)
  end function e_at_h
