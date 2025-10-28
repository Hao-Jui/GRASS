subroutine loadEos
#include "option_macro.h"
    use para_mod
    implicit none
    integer :: i

    p_at_PT = 1
    open(16,file="./eos/"//trim(adjustl(eos_file))//".dat",status="unknown",form='formatted')
    read(16,*) num_tab
    allocate(log_e(num_tab),log_p(num_tab),log_h(num_tab),log_n0(num_tab))
    do i = 1, num_tab
      read(16,*) log_e(i), log_p(i), log_h(i), log_n0(i)
      if ( i == 1 ) then
        enthalpy_min =   log( log_h(i) )
      else
        if ( log_p(i)/log_p(i-1) < 1.d0 + 1.d-15 ) then
          p_at_PT = i
          !write(*,*) log( log_p(i-1)*KSCALE ), log( log_p(i)*KSCALE ) 
          !write(*,*) log( log_e(i-1)*C*C*KSCALE ), log( log_e(i)*C*C*KSCALE ) 
        endif
      endif
    enddo
    close(16)
    
    log_e(:) = log( log_e(:)*C*C*KSCALE)
    log_p(:) = log( log_p(:)*KSCALE)
    log_h(:) = log( log( log_h(:) ) )
    log_n0(:)= log( log_n0(:))

    write(*,*) ' '
    write(*,*) ' '
    write(*,*) '# EOS: ', eos_file
    write(*,*) 'EOS data is in with log(h_min) =',enthalpy_min
    
    if ( any(isnan(log_e)) .or. any(isnan(log_p)) .or. any(isnan(log_h)) ) stop "wrong table"

  end subroutine loadEos


  real(8) function e_at_p(pp)
#include "option_macro.h"
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_e, num_tab
    implicit none
    real(8) pp,pwr
#if defined(PT)
    call interp_pt(log_p, log_e, num_tab, log(pp), pwr)
#else
    call interp(log_p, log_e, num_tab, log(pp), pwr)
#endif
    e_at_p = exp(pwr)
  end function e_at_p

  real(8) function p_at_e(ee)
#include "option_macro.h"
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_e, num_tab
    implicit none
    real(8) ee,pwr
#if defined(PT)
    call interp_pt(log_e, log_p, num_tab, log(ee), pwr)
#else
    call interp(log_e, log_p, num_tab, log(ee), pwr)
#endif
    p_at_e = exp(pwr)
  end function p_at_e

  type(dual) function p_at_e_dual(ee)
#include "option_macro.h"
    use toolkit_mod, only: interp_dual
    use ad_mod, only: dual, log, exp
    use para_mod, only : log_p, log_e, num_tab
    implicit none
    type(dual), intent(in) :: ee
    type(dual) :: pwr

    call interp_dual(log_e, log_p, num_tab, log(ee), pwr)
    p_at_e_dual = exp(pwr)
  end function p_at_e_dual

  real(8) function n0_at_e(ee)
#include "option_macro.h"
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_n0, log_e, num_tab
    implicit none
    real(8) ee,pwr
#if defined(PT)
    call interp_pt(log_e, log_n0, num_tab, log(ee), pwr)
#else
    call interp(log_e, log_n0, num_tab, log(ee), pwr)
#endif
    n0_at_e = exp(pwr)
  end function n0_at_e

  real(8) function n0_at_h(hh)
#include "option_macro.h"
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_n0, log_h, num_tab
    implicit none
    real(8) hh,pwr
#if defined(PT)
    call interp_pt(log_h, log_n0, num_tab, log(hh), pwr)
#else
    call interp(log_h, log_n0, num_tab, log(hh), pwr)
#endif
    n0_at_h = exp(pwr)
  end function n0_at_h

  real(8) function e_at_h(hh)
#include "option_macro.h"
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_e, log_h, num_tab
    implicit none
    real(8) hh,pwr
#if defined(PT)
    call interp_pt(log_h, log_e, num_tab, log(hh), pwr)
#else
    call interp(log_h, log_e, num_tab, log(hh), pwr)
#endif
    e_at_h = exp(pwr)
  end function e_at_h

  real(8) function p_at_h(hh)
#include "option_macro.h"
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_h, num_tab
    implicit none
    real(8) hh,pwr
    call interp(log_h, log_p, num_tab, log(hh), pwr)
    p_at_h = exp(pwr)
  end function p_at_h

  real(8) function h_at_p(pp)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_h, num_tab
    implicit none
    real(8) pp,pwr
    call interp(log_p, log_h, num_tab, log(pp), pwr)
    h_at_p = exp(pwr)
  end function h_at_p
