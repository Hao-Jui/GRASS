module ad_mod
  use precision_mod, only: wp
  ! Dual-number forward-mode automatic differentiation.
  ! A dual value stores (f, f'), and overloaded arithmetic/exp/log
  ! propagate derivatives automatically.
  implicit none

  type :: dual
    real(wp) :: val
    real(wp) :: der
  end type dual

  interface operator(+)
    module procedure add_dd, add_dr, add_rd
  end interface

  interface operator(-)
    module procedure sub_dd, sub_dr, sub_rd, neg_d
  end interface

  interface operator(*)
    module procedure mul_dd, mul_dr, mul_rd
  end interface

  interface operator(/)
    module procedure div_dd, div_dr, div_rd
  end interface

  interface exp
    module procedure exp_d
  end interface

  interface log
    module procedure log_d
  end interface

contains

  pure function dual_const(x) result(d)
    real(wp), intent(in) :: x
    type(dual) :: d
    d%val = x
    d%der = 0._wp
  end function dual_const

  pure function dual_var(x) result(d)
    real(wp), intent(in) :: x
    type(dual) :: d
    d%val = x
    d%der = 1._wp
  end function dual_var

  pure function add_dd(a, b) result(c)
    type(dual), intent(in) :: a, b
    type(dual) :: c
    c%val = a%val + b%val
    c%der = a%der + b%der
  end function add_dd

  pure function add_dr(a, b) result(c)
    type(dual), intent(in) :: a
    real(wp),  intent(in) :: b
    type(dual) :: c
    c%val = a%val + b
    c%der = a%der
  end function add_dr

  pure function add_rd(a, b) result(c)
    real(wp),  intent(in) :: a
    type(dual), intent(in) :: b
    type(dual) :: c
    c%val = a + b%val
    c%der = b%der
  end function add_rd

  pure function sub_dd(a, b) result(c)
    type(dual), intent(in) :: a, b
    type(dual) :: c
    c%val = a%val - b%val
    c%der = a%der - b%der
  end function sub_dd

  pure function sub_dr(a, b) result(c)
    type(dual), intent(in) :: a
    real(wp),  intent(in) :: b
    type(dual) :: c
    c%val = a%val - b
    c%der = a%der
  end function sub_dr

  pure function sub_rd(a, b) result(c)
    real(wp),  intent(in) :: a
    type(dual), intent(in) :: b
    type(dual) :: c
    c%val = a - b%val
    c%der = -b%der
  end function sub_rd

  pure function neg_d(a) result(c)
    type(dual), intent(in) :: a
    type(dual) :: c
    c%val = -a%val
    c%der = -a%der
  end function neg_d

  pure function mul_dd(a, b) result(c)
    type(dual), intent(in) :: a, b
    type(dual) :: c
    c%val = a%val * b%val
    c%der = a%der * b%val + a%val * b%der
  end function mul_dd

  pure function mul_dr(a, b) result(c)
    type(dual), intent(in) :: a
    real(wp),  intent(in) :: b
    type(dual) :: c
    c%val = a%val * b
    c%der = a%der * b
  end function mul_dr

  pure function mul_rd(a, b) result(c)
    real(wp),  intent(in) :: a
    type(dual), intent(in) :: b
    type(dual) :: c
    c%val = a * b%val
    c%der = a * b%der
  end function mul_rd

  pure function div_dd(a, b) result(c)
    type(dual), intent(in) :: a, b
    type(dual) :: c
    real(wp) :: inv
    inv = 1._wp / b%val
    c%val = a%val * inv
    c%der = (a%der - c%val * b%der) * inv
  end function div_dd

  pure function div_dr(a, b) result(c)
    type(dual), intent(in) :: a
    real(wp),  intent(in) :: b
    type(dual) :: c
    real(wp) :: inv
    inv = 1._wp / b
    c%val = a%val * inv
    c%der = a%der * inv
  end function div_dr

  pure function div_rd(a, b) result(c)
    real(wp),  intent(in) :: a
    type(dual), intent(in) :: b
    type(dual) :: c
    real(wp) :: inv
    inv = 1._wp / b%val
    c%val = a * inv
    c%der = -c%val * b%der * inv
  end function div_rd

  pure function exp_d(a) result(c)
    type(dual), intent(in) :: a
    type(dual) :: c
    real(wp) :: ev
    ev = exp(a%val)
    c%val = ev
    c%der = ev * a%der
  end function exp_d

  pure function log_d(a) result(c)
    type(dual), intent(in) :: a
    type(dual) :: c
    c%val = log(a%val)
    c%der = a%der / a%val
  end function log_d

end module ad_mod
