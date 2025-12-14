module rotation_dispatch
  use para_mod
  !use rotation_const_j,  only: spin_gr_const_j
  use rotation_uniform,  only: spin
  use rotation_massless, only: spin_massless
  use rotation_gr,       only: spin_gr
  implicit none
contains

  subroutine call_rotation_solver()
    call spin()
  end subroutine call_rotation_solver

end module rotation_dispatch
