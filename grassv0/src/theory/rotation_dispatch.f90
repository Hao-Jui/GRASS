module rotation_dispatch
  use para_mod
  use rotation_gr, only: spin_gr
  use rotation_st, only: spin_st
  implicit none
contains

  subroutine call_rotation_solver()
    select case (active_theory)
    case (THEORY_GR)
      call spin_gr()
    case (THEORY_ST)
      call spin_st()
    case default
      stop "call_rotation_solver: unknown theory"
    end select
  end subroutine call_rotation_solver

end module rotation_dispatch
