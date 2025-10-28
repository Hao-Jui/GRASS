module rotation_dispatch
  use para_mod
  use rotation_gr,       only: spin_gr
  use rotation_const_j,  only: spin_gr_const_j
  use rotation_st,       only: spin_st
  use rotation_massless, only: spin_massless
  implicit none
contains

  subroutine call_rotation_solver()
    select case (active_theory)
    case (THEORY_GR)
      call dispatch_gr_rotation()
    case (THEORY_ST)
      call dispatch_st_rotation()
    case default
      stop "call_rotation_solver: unknown theory"
    end select
  end subroutine call_rotation_solver

  subroutine dispatch_gr_rotation()
    select case (trim(adjustl(solver_type)))
    case ("uniform")
      call spin_gr()
    case ("const_j")
      call spin_gr_const_j()
    case default
      write(*,*) "Unknown GR rotation solver type:", trim(solver_type)
      stop 1
    end select
  end subroutine dispatch_gr_rotation

  subroutine dispatch_st_rotation()
    select case (trim(adjustl(solver_type)))
    case ("uniform")
      if (mphi_goal < 1.d-10) then
        call spin_massless()
      else
        call spin_st()
      end if
    case default
      write(*,*) "Unknown ST rotation solver type:", trim(solver_type)
      stop 1
    end select
  end subroutine dispatch_st_rotation

end module rotation_dispatch
