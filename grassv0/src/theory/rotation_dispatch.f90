module rotation_dispatch
  use para_mod
  !use rotation_const_j,  only: spin_gr_const_j
  use rotation_uniform,  only: spin
  use rotation_massless, only: spin_massless
  use rotation_gr,       only: spin_gr
  implicit none
contains

  subroutine call_rotation_solver()
    select case (trim(adjustl(solver_type)))
    case ("uniform")
      if ( active_theory == THEORY_ST .and. mphi_goal < 1.d-10 ) then
        call spin_massless()
      else
        call spin()
      endif
    case ("const_j")
      !call spin_gr_const_j()
    case default
      write(*,*) "Unknown GR rotation solver type:", trim(solver_type)
      stop 1
    end select
  end subroutine call_rotation_solver

end module rotation_dispatch
