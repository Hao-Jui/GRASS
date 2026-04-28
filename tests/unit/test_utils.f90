module test_utils
  use precision_mod, only: wp
  implicit none
  integer :: n_pass = 0, n_fail = 0
contains
  subroutine assert_near(label, expected, actual, tol)
    character(*), intent(in) :: label
    real(wp), intent(in) :: expected, actual, tol
    if (abs(expected - actual) > tol) then
      write(*,'(A,A,A,es15.7,A,es15.7)') "  FAIL: ", label, &
        " expected=", expected, " got=", actual
      n_fail = n_fail + 1
    else
      n_pass = n_pass + 1
    end if
  end subroutine
  subroutine assert_rel_near(label, expected, actual, rtol)
    character(*), intent(in) :: label
    real(wp), intent(in) :: expected, actual, rtol
    real(wp) :: rel_err
    if (abs(expected) < tiny(expected)) then
      rel_err = abs(actual)
    else
      rel_err = abs((expected - actual) / expected)
    end if
    if (rel_err > rtol) then
      write(*,'(A,A,A,es15.7,A,es15.7,A,es10.2)') "  FAIL: ", label, &
        " expected=", expected, " got=", actual, " rel_err=", rel_err
      n_fail = n_fail + 1
    else
      n_pass = n_pass + 1
    end if
  end subroutine
  subroutine assert_true(label, condition)
    character(*), intent(in) :: label
    logical, intent(in) :: condition
    if (.not. condition) then
      write(*,'(A,A)') "  FAIL: ", label
      n_fail = n_fail + 1
    else
      n_pass = n_pass + 1
    end if
  end subroutine
  subroutine test_summary(suite_name)
    character(*), intent(in) :: suite_name
    write(*,'(A,A,A,I0,A,I0,A)') &
      "  ", suite_name, ": ", n_pass, " passed, ", n_fail, " failed"
    if (n_fail > 0) error stop "Tests failed"
    n_pass = 0; n_fail = 0
  end subroutine
end module test_utils
