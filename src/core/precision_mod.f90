module precision_mod
  use iso_fortran_env, only: real64, real128
  implicit none
  ! Default working precision. Set to real128 for quad builds.
  integer, parameter :: wp = real64
end module precision_mod
