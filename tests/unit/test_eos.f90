program test_eos
  use precision_mod, only: wp
  use test_utils
  use para_mod, only: num_tab, eos_file, log_e
  use eos_mod, only: loadEos, p_at_e
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(wp) :: e_mid, p_mid
  integer :: mid_idx

  ! -----------------------------------------------------------------
  ! Setup: set EOS filename and load table.
  ! The working directory must be the GRASS root so ./eos/MPA1.dat
  ! is accessible.
  ! -----------------------------------------------------------------
  eos_file = "MPA1"
  call loadEos

  ! -----------------------------------------------------------------
  ! Test 1: Table loaded with a positive number of entries
  ! -----------------------------------------------------------------
  call assert_true("num_tab > 0", num_tab > 0)

  ! -----------------------------------------------------------------
  ! Test 2: p_at_e at a midpoint is finite and positive
  ! -----------------------------------------------------------------
  mid_idx = num_tab / 2
  e_mid   = exp(log_e(mid_idx))
  p_mid   = p_at_e(e_mid)
  call assert_true("p_at_e midpoint is finite",   ieee_is_finite(p_mid))
  call assert_true("p_at_e midpoint is positive", p_mid > 0.0_wp)

  call test_summary("eos_mod")
end program test_eos
