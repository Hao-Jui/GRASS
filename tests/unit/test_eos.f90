program test_eos
  use precision_mod, only: wp
  use test_utils
  use para_mod, only: num_tab, eos_file, log_e, log_h, log_p
  use eos_mod, only: loadEos, p_at_e, e_at_p, p_at_e_dual, &
                     n0_at_e, n0_at_h, e_at_h, p_at_h, h_at_p, &
                     pe_at_h, pressure_derivative_n
  use ad_mod, only: dual
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  real(wp) :: e_mid, p_mid, h_mid, e_off
  real(wp) :: e_rt, p_rt, h_rt
  real(wp) :: e_from_h, p_from_h
  real(wp) :: pe_p, pe_e
  real(wp) :: n0_e, n0_h
  real(wp) :: dp_fd, dp_pdn, delta
  type(dual) :: e_dual, p_dual
  integer :: mid_idx

  ! Setup: load EOS table
  eos_file = "MPA1"
  call loadEos

  ! Pick mid-table test points
  mid_idx = num_tab / 2
  e_mid = exp(log_e(mid_idx))
  p_mid = p_at_e(e_mid)
  h_mid = exp(log_h(mid_idx))

  ! --- Original tests (3 assertions) ---

  ! 1. Table loaded with positive entries
  call assert_true("num_tab > 0", num_tab > 0)

  ! 2. p_at_e at midpoint is finite
  call assert_true("p_at_e midpoint is finite", ieee_is_finite(p_mid))

  ! 3. p_at_e at midpoint is positive
  call assert_true("p_at_e midpoint is positive", p_mid > 0.0_wp)

  ! --- Round-trip consistency tests (5 assertions) ---

  ! 4. p -> e -> p round-trip
  e_rt = e_at_p(p_mid)
  p_rt = p_at_e(e_rt)
  call assert_rel_near("p->e->p round-trip", p_mid, p_rt, 1e-10_wp)

  ! 5. e -> p -> e round-trip
  p_rt = p_at_e(e_mid)
  e_rt = e_at_p(p_rt)
  call assert_rel_near("e->p->e round-trip", e_mid, e_rt, 1e-10_wp)

  ! 6. e_at_h is finite and positive
  e_from_h = e_at_h(h_mid)
  call assert_true("e_at_h finite+positive", &
    ieee_is_finite(e_from_h) .and. e_from_h > 0.0_wp)

  ! 7. p_at_h matches p_at_e(e_at_h(h))
  p_from_h = p_at_h(h_mid)
  call assert_rel_near("p_at_h consistency", p_at_e(e_at_h(h_mid)), p_from_h, 1e-10_wp)

  ! 8. h -> p -> h round-trip
  h_rt = h_at_p(p_at_h(h_mid))
  call assert_rel_near("h->p->h round-trip", h_mid, h_rt, 1e-10_wp)

  ! --- Bundle function test (2 assertions) ---

  ! 9-10. pe_at_h matches individual lookups
  call pe_at_h(h_mid, pe_p, pe_e)
  call assert_rel_near("pe_at_h pressure", p_at_h(h_mid), pe_p, 1e-12_wp)
  call assert_rel_near("pe_at_h energy", e_at_h(h_mid), pe_e, 1e-12_wp)

  ! --- Dual-number tests (2 assertions) ---
  ! Use an off-grid point: barycentric interpolation returns zero derivative
  ! at exact table nodes (L'Hopital shortcut in eos_mod.f90:477-479).
  e_off = exp(0.5_wp * (log_e(mid_idx) + log_e(mid_idx + 1)))

  ! 11. Dual real part matches scalar (off-grid point)
  e_dual = dual(e_off, 0.0_wp)
  p_dual = p_at_e_dual(e_dual)
  call assert_rel_near("dual real part", p_at_e(e_off), p_dual%val, 1e-14_wp)

  ! 12. Dual derivative matches finite difference (off-grid point)
  delta = e_off * 1.0e-6_wp
  dp_fd = (p_at_e(e_off + delta) - p_at_e(e_off - delta)) / (2.0_wp * delta)
  e_dual = dual(e_off, 1.0_wp)
  p_dual = p_at_e_dual(e_dual)
  call assert_rel_near("dual derivative vs FD", dp_fd, p_dual%der, 1.0e-4_wp)

  ! --- Number density tests (2 assertions) ---

  ! 13. n0_at_e positive
  n0_e = n0_at_e(e_mid)
  call assert_true("n0_at_e positive", n0_e > 0.0_wp)

  ! 14. n0_at_h positive
  n0_h = n0_at_h(h_mid)
  call assert_true("n0_at_h positive", n0_h > 0.0_wp)

  ! --- Derivative test (1 assertion) ---

  ! 15. pressure_derivative_n order 1 vs finite difference
  delta = e_mid * 1.0e-6_wp
  dp_fd = (p_at_e(e_mid + delta) - p_at_e(e_mid - delta)) / (2.0_wp * delta)
  call pressure_derivative_n(e_mid, 1, dp_pdn)
  call assert_rel_near("dP/de order 1 vs FD", dp_fd, dp_pdn, 1.0e-4_wp)

  ! --- Table edge tests (3 assertions) ---

  ! 16. Low-energy edge
  call assert_true("p_at_e low-e finite", &
    ieee_is_finite(p_at_e(exp(log_e(2)))))

  ! 17. High-energy edge
  call assert_true("p_at_e high-e finite", &
    ieee_is_finite(p_at_e(exp(log_e(num_tab - 1)))))

  ! 18. Low-enthalpy edge
  call assert_true("e_at_h low-h finite", &
    ieee_is_finite(e_at_h(exp(log_h(2)))))

  call test_summary("eos_mod")
end program test_eos
