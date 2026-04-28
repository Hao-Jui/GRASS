module eos_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use precision_mod, only: wp
  use para_mod, only: log_e, log_p, log_h, log_n0, num_tab, n_PT, phase_transition
  use ad_mod, only: dual, dual_const, operator(+), operator(-), operator(*), operator(/)
  implicit none
  private

  ! ------------------------------------------------------------------
  ! Module state: log-space tables (in para_mod) + PCHIP slopes here.
  ! Slopes are precomputed once in loadEos for every direction used by
  ! the public API. Each array stores d(out)/d(in) at table nodes.
  ! ------------------------------------------------------------------
  real(wp), allocatable, save :: e_tab(:), p_tab(:)
  real(wp), allocatable, save :: m_p_of_e(:)    ! d log_p / d log_e
  real(wp), allocatable, save :: m_e_of_p(:)    ! d log_e / d log_p
  real(wp), allocatable, save :: m_n0_of_e(:)   ! d log_n0 / d log_e
  real(wp), allocatable, save :: m_n0_of_h(:)   ! d log_n0 / d log_h
  real(wp), allocatable, save :: m_e_of_h(:)    ! d log_e  / d log_h
  real(wp), allocatable, save :: m_p_of_h(:)    ! d log_p  / d log_h
  real(wp), allocatable, save :: m_h_of_p(:)    ! d log_h  / d log_p
  real(wp), allocatable, save :: deriv_coeffs_work(:,:)

  public :: e_at_p, p_at_e, p_at_e_dual, n0_at_e, n0_at_h, e_at_h, p_at_h, h_at_p
  public :: loadEos, pressure_derivative_n, pe_at_h

contains

  ! ==================================================================
  ! Section 1: EOS table loader (also precomputes PCHIP slopes)
  ! ==================================================================
  subroutine loadEos
    use iso_fortran_env, only: output_unit
    use para_mod, only: p_at_PT, n_PT, eos_file, num_tab, log_e, log_p, log_h, log_n0, &
                        enthalpy_min, C, KSCALE
    implicit none
    integer :: i, unit, ios
    real(wp) :: dle, dlp
    real(wp), parameter :: GAMMA_PT = 6.e-2_wp  ! polytropic-index ceiling for a phase transition
    real(wp), parameter :: DLE_MIN  = 1.e-2_wp  ! minimum d(ln e) to count as sizable
    integer :: pt_buf(1024), pt_count

    pt_count = 0
    if (allocated(p_at_PT)) deallocate(p_at_PT)
    if (allocated(log_e))   deallocate(log_e, log_p, log_h, log_n0)
    if (allocated(e_tab))   deallocate(e_tab, p_tab)
    if (allocated(m_p_of_e))  deallocate(m_p_of_e, m_e_of_p, m_n0_of_e, &
                                         m_n0_of_h, m_e_of_h, m_p_of_h, m_h_of_p)

    open(newunit=unit, file="./eos/" // trim(adjustl(eos_file)) // ".dat", &
         status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "loadEos: failed to open EOS table"

    read(unit, *, iostat=ios) num_tab
    if (ios /= 0) error stop "loadEos: failed to read table size"

    allocate(log_e(num_tab), log_p(num_tab), log_h(num_tab), log_n0(num_tab))
    read(unit, *, iostat=ios) log_e(1), log_p(1), log_h(1), log_n0(1)
    if (ios /= 0) error stop "loadEos: malformed EOS row"
    enthalpy_min = log(log_h(1))
    if (log_h(1) <= 1._wp) error stop "loadEos: first enthalpy h(1) <= 1, log(log(h)) undefined"

    do i = 2, num_tab
      read(unit, *, iostat=ios) log_e(i), log_p(i), log_h(i), log_n0(i)
      if (ios /= 0) error stop "loadEos: malformed EOS row"
      dle = log(log_e(i) / log_e(i-1))
      dlp = log(log_p(i) / log_p(i-1))
      if (dle > DLE_MIN .and. abs(dlp / dle) < GAMMA_PT) then
        pt_count = pt_count + 1
        pt_buf(pt_count) = i
      end if
    end do
    close(unit)

    n_PT = pt_count
    if (n_PT > 0) then
      allocate(p_at_PT(n_PT))
      p_at_PT(:) = pt_buf(1:n_PT)
    else
      allocate(p_at_PT(0))
    end if

    log_e  = log(log_e * C * C * KSCALE)
    log_p  = log(log_p * KSCALE)
    log_h  = log(log(log_h))
    log_n0 = log(log_n0)
    allocate(e_tab(num_tab), p_tab(num_tab))
    e_tab = exp(log_e)
    p_tab = exp(log_p)

    ! PCHIP slopes for every (in, out) pair used by the public API.
    allocate(m_p_of_e(num_tab),  m_e_of_p(num_tab),  m_n0_of_e(num_tab), &
             m_n0_of_h(num_tab), m_e_of_h(num_tab),  m_p_of_h(num_tab), &
             m_h_of_p(num_tab))
    call pchip_slopes(log_e, log_p,  num_tab, m_p_of_e)
    call pchip_slopes(log_p, log_e,  num_tab, m_e_of_p)
    call pchip_slopes(log_e, log_n0, num_tab, m_n0_of_e)
    call pchip_slopes(log_h, log_n0, num_tab, m_n0_of_h)
    call pchip_slopes(log_h, log_e,  num_tab, m_e_of_h)
    call pchip_slopes(log_h, log_p,  num_tab, m_p_of_h)
    call pchip_slopes(log_p, log_h,  num_tab, m_h_of_p)

    write(output_unit, *) " "
    write(output_unit, *) "# EOS: ", eos_file
    if (n_PT > 0) then
      do i = 1, n_PT
        write(output_unit, "(A,i0,A,1X,i0)") " p_at_PT(", i, "):", p_at_PT(i)
      end do
      phase_transition = .true.
    end if
    write(output_unit, *) "EOS data is in with log(h_min) =", enthalpy_min

    if (any(ieee_is_nan(log_e)) .or. any(ieee_is_nan(log_p)) .or. any(ieee_is_nan(log_h))) &
      error stop "loadEos: wrong table"
  end subroutine loadEos

  ! ==================================================================
  ! Section 2: PCHIP slope kernel (Fritsch-Carlson, monotone)
  ! ------------------------------------------------------------------
  ! Interior slopes use the weighted harmonic mean of adjacent secants
  ! and collapse to zero when secants change sign or either is zero.
  ! That zero-slope branch is what handles flat phase-transition steps
  ! without needing a separate code path. Endpoints use the standard
  ! 3-point one-sided formula with monotonicity safeguards.
  ! ==================================================================
  pure subroutine pchip_slopes(x, y, n, m)
    integer, intent(in)  :: n
    real(wp), intent(in) :: x(n), y(n)
    real(wp), intent(out) :: m(n)
    integer :: i
    real(wp) :: h_a, h_b, d_a, d_b, w1, w2

    if (n <= 1) then
      m = 0._wp
      return
    end if
    if (n == 2) then
      m(1) = (y(2) - y(1)) / (x(2) - x(1))
      m(2) = m(1)
      return
    end if

    do i = 2, n - 1
      h_a = x(i)   - x(i-1)
      h_b = x(i+1) - x(i)
      d_a = (y(i)   - y(i-1)) / h_a
      d_b = (y(i+1) - y(i))   / h_b
      if (d_a * d_b <= 0._wp) then
        m(i) = 0._wp
      else
        w1 = 2._wp * h_b + h_a
        w2 = h_b + 2._wp * h_a
        m(i) = (w1 + w2) / (w1 / d_a + w2 / d_b)
      end if
    end do

    m(1) = endpoint_slope(x(2) - x(1), x(3) - x(2), &
                          (y(2) - y(1)) / (x(2) - x(1)), &
                          (y(3) - y(2)) / (x(3) - x(2)))
    m(n) = endpoint_slope(x(n) - x(n-1), x(n-1) - x(n-2), &
                          (y(n)   - y(n-1)) / (x(n)   - x(n-1)), &
                          (y(n-1) - y(n-2)) / (x(n-1) - x(n-2)))
  end subroutine pchip_slopes

  pure real(wp) function endpoint_slope(h_a, h_b, d_a, d_b) result(m_e)
    real(wp), intent(in) :: h_a, h_b, d_a, d_b
    m_e = ((2._wp * h_a + h_b) * d_a - h_a * d_b) / (h_a + h_b)
    if (m_e * d_a <= 0._wp) then
      m_e = 0._wp
    else if (d_a * d_b < 0._wp .and. abs(m_e) > 3._wp * abs(d_a)) then
      m_e = 3._wp * d_a
    end if
  end function endpoint_slope

  ! ==================================================================
  ! Section 3: cell finder (binary search with optional bracket walk)
  ! ==================================================================
  pure integer function find_cell(x, xb, idx_hint) result(i)
    real(wp), intent(in) :: x(:), xb
    integer,  intent(in), optional :: idx_hint
    integer :: lo, hi, mid, n, k, steps
    integer, parameter :: MAX_HINT_STEPS = 8

    n = size(x)
    if (n <= 1) then
      i = 1
      return
    end if
    if (xb <= x(1)) then
      i = 1
      return
    end if
    if (xb >= x(n)) then
      i = n - 1
      return
    end if

    if (present(idx_hint)) then
      k = min(n - 1, max(1, idx_hint))
      steps = 0
      do while (k > 1 .and. x(k) > xb .and. steps < MAX_HINT_STEPS)
        k = k - 1
        steps = steps + 1
      end do
      do while (k < n - 1 .and. x(k+1) <= xb .and. steps < MAX_HINT_STEPS)
        k = k + 1
        steps = steps + 1
      end do
      if (x(k) <= xb .and. xb < x(k+1)) then
        i = k
        return
      end if
    end if

    lo = 1
    hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (x(mid) <= xb) then
        lo = mid
      else
        hi = mid
      end if
    end do
    i = lo
  end function find_cell

  ! ==================================================================
  ! Section 4: Hermite basis + cell evaluators (scalar / pair / dual)
  ! ==================================================================
  pure subroutine hermite_basis(t, h00, h10, h01, h11)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: h00, h10, h01, h11
    real(wp) :: t2, t3
    t2 = t * t
    t3 = t2 * t
    h00 =  2._wp * t3 - 3._wp * t2 + 1._wp
    h10 =          t3 - 2._wp * t2 + t
    h01 = -2._wp * t3 + 3._wp * t2
    h11 =          t3 -          t2
  end subroutine hermite_basis

  pure subroutine hermite_basis_deriv(t, dh00, dh10, dh01, dh11)
    real(wp), intent(in)  :: t
    real(wp), intent(out) :: dh00, dh10, dh01, dh11
    real(wp) :: t2
    t2 = t * t
    dh00 =  6._wp * t2 - 6._wp * t
    dh10 =  3._wp * t2 - 4._wp * t + 1._wp
    dh01 = -6._wp * t2 + 6._wp * t
    dh11 =  3._wp * t2 - 2._wp * t
  end subroutine hermite_basis_deriv

  pure real(wp) function hermite_eval(x, y, m, i, xb) result(yb)
    real(wp), intent(in) :: x(:), y(:), m(:), xb
    integer,  intent(in) :: i
    real(wp) :: h, t, h00, h10, h01, h11
    h = x(i+1) - x(i)
    t = (xb - x(i)) / h
    call hermite_basis(t, h00, h10, h01, h11)
    yb = h00 * y(i) + h10 * h * m(i) + h01 * y(i+1) + h11 * h * m(i+1)
  end function hermite_eval

  pure subroutine hermite_eval_pair(x, y1, y2, m1, m2, i, xb, yb1, yb2)
    real(wp), intent(in)  :: x(:), y1(:), y2(:), m1(:), m2(:), xb
    integer,  intent(in)  :: i
    real(wp), intent(out) :: yb1, yb2
    real(wp) :: h, t, h00, h10, h01, h11
    h = x(i+1) - x(i)
    t = (xb - x(i)) / h
    call hermite_basis(t, h00, h10, h01, h11)
    yb1 = h00 * y1(i) + h10 * h * m1(i) + h01 * y1(i+1) + h11 * h * m1(i+1)
    yb2 = h00 * y2(i) + h10 * h * m2(i) + h01 * y2(i+1) + h11 * h * m2(i+1)
  end subroutine hermite_eval_pair

  pure subroutine hermite_eval_dual(x, y, m, i, xb, yb)
    real(wp),    intent(in)  :: x(:), y(:), m(:)
    integer,     intent(in)  :: i
    type(dual),  intent(in)  :: xb
    type(dual),  intent(out) :: yb
    real(wp) :: h, t, h00, h10, h01, h11, dh00, dh10, dh01, dh11
    real(wp) :: val, dval_dxb
    h = x(i+1) - x(i)
    t = (xb%val - x(i)) / h
    call hermite_basis(t, h00, h10, h01, h11)
    call hermite_basis_deriv(t, dh00, dh10, dh01, dh11)
    val      = h00 * y(i) + h10 * h * m(i) + h01 * y(i+1) + h11 * h * m(i+1)
    dval_dxb = (dh00 * y(i) + dh01 * y(i+1)) / h + dh10 * m(i) + dh11 * m(i+1)
    yb%val = val
    yb%der = dval_dxb * xb%der
  end subroutine hermite_eval_dual

  ! ==================================================================
  ! Section 5: log-space wrappers (input/output through exp/log)
  ! ==================================================================
  real(wp) function interp_eos(val_in, t_in, t_out, slope) result(out)
    real(wp), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out(:), slope(:)
    real(wp) :: log_in, log_out
    integer  :: i
    log_in  = log(val_in)
    i       = find_cell(t_in, log_in)
    log_out = hermite_eval(t_in, t_out, slope, i, log_in)
    out     = exp(log_out)
  end function interp_eos

  subroutine interp_eos_pair(val_in, t_in, t_out1, t_out2, slope1, slope2, &
                             val_out1, val_out2, idx_hint)
    real(wp), intent(in)  :: val_in
    real(wp), intent(in)  :: t_in(:), t_out1(:), t_out2(:), slope1(:), slope2(:)
    real(wp), intent(out) :: val_out1, val_out2
    integer,  intent(inout), optional :: idx_hint
    real(wp) :: log_in, log_out1, log_out2
    integer  :: i

    log_in = log(val_in)
    if (present(idx_hint)) then
      i = find_cell(t_in, log_in, idx_hint)
      idx_hint = i
    else
      i = find_cell(t_in, log_in)
    end if
    call hermite_eval_pair(t_in, t_out1, t_out2, slope1, slope2, i, log_in, &
                           log_out1, log_out2)
    val_out1 = exp(log_out1)
    val_out2 = exp(log_out2)
  end subroutine interp_eos_pair

  type(dual) function interp_eos_dual(val_in, t_in, t_out, slope) result(out)
    use ad_mod, only: log, exp
    type(dual), intent(in) :: val_in
    real(wp),   intent(in) :: t_in(:), t_out(:), slope(:)
    type(dual) :: log_in, log_out
    integer    :: i
    log_in = log(val_in)
    i      = find_cell(t_in, log_in%val)
    call hermite_eval_dual(t_in, t_out, slope, i, log_in, log_out)
    out    = exp(log_out)
  end function interp_eos_dual

  ! ==================================================================
  ! Section 6: Public API (thin wrappers over the log-space helpers)
  ! ==================================================================
  real(wp) function e_at_p(pp)
    real(wp), intent(in) :: pp
    e_at_p = interp_eos(pp, log_p, log_e, m_e_of_p)
  end function e_at_p

  real(wp) function p_at_e(ee)
    real(wp), intent(in) :: ee
    p_at_e = interp_eos(ee, log_e, log_p, m_p_of_e)
  end function p_at_e

  function p_at_e_dual(ee) result(res)
    type(dual), intent(in) :: ee
    type(dual) :: res
    res = interp_eos_dual(ee, log_e, log_p, m_p_of_e)
  end function p_at_e_dual

  real(wp) function n0_at_e(ee)
    real(wp), intent(in) :: ee
    n0_at_e = interp_eos(ee, log_e, log_n0, m_n0_of_e)
  end function n0_at_e

  real(wp) function n0_at_h(hh)
    real(wp), intent(in) :: hh
    n0_at_h = interp_eos(hh, log_h, log_n0, m_n0_of_h)
  end function n0_at_h

  real(wp) function e_at_h(hh)
    real(wp), intent(in) :: hh
    e_at_h = interp_eos(hh, log_h, log_e, m_e_of_h)
  end function e_at_h

  real(wp) function p_at_h(hh)
    real(wp), intent(in) :: hh
    p_at_h = interp_eos(hh, log_h, log_p, m_p_of_h)
  end function p_at_h

  real(wp) function h_at_p(pp)
    real(wp), intent(in) :: pp
    h_at_p = interp_eos(pp, log_p, log_h, m_h_of_p)
  end function h_at_p

  subroutine pe_at_h(hh, pp, ee, idx_hint)
    real(wp), intent(in)  :: hh
    real(wp), intent(out) :: pp, ee
    integer,  intent(inout), optional :: idx_hint
    call interp_eos_pair(hh, log_h, log_p, log_e, m_p_of_h, m_e_of_h, &
                         pp, ee, idx_hint)
  end subroutine pe_at_h

  ! ==================================================================
  ! Section 7: Fornberg high-order pressure derivative (n>=1)
  ! ------------------------------------------------------------------
  ! Independent of PCHIP. Uses on-the-fly Fornberg weights on raw
  ! e_tab/p_tab samples. Status flags: 1=bad order, 2=no EOS,
  ! 3=bad energy, 4=small stencil, 5=clamped to range.
  ! ==================================================================
  subroutine pressure_derivative_n(ee, n, derivative, status, idx_hint)
    implicit none
    real(wp), intent(in)  :: ee
    integer,  intent(in)  :: n
    real(wp), intent(out) :: derivative
    integer,  intent(out), optional :: status
    integer,  intent(inout), optional :: idx_hint

    integer :: info, idx, half_width, left, right, n_points
    real(wp) :: x0, min_e, max_e, ee_clamped
    logical :: fatal_error
    real(wp), parameter :: CLAMP_TOL = 1.e-12_wp

    info = 0
    fatal_error = .false.
    derivative = 0._wp

    if (n < 0) then
      info = 1; fatal_error = .true.
    else if (.not. allocated(log_e) .or. .not. allocated(log_p) .or. &
             .not. allocated(e_tab) .or. .not. allocated(p_tab) .or. num_tab <= 0) then
      info = 2; fatal_error = .true.
    else if (n >= num_tab) then
      info = 4; fatal_error = .true.
    else if (ee <= 0._wp) then
      info = 3; fatal_error = .true.
    end if

    if (fatal_error) then
      if (present(status)) status = info
      derivative = 0._wp
      return
    end if

    min_e = e_tab(1)
    max_e = e_tab(num_tab)
    ee_clamped = max(min_e, min(max_e, ee))
    if (abs(ee_clamped - ee) > CLAMP_TOL * max(1._wp, abs(ee_clamped))) info = 5

    x0  = ee_clamped
    idx = find_cell(log_e, log(x0), idx_hint)
    if (present(idx_hint)) idx_hint = idx

    half_width = max(n, 4)
    half_width = min(half_width, num_tab - 1)
    left  = idx - half_width
    right = idx + half_width
    if (left < 1) then
      right = min(num_tab, right + (1 - left))
      left  = 1
    end if
    if (right > num_tab) then
      left  = max(1, left - (right - num_tab))
      right = num_tab
    end if

    n_points = right - left + 1
    if (n_points <= n) then
      info = 4
      if (present(status)) status = info
      derivative = 0._wp
      return
    end if

    if (.not. allocated(deriv_coeffs_work)) then
      allocate(deriv_coeffs_work(n_points, n+1))
    else if (size(deriv_coeffs_work, 1) < n_points .or. size(deriv_coeffs_work, 2) < n+1) then
      deallocate(deriv_coeffs_work)
      allocate(deriv_coeffs_work(n_points, n+1))
    end if

    call fornberg_weights(x0, e_tab(left:right), n_points, n, &
                          deriv_coeffs_work(1:n_points, 1:n+1))
    derivative = dot_product(deriv_coeffs_work(1:n_points, n+1), p_tab(left:right))

    if (present(status)) status = info

  contains

    subroutine fornberg_weights(x0_local, x, m_pts, n_deriv, c)
      implicit none
      integer,  intent(in)  :: m_pts, n_deriv
      real(wp), intent(in)  :: x0_local
      real(wp), intent(in)  :: x(m_pts)
      real(wp), intent(out) :: c(m_pts, n_deriv + 1)
      integer  :: i, j, k, max_k
      real(wp) :: c1, c2, c3, c4, c5

      c(:, :) = 0._wp
      c(1,1) = 1._wp
      c1 = 1._wp
      c4 = x(1) - x0_local

      do i = 2, m_pts
        c2 = 1._wp
        c5 = c4
        c4 = x(i) - x0_local
        max_k = min(i - 1, n_deriv)
        do j = 1, i - 1
          c3 = x(i) - x(j)
          c2 = c2 * c3
          if (j == i - 1) then
            do k = max_k, 1, -1
              c(i, k+1) = c1 * (k * c(i-1, k) - c5 * c(i-1, k+1)) / c2
            end do
            c(i, 1) = -c1 * c5 * c(i-1, 1) / c2
          end if
          do k = max_k, 1, -1
            c(j, k+1) = (c4 * c(j, k+1) - k * c(j, k)) / c3
          end do
          c(j, 1) = c4 * c(j, 1) / c3
        end do
        c1 = c2
      end do
    end subroutine fornberg_weights
  end subroutine pressure_derivative_n

end module eos_mod
