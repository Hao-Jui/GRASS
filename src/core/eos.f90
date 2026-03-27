module eos_mod
  use precision_mod, only: wp
  use toolkit_mod, only: interp, interp_pt, interp_dual, interp_pt_dual
  use para_mod, only: log_e, log_p, log_h, log_n0, num_tab, n_PT, phase_transition
  use ad_mod, only: dual
  implicit none
  private
  ! Public elemental functions
  public :: e_at_p, p_at_e, p_at_e_dual, n0_at_e, n0_at_h, e_at_h, p_at_h, h_at_p
  ! Public subroutines
  public :: loadEos, pressure_derivative_n
contains

  real(wp) function interp_eos_with_phase(val_in, t_in, t_out)
    real(wp), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out(:)
    real(wp) :: res_log
    if (phase_transition) then
      call interp_pt(t_in, t_out, num_tab, log(val_in), res_log)
    else
      call interp(t_in, t_out, num_tab, log(val_in), res_log)
    end if
    interp_eos_with_phase = exp(res_log)
  end function interp_eos_with_phase

  real(wp) function interp_eos_simple(val_in, t_in, t_out)
    real(wp), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out(:)
    real(wp) :: res_log
    if (phase_transition) then
      call interp_pt(t_in, t_out, num_tab, log(val_in), res_log)
    else
      call interp(t_in, t_out, num_tab, log(val_in), res_log)
    end if
    interp_eos_simple = exp(res_log)
  end function interp_eos_simple

  type(dual) function interp_eos_dual(val_in, t_in, t_out)
    use ad_mod, only: log, exp
    type(dual), intent(in) :: val_in
    real(wp), intent(in) :: t_in(:), t_out(:)
    type(dual) :: res_log
    if (phase_transition) then
      call interp_pt_dual(t_in, t_out, num_tab, log(val_in), res_log)
    else
      call interp_dual(t_in, t_out, num_tab, log(val_in), res_log)
    end if
    interp_eos_dual = exp(res_log)
  end function interp_eos_dual

  subroutine loadEos
    use iso_fortran_env, only: output_unit
    use para_mod, only: p_at_PT, n_PT, eos_file, num_tab, log_e, log_p, log_h, log_n0, enthalpy_min, C, KSCALE
    implicit none
    integer :: i, unit, ios
    real(wp) :: dle, dlp
    real(wp), parameter :: GAMMA_PT = 6.e-2_wp  ! polytropic-index ceiling for a phase transition
    real(wp), parameter :: DLE_MIN  = 1.e-2_wp  ! minimum d(ln e) to count as sizable
    integer :: pt_buf(1024), pt_count

    pt_count = 0
    if (allocated(p_at_PT)) deallocate(p_at_PT)
    if (allocated(log_e)) deallocate(log_e, log_p, log_h, log_n0)

    open(newunit=unit, file="./eos/" // trim(adjustl(eos_file)) // ".dat", status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "loadEos: failed to open EOS table"

    read(unit, *, iostat=ios) num_tab
    if (ios /= 0) error stop "loadEos: failed to read table size"

    allocate(log_e(num_tab), log_p(num_tab), log_h(num_tab), log_n0(num_tab))
    read(unit, *, iostat=ios) log_e(1), log_p(1), log_h(1), log_n0(1)
    if (ios /= 0) error stop "loadEos: malformed EOS row"
    enthalpy_min = log(log_h(1))
    if (log_h(1) <= 1._wp) error stop "loadEos: first enthalpy h(1) <= 1 — log(log(h)) is undefined"

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

    log_e = log(log_e * C * C * KSCALE)
    log_p = log(log_p * KSCALE)
    log_h = log(log(log_h))
    log_n0 = log(log_n0)

    write(output_unit, *) " "
    write(output_unit, *) "# EOS: ", eos_file
    if (n_PT > 0) then
      do i = 1, n_PT
        write(output_unit, "(A,i0,A,X,i0)") " p_at_PT(", i, "):", p_at_PT(i)
      end do
      phase_transition = .true.
    end if
    write(output_unit, *) "EOS data is in with log(h_min) =", enthalpy_min

    if (any(isnan(log_e)) .or. any(isnan(log_p)) .or. any(isnan(log_h))) error stop "loadEos: wrong table"
  end subroutine loadEos


  real(wp) function e_at_p(pp)
    real(wp), intent(in) :: pp
    e_at_p = interp_eos_with_phase(pp, log_p, log_e)
  end function e_at_p

  real(wp) function p_at_e(ee)
    real(wp), intent(in) :: ee
    p_at_e = interp_eos_with_phase(ee, log_e, log_p)
  end function p_at_e

  function p_at_e_dual(ee) result(res)
    type(dual), intent(in) :: ee
    type(dual) :: res
    res = interp_eos_dual(ee, log_e, log_p)
  end function p_at_e_dual

  real(wp) function n0_at_e(ee)
    real(wp), intent(in) :: ee
    n0_at_e = interp_eos_with_phase(ee, log_e, log_n0)
  end function n0_at_e

  real(wp) function n0_at_h(hh)
    real(wp), intent(in) :: hh
    n0_at_h = interp_eos_simple(hh, log_h, log_n0)
  end function n0_at_h

  real(wp) function e_at_h(hh)
    real(wp), intent(in) :: hh
    e_at_h = interp_eos_simple(hh, log_h, log_e)
  end function e_at_h

  real(wp) function p_at_h(hh)
    real(wp), intent(in) :: hh
    p_at_h = interp_eos_simple(hh, log_h, log_p)
  end function p_at_h

  real(wp) function h_at_p(pp)
    real(wp), intent(in) :: pp
    h_at_p = interp_eos_simple(pp, log_p, log_h)
  end function h_at_p

  ! **********************************************************************
  ! Selects a stencil around the requested energy clamps out-of-range 
  !  queries, and evaluates d^n p / d e^n
  !    ifail  (output) - Error flag (1=invalid order, 2=EOS unavailable, 
  !                     3=nonpositive energy, 4=insufficient stencil, 
  !                     5=input clamped to table range)
  ! **********************************************************************
  subroutine pressure_derivative_n(ee, n, derivative, status)
    use para_mod, only: log_e, log_p, num_tab
    implicit none
    real(wp), intent(in) :: ee
    integer, intent(in) :: n
    real(wp), intent(out) :: derivative
    integer, intent(out), optional :: status

    integer :: info, idx, half_width, left, right, n_points, i
    real(wp) :: x0, min_e, max_e
    real(wp), allocatable :: nodes(:), values(:), coeffs(:, :)
    real(wp) :: ee_clamped
    logical :: fatal_error
    real(wp), parameter :: clamp_tol = 1.d-12

    info = 0
    fatal_error = .false.
    derivative = 0.d0

    if (n < 0) then
      info = 1
      fatal_error = .true.
    else if (.not. allocated(log_e) .or. .not. allocated(log_p) .or. num_tab <= 0) then
      info = 2
      fatal_error = .true.
    else if (n >= num_tab) then
      info = 4
      fatal_error = .true.
    else if (ee <= 0.d0) then
      info = 3
      fatal_error = .true.
    end if

    if (fatal_error) then
      if (present(status)) status = info
      derivative = 0.d0
      return
    end if

    min_e = exp(log_e(1))
    max_e = exp(log_e(num_tab))

    ee_clamped = max(min_e, min(max_e, ee))
    if (abs(ee_clamped - ee) > clamp_tol * max(1.d0, abs(ee_clamped))) info = 5

    x0 = ee_clamped

    idx = minloc(abs(log(x0) - log_e), 1)

    half_width = max(n, 4)
    half_width = min(half_width, num_tab - 1)

    left = idx - half_width
    right = idx + half_width

    if (left < 1) then
      right = min(num_tab, right + (1 - left))
      left = 1
    end if
    if (right > num_tab) then
      left = max(1, left - (right - num_tab))
      right = num_tab
    end if

    n_points = right - left + 1

    if (n_points <= n) then
      info = 4
      fatal_error = .true.
      if (present(status)) status = info
      derivative = 0.d0
      return
    end if
    
    allocate(nodes(n_points), values(n_points), coeffs(n_points, n+1))

    do i = 1, n_points
      nodes(i) = exp(log_e(left + i - 1))
      values(i) = exp(log_p(left + i - 1))
    end do

    call fornberg_weights(x0, nodes, n_points, n, coeffs)

    do i = 1, n_points
      derivative = derivative + coeffs(i, n+1) * values(i)
    end do

    if (abs(ee_clamped - ee) > clamp_tol * max(1.d0, abs(ee_clamped))) &
      info = max(info, 5)

    deallocate(nodes, values, coeffs)

    if (present(status)) status = info
    if (fatal_error) derivative = 0.d0

  contains
    ! helper that constructs finite-difference weights for arbitrary 
    ! derivative order, allowing the outer routine to reuse the existing 
    ! EOS samples without new interpolation logic.
    subroutine fornberg_weights(x0_local, x, m, n_deriv, c)
      implicit none
      integer, intent(in) :: m, n_deriv
      real(wp), intent(in) :: x0_local
      real(wp), intent(in) :: x(m)
      real(wp), intent(out) :: c(m, n_deriv + 1)

      integer :: i, j, k, max_k
      real(wp) :: c1, c2, c3, c4, c5

      c(:, :) = 0.d0
      c(1,1) = 1.d0
      c1 = 1.d0
      c4 = x(1) - x0_local

      do i = 2, m
        c2 = 1.d0
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
            c(i,1) = -c1 * c5 * c(i-1,1) / c2
          end if
          do k = max_k, 1, -1
            c(j, k+1) = (c4 * c(j, k+1) - k * c(j, k)) / c3
          end do
          c(j,1) = c4 * c(j,1) / c3
        end do
        c1 = c2
      end do
    end subroutine fornberg_weights
  end subroutine pressure_derivative_n

end module eos_mod
