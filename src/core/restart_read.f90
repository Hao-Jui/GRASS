module regrid_mod
  use iso_fortran_env, only: output_unit
  use para_mod, only: wp
  use grid_mod, only: make_grid, GridTrig
  implicit none
  private

  public :: regrid_read

  integer, parameter :: restart_value_count = 13
  character(len=*), parameter :: restart_file_path = "./Res/res.dat"

  integer, parameter :: NFIELDS = 10
  integer, parameter :: F_ALPHA = 1, F_GAMA = 2, F_RHO = 3, F_WW = 4, F_PRESSURE = 5
  integer, parameter :: F_ENERGY = 6, F_ENTHALPY = 7, F_VELOCITY_SQ = 8, F_OMG = 9, F_SPHI = 10

  type :: restart_meta_t
    integer :: sdiv = 0
    integer :: mdiv = 0
    integer :: spwr = 0
    real(wp) :: r_e = 0._wp
    real(wp) :: e_center = 0._wp
    real(wp) :: r_ratio = 0._wp
    real(wp) :: omega_e = 0._wp
    real(wp) :: omega_c = 0._wp
  end type restart_meta_t

contains

  subroutine regrid_read(new_sdiv, new_mdiv, interpolation_order, ierr, errmsg)
    use para_mod, only: s_pwr, s_gp, mu
    implicit none
    integer, intent(in) :: new_sdiv, new_mdiv
    integer, intent(in), optional :: interpolation_order
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg
    type(restart_meta_t) :: old_meta
    real(wp), allocatable :: old_data(:,:,:), new_data(:,:,:)
    real(wp), allocatable :: old_s(:), old_m(:)
    integer :: interp_order, status
    real(wp) :: t0, t1
    character(len=256) :: message

    status = 0
    message = ""
    if (present(ierr)) ierr = 0
    if (present(errmsg)) errmsg = ""

    if (new_sdiv < 2 .or. new_mdiv < 2) then
      status = 1
      message = "regrid_read: target resolution must be >= 2"
    else
      interp_order = 1
      if (present(interpolation_order)) interp_order = interpolation_order
      interp_order = max(1, min(2, interp_order))

      call cpu_time(t0)

      call read_restart_file(restart_file_path, old_meta, old_s, old_m, old_data, status, message)
      if (status == 0) then
        call normalize_restart_meta(old_meta)

        if (old_meta%spwr /= s_pwr) then
          status = 2
          message = "regrid_read: s-grid power doesn't match."
        else
          call ensure_runtime_grid(new_sdiv, new_mdiv)

          if (old_meta%sdiv == new_sdiv .and. old_meta%mdiv == new_mdiv .and. &
              maxval(abs(old_s - s_gp)) <= 32.0_wp * epsilon(1.0_wp) .and. &
              maxval(abs(old_m - mu)) <= 32.0_wp * epsilon(1.0_wp)) then
            call apply_runtime_meta(old_meta)
            call apply_runtime_fields(old_data)
          else
            call regrid_fields(old_data, old_s, old_m, s_gp, mu, interp_order, new_data)
            call apply_runtime_meta(old_meta)
            call apply_runtime_fields(new_data)
          end if

          call cpu_time(t1)
          write(output_unit, fmt=*) " "
          write(output_unit, fmt='(A,f12.6,A)') "Regrid-read ---", t1 - t0, " [s]"
        end if
      end if
    end if

    if (present(ierr)) ierr = status
    if (present(errmsg)) errmsg = trim(message)
  end subroutine regrid_read

  subroutine read_restart_file(path, old_meta, old_s, old_m, old_data, ierr, errmsg)
    use para_mod, only: C, KSCALE
    implicit none
    character(len=*), intent(in) :: path
    type(restart_meta_t), intent(out) :: old_meta
    real(wp), allocatable, intent(out) :: old_s(:), old_m(:)
    real(wp), allocatable, intent(out) :: old_data(:,:,:)
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: errmsg
    integer :: unit, ios, s, m
    character(len=512) :: line
    real(wp) :: vals(restart_value_count)
    integer, parameter :: val_idx(NFIELDS) = [3, 4, 5, 6, 7, 8, 9, 11, 12, 13]
    real(wp) :: scale(NFIELDS)

    ierr = 0
    errmsg = ""

    open(newunit=unit, file=path, status="old", action="read", iostat=ios)
    if (ios /= 0) then
      ierr = 1; errmsg = "regrid_read: failed to open restart file"; return
    end if

    read(unit, '(A)', iostat=ios) line
    if (ios /= 0) then
      ierr = 2; errmsg = "regrid_read: failed to read header"; close(unit); return
    end if

    read(line, *, iostat=ios) old_meta%sdiv, old_meta%mdiv, old_meta%spwr, &
                              old_meta%r_e, old_meta%e_center, old_meta%r_ratio, &
                              old_meta%omega_e, old_meta%omega_c
    if (ios /= 0) then
      ierr = 3; errmsg = "regrid_read: malformed header"; close(unit); return
    end if

    allocate(old_s(old_meta%sdiv), old_m(old_meta%mdiv), source=0._wp)
    allocate(old_data(NFIELDS, old_meta%sdiv, old_meta%mdiv), source=0._wp)

    scale = 1._wp
    scale(F_PRESSURE) = KSCALE
    scale(F_ENERGY) = C * C * KSCALE

    do s = 1, old_meta%sdiv
      do m = 1, old_meta%mdiv
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) then
          ierr = 4; errmsg = "regrid_read: unexpected end of file"; close(unit); return
        end if

        vals = 0._wp
        read(line, *, iostat=ios) vals
        if (ios /= 0) then
          ierr = 1; errmsg = "restart_read: malformed data"; close(unit); return
        end if

        old_s(s) = vals(1)
        old_m(m) = vals(2)
        old_data(:, s, m) = vals(val_idx) * scale
      end do
    end do

    close(unit)
  end subroutine read_restart_file

  subroutine normalize_restart_meta(old_meta)
    use para_mod, only: C, KAPPA
    implicit none
    type(restart_meta_t), intent(inout) :: old_meta

    old_meta%r_e = old_meta%r_e / (sqrt(KAPPA) / 1.e5_wp)
    old_meta%omega_e = old_meta%omega_e / (C / sqrt(KAPPA)) * old_meta%r_e
    old_meta%omega_c = old_meta%omega_c / (C / sqrt(KAPPA)) * old_meta%r_e
  end subroutine normalize_restart_meta

  subroutine apply_runtime_meta(old_meta)
    use para_mod, only: e_center, Omega_e, Omega_c, r_e, r_ratio
    implicit none
    type(restart_meta_t), intent(in) :: old_meta

    r_e = old_meta%r_e
    e_center = old_meta%e_center
    r_ratio = old_meta%r_ratio
    Omega_e = old_meta%omega_e
    Omega_c = old_meta%omega_c
  end subroutine apply_runtime_meta

  subroutine apply_runtime_fields(new_data)
    use para_mod, only: B_goal, C, KAPPA, alpha, gama, rho, ww, pressure, energy, &
                        enthalpy, velocity_sq, omg, sphi, has_scalar
    implicit none
    real(wp), intent(in) :: new_data(:,:,:)
    real(wp) :: angular_scale

    angular_scale = C / sqrt(KAPPA)

    alpha = new_data(F_ALPHA,:,:)
    gama = new_data(F_GAMA,:,:)
    rho = new_data(F_RHO,:,:)
    ww = new_data(F_WW,:,:) / angular_scale
    pressure = new_data(F_PRESSURE,:,:)
    energy = new_data(F_ENERGY,:,:)
    enthalpy = new_data(F_ENTHALPY,:,:)
    velocity_sq = new_data(F_VELOCITY_SQ,:,:)
    omg = new_data(F_OMG,:,:) / angular_scale

    if (has_scalar) then
      sphi = new_data(F_SPHI,:,:) / sqrt(max(B_goal, 1.e-30_wp))
    else
      sphi = 0._wp
      has_scalar = .false.
    end if
  end subroutine apply_runtime_fields

  subroutine ensure_runtime_grid(new_sdiv, new_mdiv)
    use para_mod, only: MDIV, SDIV, SMAX, DS, DM, allocate_fields, &
                        alpha, gama, s_gp, mu, sin_theta, P_2n, P1_2n_1, sin_2n_1_theta
    implicit none
    integer, intent(in) :: new_sdiv, new_mdiv
    logical :: need_reinit, need_grid_rebuild

    need_reinit = (.not. allocated(alpha)) .or. (.not. allocated(gama)) .or. &
                  size(alpha, 1) /= new_sdiv .or. size(alpha, 2) /= new_mdiv

    need_grid_rebuild = (.not. allocated(s_gp)) .or. (.not. allocated(mu)) .or. &
                        (.not. allocated(sin_theta)) .or. (.not. allocated(P_2n)) .or. &
                        (.not. allocated(P1_2n_1)) .or. (.not. allocated(sin_2n_1_theta)) .or. &
                        size(s_gp) /= new_sdiv .or. size(mu) /= new_mdiv .or. &
                        size(P_2n, 1) /= new_mdiv

    if (need_reinit .or. SDIV /= new_sdiv .or. MDIV /= new_mdiv) then
      SDIV = new_sdiv
      MDIV = new_mdiv
      DS = SMAX / (real(SDIV, wp) - 1._wp)
      DM = 1._wp / (real(MDIV, wp) - 1._wp)
      call allocate_fields()
      need_grid_rebuild = .true.
    end if

    if (need_grid_rebuild) then
      call make_grid
      call GridTrig
    end if
  end subroutine ensure_runtime_grid

  subroutine regrid_fields(old_data, old_s, old_m, new_s, new_m, interp_order, new_data)
    implicit none
    real(wp), intent(in) :: old_data(:,:,:)
    real(wp), intent(in) :: old_s(:), old_m(:), new_s(:), new_m(:)
    integer, intent(in) :: interp_order
    real(wp), allocatable, intent(out) :: new_data(:,:,:)
    integer :: s, m, i, j, i0, i1, j0, j1
    real(wp) :: c00, c10, c01, c11, w_s, w_m, wij
    real(wp) :: acc(NFIELDS)
    logical :: use_quadratic
    integer :: old_sdiv, old_mdiv, new_sdiv, new_mdiv
    integer, allocatable :: s_lo(:), s_hi(:), m_lo(:), m_hi(:)
    real(wp), allocatable :: s_wt(:), m_wt(:)
    integer, allocatable :: sq_idx(:,:), mq_idx(:,:)
    real(wp), allocatable :: sq_wt(:,:), mq_wt(:,:)

    old_sdiv = size(old_s)
    old_mdiv = size(old_m)
    new_sdiv = size(new_s)
    new_mdiv = size(new_m)

    allocate(new_data(NFIELDS, new_sdiv, new_mdiv), source=0._wp)

    use_quadratic = (interp_order == 2) .and. old_sdiv >= 3 .and. old_mdiv >= 3

    allocate(s_lo(new_sdiv), s_hi(new_sdiv), s_wt(new_sdiv))
    allocate(m_lo(new_mdiv), m_hi(new_mdiv), m_wt(new_mdiv))

    if (use_quadratic) then
      allocate(sq_idx(3, new_sdiv), sq_wt(3, new_sdiv))
      allocate(mq_idx(3, new_mdiv), mq_wt(3, new_mdiv))
    end if

    do s = 1, new_sdiv
      call locate(old_s, old_sdiv, new_s(s), s_lo(s), s_hi(s), s_wt(s))
      if (use_quadratic) call quadratic_weights(old_s, old_sdiv, new_s(s), sq_idx(:,s), sq_wt(:,s))
    end do

    do m = 1, new_mdiv
      call locate(old_m, old_mdiv, new_m(m), m_lo(m), m_hi(m), m_wt(m))
      if (use_quadratic) call quadratic_weights(old_m, old_mdiv, new_m(m), mq_idx(:,m), mq_wt(:,m))
    end do

    if (use_quadratic) then
      do m = 1, new_mdiv
        do s = 1, new_sdiv
          acc = 0._wp
          do j = 1, 3
            do i = 1, 3
              wij = sq_wt(i, s) * mq_wt(j, m)
              acc = acc + wij * old_data(:, sq_idx(i,s), mq_idx(j,m))
            end do
          end do
          new_data(:, s, m) = acc
        end do
      end do
    else
      do m = 1, new_mdiv
        j0 = m_lo(m); j1 = m_hi(m)
        w_m = merge(0._wp, m_wt(m), j0 == j1)
        do s = 1, new_sdiv
          i0 = s_lo(s); i1 = s_hi(s)
          w_s = merge(0._wp, s_wt(s), i0 == i1)
          c00 = (1._wp - w_s) * (1._wp - w_m)
          c10 = w_s * (1._wp - w_m)
          c01 = (1._wp - w_s) * w_m
          c11 = w_s * w_m
          new_data(:, s, m) = c00 * old_data(:, i0, j0) + c10 * old_data(:, i1, j0) + &
                              c01 * old_data(:, i0, j1) + c11 * old_data(:, i1, j1)
        end do
      end do
    end if
  end subroutine regrid_fields

  subroutine locate(grid, n, value, idx_low, idx_high, weight)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: grid(n), value
    integer, intent(out) :: idx_low, idx_high
    real(wp), intent(out) :: weight
    integer :: lo, hi, mid

    if (n < 2) then
      idx_low = 1; idx_high = 1; weight = 0._wp; return
    end if

    if (value <= grid(1)) then
      idx_low = 1; idx_high = 1; weight = 0._wp
    else if (value >= grid(n)) then
      idx_low = n; idx_high = n; weight = 0._wp
    else
      lo = 1
      hi = n
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (grid(mid) <= value) then
          lo = mid
        else
          hi = mid
        end if
      end do
      idx_low = lo
      idx_high = hi
      weight = (value - grid(lo)) / max(grid(hi) - grid(lo), tiny(1.0_wp))
    end if
  end subroutine locate

  subroutine quadratic_weights(grid, n, value, idx, weights)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: grid(n), value
    integer, intent(out) :: idx(3)
    real(wp), intent(out) :: weights(3)
    integer :: k, i, j
    real(wp) :: x(3), denom, sum_w

    if (n < 3) then
      idx(1) = 1; idx(2) = min(2, n); idx(3) = max(idx(2), 1)
      weights = 0._wp; weights(1) = 1._wp
      if (idx(2) > 1) then
        weights(1) = (grid(idx(2)) - value) / (grid(idx(2)) - grid(idx(1)))
        weights(2) = 1._wp - weights(1)
      end if
      return
    end if

    if (value <= grid(2)) then
      idx = [1, 2, 3]
    else if (value >= grid(n - 1)) then
      idx = [n - 2, n - 1, n]
    else
      idx = 0
      do k = 2, n - 2
        if (value <= grid(k + 1)) then
          idx = [k - 1, k, k + 1]; exit
        end if
      end do
      if (idx(1) == 0) idx = [n - 2, n - 1, n]
    end if

    x = grid(idx)

    do i = 1, 3
      weights(i) = 1._wp
      do j = 1, 3
        if (j == i) cycle
        denom = x(i) - x(j)
        if (abs(denom) <= 1.e-14_wp) then
          weights = 0._wp; weights(i) = 1._wp; return
        end if
        weights(i) = weights(i) * (value - x(j)) / denom
      end do
    end do

    sum_w = sum(weights)
    if (abs(sum_w) > 0._wp) weights = weights / sum_w
  end subroutine quadratic_weights

end module regrid_mod
