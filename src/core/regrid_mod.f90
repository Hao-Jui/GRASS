module regrid_mod
  use iso_fortran_env, only: output_unit, int32
  use para_mod, only: wp
  use grid_mod, only: make_grid, GridTrig
  implicit none
  private
  public :: regrid_read

  character(len=*), parameter :: restart_binary_path = "./Res/res.rst"
  character(len=8), parameter :: restart_magic = "GRASSRST01"
  integer(int32), parameter :: restart_format_version = 1_int32

  integer, parameter :: NFIELDS = 10
  integer, parameter :: F_ALPHA = 1, F_GAMA = 2, F_RHO = 3, F_WW = 4, F_PRESSURE = 5
  integer, parameter :: F_ENERGY = 6, F_ENTHALPY = 7, F_VELOCITY_SQ = 8, F_OMG = 9, F_SPHI = 10

  type :: restart_meta_t
    integer :: sdiv = 0, mdiv = 0, spwr = 0
    real(wp) :: r_e = 0._wp, e_center = 0._wp, r_ratio = 0._wp
    real(wp) :: omega_e = 0._wp, omega_c = 0._wp
  end type restart_meta_t

contains

  subroutine regrid_read(new_sdiv, new_mdiv, interpolation_order, ierr, errmsg)
    use para_mod, only: s_pwr, s_gp, mu, e_center, Omega_e, Omega_c, r_e, r_ratio, &
                        alpha, gama, rho, ww, pressure, energy, &
                        enthalpy, velocity_sq, omg, sphi, has_scalar
    integer, intent(in) :: new_sdiv, new_mdiv
    integer, intent(in), optional :: interpolation_order
    integer, intent(out), optional :: ierr
    character(len=*), intent(out), optional :: errmsg
    type(restart_meta_t) :: meta
    real(wp), allocatable :: old_data(:,:,:), new_data(:,:,:), old_s(:), old_m(:)
    real(wp) :: t0, t1
    integer :: interp_order, status
    character(len=256) :: message

    status = 0; message = ""

    if (new_sdiv < 2 .or. new_mdiv < 2) then
      status = 1; message = "regrid_read: target resolution must be >= 2"
      goto 99
    end if

    interp_order = 1
    if (present(interpolation_order)) interp_order = max(1, min(2, interpolation_order))

    call cpu_time(t0)

    call read_binary_restart(restart_binary_path, meta, old_s, old_m, old_data, status, message)
    if (status /= 0) goto 99

    if (meta%spwr /= s_pwr) then
      status = 2; message = "regrid_read: s-grid power doesn't match."
      goto 99
    end if

    call ensure_runtime_grid(new_sdiv, new_mdiv)

    if (meta%sdiv == new_sdiv .and. meta%mdiv == new_mdiv .and. &
        maxval(abs(old_s - s_gp)) <= 32._wp * epsilon(1._wp) .and. &
        maxval(abs(old_m - mu))   <= 32._wp * epsilon(1._wp)) then
      new_data = old_data
    else
      call regrid_fields(old_data, old_s, old_m, s_gp, mu, interp_order, new_data)
    end if

    r_e = meta%r_e; e_center = meta%e_center; r_ratio = meta%r_ratio
    Omega_e = meta%omega_e; Omega_c = meta%omega_c

    alpha       = new_data(F_ALPHA,:,:)
    gama        = new_data(F_GAMA,:,:)
    rho         = new_data(F_RHO,:,:)
    ww          = new_data(F_WW,:,:)
    pressure    = new_data(F_PRESSURE,:,:)
    energy      = new_data(F_ENERGY,:,:)
    enthalpy    = new_data(F_ENTHALPY,:,:)
    velocity_sq = new_data(F_VELOCITY_SQ,:,:)
    omg         = new_data(F_OMG,:,:)
    if (has_scalar) then
      sphi = new_data(F_SPHI,:,:)
    else
      sphi = 0._wp
    end if

    call cpu_time(t1)
    write(output_unit, '(A)') " "
    write(output_unit, '(A,f12.6,A)') "Regrid-read ---", t1 - t0, " [s]"

    99 continue
    if (present(ierr)) ierr = status
    if (present(errmsg)) errmsg = trim(message)
  end subroutine regrid_read

  subroutine read_binary_restart(path, meta, old_s, old_m, old_data, ierr, errmsg)
    use para_mod, only: C, KSCALE
    character(len=*), intent(in) :: path
    type(restart_meta_t), intent(out) :: meta
    real(wp), allocatable, intent(out) :: old_s(:), old_m(:), old_data(:,:,:)
    integer, intent(out) :: ierr
    character(len=*), intent(out) :: errmsg
    integer :: unit, ios
    integer(int32) :: hi(6)
    real(wp) :: hm(5)
    character(len=len(restart_magic)) :: magic

    ierr = 0; errmsg = ""
    open(newunit=unit, file=path, status="old", action="read", access="stream", &
         form="unformatted", iostat=ios)
    if (ios /= 0) then
      ierr = 1; errmsg = "regrid_read: failed to open " // trim(path); return
    end if

    read(unit, iostat=ios) magic
    if (ios /= 0 .or. magic /= restart_magic) goto 90
    read(unit, iostat=ios) hi
    if (ios /= 0) goto 90
    if (hi(1) /= restart_format_version .or. hi(2) /= int(storage_size(1.0_wp), int32) .or. &
        hi(3) /= int(NFIELDS, int32) .or. hi(4) < 2 .or. hi(5) < 2) goto 90
    read(unit, iostat=ios) hm
    if (ios /= 0) goto 90

    meta = restart_meta_t(sdiv=hi(4), mdiv=hi(5), spwr=hi(6), &
                          r_e=hm(1), e_center=hm(2) / (C * C * KSCALE), r_ratio=hm(3), &
                          omega_e=hm(4) * hm(1), omega_c=hm(5) * hm(1))
    allocate(old_s(meta%sdiv), old_m(meta%mdiv), old_data(NFIELDS, meta%sdiv, meta%mdiv))

    read(unit, iostat=ios) old_s
    if (ios == 0) read(unit, iostat=ios) old_m
    if (ios == 0) read(unit, iostat=ios) old_data
    if (ios /= 0) goto 90

    close(unit); return

    90 ierr = 3; errmsg = "regrid_read: invalid binary restart file"; close(unit)
  end subroutine read_binary_restart

  subroutine ensure_runtime_grid(new_sdiv, new_mdiv)
    use para_mod, only: MDIV, SDIV, SMAX, DS, DM, allocate_fields, &
                        alpha, s_gp, mu, sin_theta, P_2n, P1_2n_1, sin_2n_1_theta
    integer, intent(in) :: new_sdiv, new_mdiv

    if (SDIV /= new_sdiv .or. MDIV /= new_mdiv .or. &
        .not. allocated(alpha) .or. size(alpha, 1) /= new_sdiv) then
      SDIV = new_sdiv; MDIV = new_mdiv
      DS = SMAX / (real(SDIV, wp) - 1._wp)
      DM = 1._wp / (real(MDIV, wp) - 1._wp)
      call allocate_fields()
    end if

    if (.not. allocated(s_gp) .or. size(s_gp) /= new_sdiv .or. &
        .not. allocated(mu) .or. size(mu) /= new_mdiv .or. &
        .not. allocated(P_2n) .or. size(P_2n, 1) /= new_mdiv) then
      call make_grid; call GridTrig
    end if
  end subroutine ensure_runtime_grid

  subroutine regrid_fields(old_data, old_s, old_m, new_s, new_m, interp_order, new_data)
    real(wp), intent(in) :: old_data(:,:,:), old_s(:), old_m(:), new_s(:), new_m(:)
    integer, intent(in) :: interp_order
    real(wp), allocatable, intent(out) :: new_data(:,:,:)
    integer :: s, m, i, j, i0, i1, j0, j1, new_sdiv, new_mdiv
    real(wp) :: w_s, w_m, wij, c00, c10, c01, c11
    logical :: use_quad
    integer, allocatable :: s_lo(:), s_hi(:), m_lo(:), m_hi(:), sq_idx(:,:), mq_idx(:,:)
    real(wp), allocatable :: s_wt(:), m_wt(:), sq_wt(:,:), mq_wt(:,:)

    new_sdiv = size(new_s); new_mdiv = size(new_m)
    allocate(new_data(NFIELDS, new_sdiv, new_mdiv), source=0._wp)
    allocate(s_lo(new_sdiv), s_hi(new_sdiv), s_wt(new_sdiv))
    allocate(m_lo(new_mdiv), m_hi(new_mdiv), m_wt(new_mdiv))

    use_quad = (interp_order == 2) .and. size(old_s) >= 3 .and. size(old_m) >= 3
    if (use_quad) then
      allocate(sq_idx(3, new_sdiv), sq_wt(3, new_sdiv))
      allocate(mq_idx(3, new_mdiv), mq_wt(3, new_mdiv))
    end if

    do s = 1, new_sdiv
      call locate(old_s, new_s(s), s_lo(s), s_hi(s), s_wt(s))
      if (use_quad) call quad_weights(old_s, new_s(s), s_lo(s), sq_idx(:,s), sq_wt(:,s))
    end do
    do m = 1, new_mdiv
      call locate(old_m, new_m(m), m_lo(m), m_hi(m), m_wt(m))
      if (use_quad) call quad_weights(old_m, new_m(m), m_lo(m), mq_idx(:,m), mq_wt(:,m))
    end do

    if (use_quad) then
      do m = 1, new_mdiv; do s = 1, new_sdiv
        do j = 1, 3; do i = 1, 3
          new_data(:,s,m) = new_data(:,s,m) + sq_wt(i,s) * mq_wt(j,m) * old_data(:, sq_idx(i,s), mq_idx(j,m))
        end do; end do
      end do; end do
    else
      do m = 1, new_mdiv
        j0 = m_lo(m); j1 = m_hi(m); w_m = merge(0._wp, m_wt(m), j0 == j1)
        do s = 1, new_sdiv
          i0 = s_lo(s); i1 = s_hi(s); w_s = merge(0._wp, s_wt(s), i0 == i1)
          c00 = (1._wp - w_s) * (1._wp - w_m); c10 = w_s * (1._wp - w_m)
          c01 = (1._wp - w_s) * w_m;           c11 = w_s * w_m
          new_data(:,s,m) = c00*old_data(:,i0,j0) + c10*old_data(:,i1,j0) + &
                            c01*old_data(:,i0,j1) + c11*old_data(:,i1,j1)
        end do
      end do
    end if
  end subroutine regrid_fields

  subroutine locate(grid, value, idx_low, idx_high, weight)
    real(wp), intent(in) :: grid(:), value
    integer, intent(out) :: idx_low, idx_high
    real(wp), intent(out) :: weight
    integer :: lo, hi, mid, n

    n = size(grid)
    if (n < 2 .or. value <= grid(1)) then
      idx_low = 1; idx_high = 1; weight = 0._wp
    else if (value >= grid(n)) then
      idx_low = n; idx_high = n; weight = 0._wp
    else
      lo = 1; hi = n
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (grid(mid) <= value) then; lo = mid; else; hi = mid; end if
      end do
      idx_low = lo; idx_high = hi
      weight = (value - grid(lo)) / max(grid(hi) - grid(lo), tiny(1.0_wp))
    end if
  end subroutine locate

  subroutine quad_weights(grid, value, idx_low, idx, wt)
    real(wp), intent(in) :: grid(:), value
    integer, intent(in) :: idx_low
    integer, intent(out) :: idx(3)
    real(wp), intent(out) :: wt(3)
    integer :: n, lo, i, j
    real(wp) :: x(3), d

    n = size(grid)
    if (n < 3) then
      idx = [1, min(2,n), min(2,n)]; wt = 0._wp; wt(1) = 1._wp; return
    end if

    if (value <= grid(2)) then;      idx = [1, 2, 3]
    else if (value >= grid(n-1)) then; idx = [n-2, n-1, n]
    else; lo = min(max(idx_low, 2), n-2); idx = [lo-1, lo, lo+1]
    end if

    x = grid(idx)
    do i = 1, 3
      wt(i) = 1._wp
      do j = 1, 3
        if (j == i) cycle
        d = x(i) - x(j)
        if (abs(d) <= 1.e-14_wp) then; wt = 0._wp; wt(i) = 1._wp; return; end if
        wt(i) = wt(i) * (value - x(j)) / d
      end do
    end do
    wt = wt / sum(wt)
  end subroutine quad_weights

end module regrid_mod
