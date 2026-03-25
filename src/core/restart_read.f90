module regrid_mod
  use para_mod, only: wp
  use grid_mod, only: make_grid, GridTrig
contains

subroutine parse_restart_line(line, vals)
  implicit none
  character(len=*), intent(in) :: line
  real(wp), intent(out) :: vals(13)
  integer :: ios

  vals = 0._wp

  read(line,*,iostat=ios) vals
  if (ios /= 0) stop "restart_read: malformed data"

end subroutine parse_restart_line

subroutine regrid_read(target_sdiv, target_mdiv, interpolation_order)
  use para_mod, only: r_e, e_center, r_ratio, Omega_e, Omega_c, &
                      KAPPA, C, KSCALE, s_pwr, &
                      SDIV, MDIV, DS, DM, SMAX, &
                      alpha, gama, rho, ww, pressure, energy, &
                      enthalpy, velocity_sq, omg, sphi, &
                      s_gp, mu, sin_theta, P_2n, P1_2n_1, sin_2n_1_theta, &
                      B_goal, has_scalar, allocate_fields
  implicit none
  integer, intent(in) :: target_sdiv, target_mdiv
  integer, intent(in), optional :: interpolation_order
  integer :: res_r, res_t, spwr, ios
  integer :: s, m
  logical :: need_reinit, need_grid_rebuild
  character(len=512) :: line
  real(wp) :: vals(13)
  real(wp), allocatable :: alpha_src(:,:), gama_src(:,:), rho_src(:,:), &
                          ww_src(:,:), pressure_src(:,:), energy_src(:,:), &
                          enthalpy_src(:,:), velocity_sq_src(:,:), omg_src(:,:), &
                          sphi_src(:,:)
  real(wp), allocatable :: s_source(:), m_source(:), s_target(:), m_target(:)
  integer, allocatable :: s_low(:), s_high(:), m_low(:), m_high(:)
  real(wp), allocatable :: s_weight(:), m_weight(:)
  integer, allocatable :: s_quad_idx(:,:), m_quad_idx(:,:)
  real(wp), allocatable :: s_quad_weight(:,:), m_quad_weight(:,:)
  real(wp) :: ds_source, dm_source, ds_target, dm_target
  integer :: i0, i1, j0, j1
  real(wp) :: ws, wm, t0, t1
  integer :: interp_order, unit
  logical :: use_quadratic
  if (target_sdiv < 2 .or. target_mdiv < 2) then
    stop "regrid_read: target resolution must be >= 2"
  end if

  interp_order = 1
  if (present(interpolation_order)) then
    interp_order = interpolation_order
  end if
  if (interp_order < 1) interp_order = 1
  call cpu_time(t0)
  open(newunit=unit, file="./Res/res.dat")
  read(unit,'(A)',iostat=ios) line
  if (ios /= 0) stop "regrid_read: failed to read header"
  read(line,*,iostat=ios) res_r, res_t, spwr, r_e, e_center, r_ratio, Omega_e, Omega_c
  if (ios /= 0) stop "regrid_read: malformed header"

  r_e = r_e / (sqrt(KAPPA)/1.e5_wp)
  Omega_e = Omega_e / (C/sqrt(kappa)) * r_e
  Omega_c = Omega_c / (C/sqrt(kappa)) * r_e
  
  if (spwr /= s_pwr) stop "s-grid power doesn't match."

  if (interp_order > 2) interp_order = 2
  use_quadratic = (interp_order == 2) .and. res_r >= 3 .and. res_t >= 3
  if (.not. use_quadratic) interp_order = 1  ! fallback to bilinear when not enough source points

  allocate(alpha_src(res_r,res_t), gama_src(res_r,res_t), rho_src(res_r,res_t))
  allocate(ww_src(res_r,res_t), pressure_src(res_r,res_t), energy_src(res_r,res_t))
  allocate(enthalpy_src(res_r,res_t), velocity_sq_src(res_r,res_t), omg_src(res_r,res_t))
  allocate(sphi_src(res_r,res_t))

  do s = 1, res_r
    do m = 1, res_t
      read(unit,'(A)',iostat=ios) line
      if (ios /= 0) stop "regrid_read: unexpected end of file"
      call parse_restart_line(line, vals)

      alpha_src(s,m)       = vals(3)
      gama_src (s,m)       = vals(4)
      rho_src  (s,m)       = vals(5)
      ww_src   (s,m)       = vals(6)
      pressure_src(s,m)    = vals(7) * KSCALE
      energy_src  (s,m)    = vals(8) * (C*C*KSCALE)
      enthalpy_src(s,m)    = vals(9)
      velocity_sq_src(s,m) = vals(11)
      omg_src     (s,m)    = vals(12)
      sphi_src(s,m) = vals(13)
    end do
  end do
  close(unit)

  need_reinit = (.not. allocated(alpha)) .or. (.not. allocated(gama)) .or. &
                size(alpha,1) /= target_sdiv .or. size(alpha,2) /= target_mdiv
  need_grid_rebuild = (.not. allocated(s_gp)) .or. (.not. allocated(mu)) .or. (.not. allocated(sin_theta)) .or. &
                      (.not. allocated(P_2n)) .or. (.not. allocated(P1_2n_1)) .or. (.not. allocated(sin_2n_1_theta)) .or. &
                      size(s_gp) /= target_sdiv .or. size(mu) /= target_mdiv .or. size(P_2n,1) /= target_mdiv

  if (need_reinit .or. SDIV /= target_sdiv .or. MDIV /= target_mdiv) then
    SDIV = target_sdiv
    MDIV = target_mdiv
    DS   = SMAX / (real(SDIV, wp) - 1._wp)
    DM   = 1._wp  / (real(MDIV, wp) - 1._wp)
    call allocate_fields()
    need_grid_rebuild = .true.
  end if

  if (need_grid_rebuild) then
    call make_grid
    call GridTrig
  end if

  if (res_r == target_sdiv .and. res_t == target_mdiv) then
    alpha       = alpha_src
    gama        = gama_src
    rho         = rho_src
    ww          = ww_src
    pressure    = pressure_src
    energy      = energy_src
    enthalpy    = enthalpy_src
    velocity_sq = velocity_sq_src
    omg         = omg_src
    sphi        = sphi_src
  else
    ds_source = SMAX / (real(res_r, wp) - 1._wp)
    dm_source = 1._wp  / (real(res_t, wp) - 1._wp)
    ds_target = SMAX / (real(target_sdiv, wp) - 1._wp)
    dm_target = 1._wp  / (real(target_mdiv, wp) - 1._wp)

    allocate(s_source(res_r), m_source(res_t))
    allocate(s_target(target_sdiv), m_target(target_mdiv))
    allocate(s_low(target_sdiv), s_high(target_sdiv), s_weight(target_sdiv))
    allocate(m_low(target_mdiv), m_high(target_mdiv), m_weight(target_mdiv))
    if (use_quadratic) then
      allocate(s_quad_idx(3,target_sdiv), s_quad_weight(3,target_sdiv))
      allocate(m_quad_idx(3,target_mdiv), m_quad_weight(3,target_mdiv))
    end if

    do s = 1, res_r
      s_source(s) = (real(s, wp) - 1._wp) * ds_source
    end do
    s_source(res_r) = SMAX

    do m = 1, res_t
      m_source(m) = (real(m, wp) - 1._wp) * dm_source
    end do
    m_source(res_t) = 1._wp

    do s = 1, target_sdiv
      if (s == target_sdiv) then
        s_target(s) = SMAX
      else
        s_target(s) = (real(s, wp) - 1._wp) * ds_target
      end if
      call locate(s_source, res_r, s_target(s), s_low(s), s_high(s), s_weight(s))
      if (use_quadratic) then
        call quadratic_weights(s_source, res_r, s_target(s), s_quad_idx(:,s), s_quad_weight(:,s))
      end if
    end do

    do m = 1, target_mdiv
      if (m == target_mdiv) then
        m_target(m) = 1._wp
      else
        m_target(m) = (real(m, wp) - 1._wp) * dm_target
      end if
      call locate(m_source, res_t, m_target(m), m_low(m), m_high(m), m_weight(m))
      if (use_quadratic) then
        call quadratic_weights(m_source, res_t, m_target(m), m_quad_idx(:,m), m_quad_weight(:,m))
      end if
    end do

    if (use_quadratic) then
      do s = 1, target_sdiv
        do m = 1, target_mdiv
          call biquadratic_all( &
            s_quad_idx(:,s), m_quad_idx(:,m), s_quad_weight(:,s), m_quad_weight(:,m), &
            alpha(s,m), gama(s,m), rho(s,m), ww(s,m), pressure(s,m), &
            energy(s,m), enthalpy(s,m), velocity_sq(s,m), omg(s,m), sphi(s,m))
        end do
      end do
    else
      do s = 1, target_sdiv
        i0 = s_low(s)
        i1 = s_high(s)
        ws = s_weight(s)
        do m = 1, target_mdiv
          j0 = m_low(m)
          j1 = m_high(m)
          wm = m_weight(m)

          call bilinear_all(i0, i1, j0, j1, ws, wm, &
            alpha(s,m), gama(s,m), rho(s,m), ww(s,m), pressure(s,m), &
            energy(s,m), enthalpy(s,m), velocity_sq(s,m), omg(s,m), sphi(s,m))
        end do
      end do
    end if
  end if

  ww  = ww  / (C/sqrt(kappa))
  omg = omg / (C/sqrt(kappa))

  if (has_scalar) then
    sphi = sphi / sqrt(max(B_goal, 1.e-30_wp))
  else
    sphi = 0._wp
    has_scalar = .false.
  end if

  call cpu_time(t1)
  write(*, fmt=*) " "
  write(*, fmt='(A,f12.6,A)') "Regrid-read ---", t1-t0, " [s]"

contains

  subroutine locate(grid, n, value, idx_low, idx_high, weight)
    implicit none
    integer, intent(in) :: n
    real(wp), intent(in) :: grid(n), value
    integer, intent(out) :: idx_low, idx_high
    real(wp), intent(out) :: weight
    real(wp) :: inv_h, t

    if (n < 2) then
      idx_low = 1; idx_high = 1; weight = 0._wp
      return
    end if

    inv_h = real(n - 1, wp) / (grid(n) - grid(1))
    t = (value - grid(1)) * inv_h

    if (t <= 0._wp) then
      idx_low = 1; idx_high = 1; weight = 0._wp
    else if (t >= real(n - 1, wp)) then
      idx_low = n; idx_high = n; weight = 0._wp
    else
      idx_low = 1 + int(t)
      idx_high = idx_low + 1
      weight = t - real(idx_low - 1, wp)
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
      idx(1) = 1
      idx(2) = min(2, n)
      idx(3) = max(idx(2), 1)
      weights = 0._wp
      weights(1) = 1._wp
      if (idx(2) > 1) then
        weights(1) = (grid(idx(2)) - value) / (grid(idx(2)) - grid(idx(1)))
        weights(2) = 1._wp - weights(1)
      end if
      return
    end if

    if (value <= grid(2)) then
      idx = (/1, 2, 3/)
    else if (value >= grid(n-1)) then
      idx = (/n-2, n-1, n/)
    else
      idx = 0
      do k = 2, n - 2
        if (value <= grid(k+1)) then
          idx = (/k-1, k, k+1/)
          exit
        end if
      end do
      if (idx(1) == 0) then
        idx = (/n-2, n-1, n/)
      end if
    end if

    x(1) = grid(idx(1))
    x(2) = grid(idx(2))
    x(3) = grid(idx(3))

    do i = 1, 3
      weights(i) = 1._wp
      do j = 1, 3
        if (j == i) cycle
        denom = x(i) - x(j)
        if (abs(denom) <= 1.e-14_wp) then
          weights = 0._wp
          weights(i) = 1._wp
          return
        end if
        weights(i) = weights(i) * (value - x(j)) / denom
      end do
    end do

    sum_w = weights(1) + weights(2) + weights(3)
    if (abs(sum_w) > 0._wp) then
      weights = weights / sum_w
    end if
  end subroutine quadratic_weights

  real(wp) function biquadratic(field, idx_s, idx_m, ws, wm)
    implicit none
    real(wp), intent(in) :: field(:,:)
    integer, intent(in) :: idx_s(3), idx_m(3)
    real(wp), intent(in) :: ws(3), wm(3)
    real(wp) :: interp_s(3)
    integer :: j, i

    do j = 1, 3
      interp_s(j) = 0._wp
      do i = 1, 3
        interp_s(j) = interp_s(j) + ws(i) * field(idx_s(i), idx_m(j))
      end do
    end do

    biquadratic = wm(1) * interp_s(1) + wm(2) * interp_s(2) + wm(3) * interp_s(3)
  end function biquadratic

  subroutine biquadratic_all(idx_s, idx_m, ws, wm, alpha_out, gama_out, rho_out, ww_out, &
                            pressure_out, energy_out, enthalpy_out, velocity_sq_out, omg_out, sphi_out)
    implicit none
    integer, intent(in) :: idx_s(3), idx_m(3)
    real(wp), intent(in) :: ws(3), wm(3)
    real(wp), intent(out) :: alpha_out, gama_out, rho_out, ww_out, pressure_out
    real(wp), intent(out) :: energy_out, enthalpy_out, velocity_sq_out, omg_out, sphi_out
    real(wp) :: interp_col(3,10)
    integer :: i, j

    interp_col = 0._wp
    do j = 1, 3
      associate(is => idx_s, jm => idx_m(j))
        do i = 1, 3
          interp_col(j,1)  = interp_col(j,1)  + ws(i) * alpha_src(is(i), jm)
          interp_col(j,2)  = interp_col(j,2)  + ws(i) * gama_src(is(i), jm)
          interp_col(j,3)  = interp_col(j,3)  + ws(i) * rho_src(is(i), jm)
          interp_col(j,4)  = interp_col(j,4)  + ws(i) * ww_src(is(i), jm)
          interp_col(j,5)  = interp_col(j,5)  + ws(i) * pressure_src(is(i), jm)
          interp_col(j,6)  = interp_col(j,6)  + ws(i) * energy_src(is(i), jm)
          interp_col(j,7)  = interp_col(j,7)  + ws(i) * enthalpy_src(is(i), jm)
          interp_col(j,8)  = interp_col(j,8)  + ws(i) * velocity_sq_src(is(i), jm)
          interp_col(j,9)  = interp_col(j,9)  + ws(i) * omg_src(is(i), jm)
          interp_col(j,10) = interp_col(j,10) + ws(i) * sphi_src(is(i), jm)
        end do
      end associate
    end do

    alpha_out       = dot_product(wm, interp_col(:,1))
    gama_out        = dot_product(wm, interp_col(:,2))
    rho_out         = dot_product(wm, interp_col(:,3))
    ww_out          = dot_product(wm, interp_col(:,4))
    pressure_out    = dot_product(wm, interp_col(:,5))
    energy_out      = dot_product(wm, interp_col(:,6))
    enthalpy_out    = dot_product(wm, interp_col(:,7))
    velocity_sq_out = dot_product(wm, interp_col(:,8))
    omg_out         = dot_product(wm, interp_col(:,9))
    sphi_out        = dot_product(wm, interp_col(:,10))
  end subroutine biquadratic_all

  real(wp) function bilinear(field, i0, i1, j0, j1, ws, wm)
    implicit none
    real(wp), intent(in) :: field(:,:)
    integer, intent(in) :: i0, i1, j0, j1
    real(wp), intent(in) :: ws, wm
    real(wp) :: w_s, w_m
    real(wp) :: f00, f10, f01, f11

    w_s = merge(0._wp, ws, i0 == i1)
    w_m = merge(0._wp, wm, j0 == j1)

    f00 = field(i0,j0)
    f10 = field(i1,j0)
    f01 = field(i0,j1)
    f11 = field(i1,j1)

    bilinear = f00*(1._wp-w_s)*(1._wp-w_m) + f10*w_s*(1._wp-w_m) + &
              f01*(1._wp-w_s)*w_m + f11*w_s*w_m
  end function bilinear

  subroutine bilinear_all(i0, i1, j0, j1, ws, wm, alpha_out, gama_out, rho_out, ww_out, &
                          pressure_out, energy_out, enthalpy_out, velocity_sq_out, omg_out, sphi_out)
    implicit none
    integer, intent(in) :: i0, i1, j0, j1
    real(wp), intent(in) :: ws, wm
    real(wp), intent(out) :: alpha_out, gama_out, rho_out, ww_out, pressure_out
    real(wp), intent(out) :: energy_out, enthalpy_out, velocity_sq_out, omg_out, sphi_out
    real(wp) :: w_s, w_m
    real(wp) :: c00, c10, c01, c11

    w_s = merge(0._wp, ws, i0 == i1)
    w_m = merge(0._wp, wm, j0 == j1)

    c00 = (1._wp - w_s) * (1._wp - w_m)
    c10 = w_s * (1._wp - w_m)
    c01 = (1._wp - w_s) * w_m
    c11 = w_s * w_m

    alpha_out       = c00*alpha_src(i0,j0)       + c10*alpha_src(i1,j0)       + c01*alpha_src(i0,j1)       + c11*alpha_src(i1,j1)
    gama_out        = c00*gama_src(i0,j0)        + c10*gama_src(i1,j0)        + c01*gama_src(i0,j1)        + c11*gama_src(i1,j1)
    rho_out         = c00*rho_src(i0,j0)         + c10*rho_src(i1,j0)         + c01*rho_src(i0,j1)         + c11*rho_src(i1,j1)
    ww_out          = c00*ww_src(i0,j0)          + c10*ww_src(i1,j0)          + c01*ww_src(i0,j1)          + c11*ww_src(i1,j1)
    pressure_out    = c00*pressure_src(i0,j0)    + c10*pressure_src(i1,j0)    + c01*pressure_src(i0,j1)    + c11*pressure_src(i1,j1)
    energy_out      = c00*energy_src(i0,j0)      + c10*energy_src(i1,j0)      + c01*energy_src(i0,j1)      + c11*energy_src(i1,j1)
    enthalpy_out    = c00*enthalpy_src(i0,j0)    + c10*enthalpy_src(i1,j0)    + c01*enthalpy_src(i0,j1)    + c11*enthalpy_src(i1,j1)
    velocity_sq_out = c00*velocity_sq_src(i0,j0) + c10*velocity_sq_src(i1,j0) + c01*velocity_sq_src(i0,j1) + c11*velocity_sq_src(i1,j1)
    omg_out         = c00*omg_src(i0,j0)         + c10*omg_src(i1,j0)         + c01*omg_src(i0,j1)         + c11*omg_src(i1,j1)
    sphi_out        = c00*sphi_src(i0,j0)        + c10*sphi_src(i1,j0)        + c01*sphi_src(i0,j1)        + c11*sphi_src(i1,j1)
  end subroutine bilinear_all

end subroutine regrid_read

end module regrid_mod
