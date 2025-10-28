subroutine restart_read
  use para_mod
  implicit none
  integer :: res_r, res_t, spwr, ios
  integer :: s, m
  logical :: includes_scalar
  character(len=512) :: line
  real(8) :: vals(13)

  sphi = 0.d0

  open(89, file="./Res/res.dat" )
  read(89,"(3i5,99es27.17)") res_r, res_t, spwr, r_e, e_center, r_ratio, Omega_e, Omega_c
  r_e = r_e / (sqrt(KAPPA)/1.d5)
  Omega_e = Omega_e / (C/sqrt(kappa)) * r_e
  Omega_c = Omega_c / (C/sqrt(kappa)) * r_e
  if (res_r .ne. SDIV) stop "Difference in the resolution."
  if ( spwr .ne. s_pwr ) stop "s-grid power doesn't match."

  includes_scalar = .false.
  has_scalar = .false.

  do s = 1, SDIV
    do m = 1, MDIV
      read(89,'(A)',iostat=ios) line
      if (ios /= 0) stop "restart_read: unexpected end of file"
      call parse_restart_line(line, vals, includes_scalar)

      alpha(s,m)       = vals(3)
      gama (s,m)       = vals(4)
      rho  (s,m)       = vals(5)
      ww   (s,m)       = vals(6)
      pressure(s,m)   = vals(7) * KSCALE
      energy(s,m)     = vals(8) * (C*C*KSCALE)
      enthalpy(s,m)   = vals(9)
      velocity_sq(s,m)= vals(11)
      omg(s,m)        = vals(12)
      if (includes_scalar) then
        sphi(s,m) = vals(13)
      else
        sphi(s,m) = 0.d0
      end if
    end do
  end do

  if (includes_scalar) then
    !sphi = sphi / sqrt(max(B_goal, 1.d-30))
    has_scalar = .true.
  end if

  do s = 1, SDIV
    do m = 1, MDIV
      ww(s,m)  = ww(s,m)  / (C/sqrt(kappa))
      omg(s,m) = omg(s,m) / (C/sqrt(kappa))
    end do
  end do

  close(89)
  write(*, fmt=*) " "
  write(*, fmt=*) "Restart OK!"

end subroutine restart_read

subroutine parse_restart_line(line, vals, includes_scalar)
  implicit none
  character(len=*), intent(in) :: line
  real(8), intent(out) :: vals(13)
  logical, intent(inout) :: includes_scalar
  integer :: ios

  vals = 0.d0

  if (includes_scalar) then
    read(line,*,iostat=ios) vals
    if (ios /= 0) stop "restart_read: inconsistent scalar data"
    return
  end if

  read(line,*,iostat=ios) vals
  if (ios == 0) then
    includes_scalar = .true.
  else
    vals = 0.d0
    read(line,*,iostat=ios) vals(1:12)
    if (ios /= 0) stop "restart_read: malformed data"
    includes_scalar = .false.
    vals(13) = 0.d0
  end if

end subroutine parse_restart_line

subroutine regrid_read(target_sdiv, target_mdiv, source_sdiv, source_mdiv)
  use para_mod
  implicit none
  integer, intent(in) :: target_sdiv, target_mdiv
  integer, intent(out), optional :: source_sdiv, source_mdiv
  integer :: res_r, res_t, spwr, ios
  integer :: s, m
  logical :: includes_scalar
  logical :: need_reinit
  character(len=512) :: line
  real(8) :: vals(13)
  real(8), allocatable :: alpha_src(:,:), gama_src(:,:), rho_src(:,:), &
                          ww_src(:,:), pressure_src(:,:), energy_src(:,:), &
                          enthalpy_src(:,:), velocity_sq_src(:,:), omg_src(:,:), &
                          sphi_src(:,:)
  real(8), allocatable :: s_source(:), m_source(:), s_target(:), m_target(:)
  integer, allocatable :: s_low(:), s_high(:), m_low(:), m_high(:)
  real(8), allocatable :: s_weight(:), m_weight(:)
  real(8) :: ds_source, dm_source, ds_target, dm_target
  integer :: i0, i1, j0, j1
  real(8) :: ws, wm
  external :: make_grid, GridTrig

  if (target_sdiv < 2 .or. target_mdiv < 2) then
    stop "regrid_read: target resolution must be >= 2"
  end if

  open(89, file="./Res/res.dat")
  read(89,"(3i5,99es27.17)") res_r, res_t, spwr, r_e, e_center, r_ratio, Omega_e, Omega_c
  r_e = r_e / (sqrt(KAPPA)/1.d5)
  Omega_e = Omega_e / (C/sqrt(kappa)) * r_e
  Omega_c = Omega_c / (C/sqrt(kappa)) * r_e
  if (spwr /= s_pwr) stop "s-grid power doesn't match."

  includes_scalar = .false.

  allocate(alpha_src(res_r,res_t), gama_src(res_r,res_t), rho_src(res_r,res_t))
  allocate(ww_src(res_r,res_t), pressure_src(res_r,res_t), energy_src(res_r,res_t))
  allocate(enthalpy_src(res_r,res_t), velocity_sq_src(res_r,res_t), omg_src(res_r,res_t))
  allocate(sphi_src(res_r,res_t))

  do s = 1, res_r
    do m = 1, res_t
      read(89,'(A)',iostat=ios) line
      if (ios /= 0) stop "regrid_read: unexpected end of file"
      call parse_restart_line(line, vals, includes_scalar)

      alpha_src(s,m)       = vals(3)
      gama_src (s,m)       = vals(4)
      rho_src  (s,m)       = vals(5)
      ww_src   (s,m)       = vals(6)
      pressure_src(s,m)    = vals(7) * KSCALE
      energy_src  (s,m)    = vals(8) * (C*C*KSCALE)
      enthalpy_src(s,m)    = vals(9)
      velocity_sq_src(s,m) = vals(11)
      omg_src     (s,m)    = vals(12)
      if (includes_scalar) then
        sphi_src(s,m) = vals(13)
      else
        sphi_src(s,m) = 0.d0
      end if
    end do
  end do
  close(89)

  if (present(source_sdiv)) source_sdiv = res_r
  if (present(source_mdiv)) source_mdiv = res_t

  need_reinit = (.not. allocated(alpha)) .or. (.not. allocated(gama)) .or. &
                size(alpha,1) /= target_sdiv .or. size(alpha,2) /= target_mdiv

  if (need_reinit .or. SDIV /= target_sdiv .or. MDIV /= target_mdiv) then
    SDIV = target_sdiv
    MDIV = target_mdiv
    DS   = SMAX / (dble(SDIV) - 1.d0)
    DM   = 1.d0  / (dble(MDIV) - 1.d0)
    call allocate_fields()
    call make_grid
    call GridTrig
  end if

  ds_source = SMAX / (dble(res_r) - 1.d0)
  dm_source = 1.d0  / (dble(res_t) - 1.d0)
  ds_target = SMAX / (dble(target_sdiv) - 1.d0)
  dm_target = 1.d0  / (dble(target_mdiv) - 1.d0)

  allocate(s_source(res_r), m_source(res_t))
  allocate(s_target(target_sdiv), m_target(target_mdiv))
  allocate(s_low(target_sdiv), s_high(target_sdiv), s_weight(target_sdiv))
  allocate(m_low(target_mdiv), m_high(target_mdiv), m_weight(target_mdiv))

  do s = 1, res_r
    s_source(s) = (dble(s) - 1.d0) * ds_source
  end do
  s_source(res_r) = SMAX

  do m = 1, res_t
    m_source(m) = (dble(m) - 1.d0) * dm_source
  end do
  m_source(res_t) = 1.d0

  do s = 1, target_sdiv
    if (s == target_sdiv) then
      s_target(s) = SMAX
    else
      s_target(s) = (dble(s) - 1.d0) * ds_target
    end if
    call locate(s_source, res_r, s_target(s), s_low(s), s_high(s), s_weight(s))
  end do

  do m = 1, target_mdiv
    if (m == target_mdiv) then
      m_target(m) = 1.d0
    else
      m_target(m) = (dble(m) - 1.d0) * dm_target
    end if
    call locate(m_source, res_t, m_target(m), m_low(m), m_high(m), m_weight(m))
  end do

  do s = 1, target_sdiv
    i0 = s_low(s)
    i1 = s_high(s)
    ws = s_weight(s)
    do m = 1, target_mdiv
      j0 = m_low(m)
      j1 = m_high(m)
      wm = m_weight(m)

      alpha(s,m)       = bilinear(alpha_src,       i0, i1, j0, j1, ws, wm)
      gama (s,m)       = bilinear(gama_src,        i0, i1, j0, j1, ws, wm)
      rho  (s,m)       = bilinear(rho_src,         i0, i1, j0, j1, ws, wm)
      ww   (s,m)       = bilinear(ww_src,          i0, i1, j0, j1, ws, wm)
      pressure(s,m)    = bilinear(pressure_src,    i0, i1, j0, j1, ws, wm)
      energy(s,m)      = bilinear(energy_src,      i0, i1, j0, j1, ws, wm)
      enthalpy(s,m)    = bilinear(enthalpy_src,    i0, i1, j0, j1, ws, wm)
      velocity_sq(s,m) = bilinear(velocity_sq_src, i0, i1, j0, j1, ws, wm)
      omg(s,m)         = bilinear(omg_src,         i0, i1, j0, j1, ws, wm)
      if (includes_scalar) then
        sphi(s,m) = bilinear(sphi_src, i0, i1, j0, j1, ws, wm)
      end if
    end do
  end do

  ww  = ww  / (C/sqrt(kappa))
  omg = omg / (C/sqrt(kappa))

  if (includes_scalar) then
    !sphi = sphi / sqrt(max(B_goal, 1.d-30))
    has_scalar = .true.
  else
    sphi = 0.d0
    has_scalar = .false.
  end if

  write(*, fmt=*) " "
  write(*, fmt=*) "Regrid-read OK!"

contains

  subroutine locate(grid, n, value, idx_low, idx_high, weight)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: grid(n), value
    integer, intent(out) :: idx_low, idx_high
    real(8), intent(out) :: weight
    integer :: k

    if (value <= grid(1)) then
      idx_low = 1
      idx_high = 1
      weight = 0.d0
      return
    end if

    if (value >= grid(n)) then
      idx_low = n
      idx_high = n
      weight = 0.d0
      return
    end if

    do k = 1, n - 1
      if (value <= grid(k+1)) then
        idx_low = k
        idx_high = k + 1
        weight = (value - grid(k)) / (grid(k+1) - grid(k))
        return
      end if
    end do

    idx_low = n - 1
    idx_high = n
    weight = 1.d0
  end subroutine locate

  real(8) function bilinear(field, i0, i1, j0, j1, ws, wm)
    implicit none
    real(8), intent(in) :: field(:,:)
    integer, intent(in) :: i0, i1, j0, j1
    real(8), intent(in) :: ws, wm
    real(8) :: w_s, w_m
    real(8) :: f00, f10, f01, f11

    w_s = merge(0.d0, ws, i0 == i1)
    w_m = merge(0.d0, wm, j0 == j1)

    f00 = field(i0,j0)
    f10 = field(i1,j0)
    f01 = field(i0,j1)
    f11 = field(i1,j1)

    bilinear = f00*(1.d0-w_s)*(1.d0-w_m) + f10*w_s*(1.d0-w_m) + &
               f01*(1.d0-w_s)*w_m + f11*w_s*w_m
  end function bilinear

end subroutine regrid_read
