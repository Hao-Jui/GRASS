module donu_mod
  use precision_mod, only: wp
  use para_mod, only: r_e
  implicit none
  private
  public :: donutization_number

  integer, parameter :: default_eta_size = 600
  real(wp), parameter :: ONE_THIRD = 1.0_wp / 3.0_wp

contains

  function donutization_number(phi_jordan, volume_density, radial_weights, angular_weights, eta, eps) result(D)
    real(wp), intent(in) :: phi_jordan(:,:), volume_density(:,:), radial_weights(:), angular_weights(:)
    real(wp), intent(in), optional :: eta(:)
    real(wp), intent(in), optional :: eps
    real(wp) :: D
    real(wp), allocatable :: phi_sq(:,:), u(:,:), eta_work(:), h(:), inner_radius(:), outer_radius(:)
    real(wp), allocatable :: cell_volume(:,:)
    logical, allocatable :: active_mask(:,:), hole_mask(:,:)
    integer, allocatable :: q_i(:), q_j(:), peak_idx(:)
    real(wp) :: eps_use, umax, active_vol, hole_vol, outer_vol, total_volume, active_volume
    real(wp) :: weight_norm, f_hole, f_extent
    integer :: n_s, n_m, n_cells, n_eta, j
    
    n_s = size(phi_jordan, 1)
    n_m = size(phi_jordan, 2)
    if (n_s < 3 .or. n_m < 2) error stop "donutization_number_axisym: need at least a 3x2 grid"
    if (size(volume_density, 1) /= n_s .or. size(volume_density, 2) /= n_m) &
      error stop "donutization_number_axisym: volume_density size mismatch"
    if (size(radial_weights) /= n_s) error stop "donutization_number_axisym: radial_weights size mismatch"
    if (size(angular_weights) /= n_m) error stop "donutization_number_axisym: angular_weights size mismatch"
    if (any(radial_weights < 0.0_wp) .or. any(angular_weights < 0.0_wp) .or. any(volume_density < 0.0_wp)) &
      error stop "donutization_number_axisym: weights and density must be non-negative"

    if (present(eps)) then
      eps_use = eps
    else
      eps_use = 1.0e-5_wp
    end if
    
    D = 0.0_wp
    allocate(phi_sq(n_s, n_m), u(n_s, n_m))
    phi_sq = abs(phi_jordan)**2
    u = phi_sq
    allocate(peak_idx(2))
    peak_idx = maxloc(u)
    umax = u(peak_idx(1), peak_idx(2))
    if (umax < eps_use) then
      deallocate(phi_sq, u, peak_idx)
      return
    end if
    u = u / umax 

    call build_eta_grid(eta, eta_work)
    n_eta = size(eta_work)
    allocate(h(n_eta), source=0.0_wp)
    allocate(inner_radius(n_eta), outer_radius(n_eta), source=0.0_wp)
    allocate(cell_volume(n_s, n_m))

    cell_volume = volume_density * spread(radial_weights, 2, n_m) * spread(angular_weights, 1, n_s)
    total_volume = sum(cell_volume)
    if (total_volume <= 0.0_wp) then
      deallocate(phi_sq, u, eta_work, h, inner_radius, outer_radius, cell_volume, peak_idx)
      return
    end if

    active_volume = sum(cell_volume, mask=(u > eps_use))
    if (active_volume > 0.0_wp) then
      f_extent = sum(sqrt(u) * cell_volume, mask=(u > eps_use)) !/ total_volume
    else
      f_extent = 0.0_wp
    end if

    n_cells = n_s * n_m
    allocate(active_mask(n_s, n_m), hole_mask(n_s, n_m))
    allocate(q_i(n_cells), q_j(n_cells))
    do j = 1, n_eta
      call level_volumes_axisym(u, cell_volume, peak_idx(1), peak_idx(2), eta_work(j), &
                                active_mask, hole_mask, q_i, q_j, active_vol, hole_vol)
      outer_vol = active_vol + hole_vol
      if (outer_vol > 0.0_wp) then
        outer_radius(j) = outer_vol**one_third
        inner_radius(j) = max(hole_vol, 0.0_wp)**one_third
      end if
    end do

    do j = 1, n_eta
      if (inner_radius(j) > 0.0_wp) h(j) = 1.0_wp
    end do

    weight_norm = trapezoid(eta_work, eta_work)
    if (weight_norm > 0.0_wp) then
      f_hole = trapezoid(eta_work, eta_work * h) / weight_norm
    else
      f_hole = 0.0_wp
    end if
    D = f_hole !1.0_wp - (1.0_wp - f_hole) * (1.0_wp - f_extent)
    !write(*,'("f_hole:", f12.6, "   f_extent", f12.6)') f_hole, f_extent

    deallocate(phi_sq, u, eta_work, h, inner_radius, outer_radius, cell_volume, &
               active_mask, hole_mask, q_i, q_j, peak_idx)
  end function donutization_number

  subroutine build_eta_grid(eta_in, eta_out)
    real(wp), intent(in), optional :: eta_in(:)
    real(wp), allocatable, intent(out) :: eta_out(:)
    real(wp), allocatable :: tmp(:)
    integer :: i, n_valid

    if (.not. present(eta_in)) then
      allocate(eta_out(default_eta_size))
      do i = 1, default_eta_size
        eta_out(i) = 1.0e-3_wp + (1.0_wp - 1.0e-3_wp) * real(i - 1, wp) / real(default_eta_size - 1, wp)
      end do
      return
    end if

    n_valid = count(eta_in > 0.0_wp .and. eta_in <= 1.0_wp)
    if (n_valid == 0) error stop "donutization_number: eta must contain values in (0,1]"

    allocate(tmp(n_valid))
    n_valid = 0
    do i = 1, size(eta_in)
      if (eta_in(i) > 0.0_wp .and. eta_in(i) <= 1.0_wp) then
        n_valid = n_valid + 1
        tmp(n_valid) = eta_in(i)
      end if
    end do
    call sort_increasing(tmp)
    allocate(eta_out(size(tmp)))
    eta_out = tmp
    deallocate(tmp)
  end subroutine build_eta_grid

  subroutine level_volumes_axisym(u, cell_volume, peak_s, peak_m, level, active_mask, hole_mask, q_i, q_j, active_vol, hole_vol)
    real(wp), intent(in) :: u(:,:), level
    real(wp), intent(in) :: cell_volume(:,:)
    integer, intent(in) :: peak_s, peak_m
    logical, intent(out) :: active_mask(:,:), hole_mask(:,:)
    integer, intent(inout) :: q_i(:), q_j(:)
    real(wp), intent(out) :: active_vol, hole_vol
    integer :: m

    active_mask = .false.
    hole_mask = .false.
    active_vol = 0.0_wp
    hole_vol = 0.0_wp
    if (u(peak_s, peak_m) < level) return

    call flood_fill_active(u, cell_volume, level, peak_s, peak_m, active_mask, q_i, q_j, active_vol)
    do m = 1, size(u, 2)
      if (.not. active_mask(1, m)) call flood_fill_hole(cell_volume, active_mask, 1, m, hole_mask, q_i, q_j, hole_vol)
    end do
  end subroutine level_volumes_axisym

  subroutine flood_fill_active(u, cell_volume, level, seed_s, seed_m, active_mask, q_i, q_j, volume_sum)
    real(wp), intent(in) :: u(:,:), cell_volume(:,:), level
    integer, intent(in) :: seed_s, seed_m
    logical, intent(inout) :: active_mask(:,:)
    integer, intent(inout) :: q_i(:), q_j(:)
    real(wp), intent(out) :: volume_sum
    integer :: head, tail, i, j

    head = 1
    tail = 1
    q_i(1) = seed_s
    q_j(1) = seed_m
    active_mask(seed_s, seed_m) = .true.
    volume_sum = cell_volume(seed_s, seed_m)

    do while (head <= tail)
      i = q_i(head)
      j = q_j(head)
      head = head + 1
      call try_push_active(i - 1, j, u, cell_volume, level, active_mask, q_i, q_j, tail, volume_sum)
      call try_push_active(i + 1, j, u, cell_volume, level, active_mask, q_i, q_j, tail, volume_sum)
      call try_push_active(i, j - 1, u, cell_volume, level, active_mask, q_i, q_j, tail, volume_sum)
      call try_push_active(i, j + 1, u, cell_volume, level, active_mask, q_i, q_j, tail, volume_sum)
    end do
  end subroutine flood_fill_active

  subroutine flood_fill_hole(cell_volume, active_mask, seed_s, seed_m, hole_mask, q_i, q_j, volume_sum)
    real(wp), intent(in) :: cell_volume(:,:)
    logical, intent(in) :: active_mask(:,:)
    integer, intent(in) :: seed_s, seed_m
    logical, intent(inout) :: hole_mask(:,:)
    integer, intent(inout) :: q_i(:), q_j(:)
    real(wp), intent(inout) :: volume_sum
    integer :: head, tail, i, j

    if (active_mask(seed_s, seed_m) .or. hole_mask(seed_s, seed_m)) return

    head = 1
    tail = 1
    q_i(1) = seed_s
    q_j(1) = seed_m
    hole_mask(seed_s, seed_m) = .true.
    volume_sum = volume_sum + cell_volume(seed_s, seed_m)

    do while (head <= tail)
      i = q_i(head)
      j = q_j(head)
      head = head + 1
      call try_push_hole(i - 1, j, cell_volume, active_mask, hole_mask, q_i, q_j, tail, volume_sum)
      call try_push_hole(i + 1, j, cell_volume, active_mask, hole_mask, q_i, q_j, tail, volume_sum)
      call try_push_hole(i, j - 1, cell_volume, active_mask, hole_mask, q_i, q_j, tail, volume_sum)
      call try_push_hole(i, j + 1, cell_volume, active_mask, hole_mask, q_i, q_j, tail, volume_sum)
    end do
  end subroutine flood_fill_hole

  subroutine try_push_active(i, j, u, cell_volume, level, active_mask, q_i, q_j, tail, volume_sum)
    integer, intent(in) :: i, j
    real(wp), intent(in) :: u(:,:), cell_volume(:,:), level
    logical, intent(inout) :: active_mask(:,:)
    integer, intent(inout) :: q_i(:), q_j(:), tail
    real(wp), intent(inout) :: volume_sum

    if (i < 1 .or. i > size(u, 1) .or. j < 1 .or. j > size(u, 2)) return
    if (active_mask(i, j)) return
    if (u(i, j) < level) return
    tail = tail + 1
    q_i(tail) = i
    q_j(tail) = j
    active_mask(i, j) = .true.
    volume_sum = volume_sum + cell_volume(i, j)
  end subroutine try_push_active

  subroutine try_push_hole(i, j, cell_volume, active_mask, hole_mask, q_i, q_j, tail, volume_sum)
    integer, intent(in) :: i, j
    real(wp), intent(in) :: cell_volume(:,:)
    logical, intent(in) :: active_mask(:,:)
    logical, intent(inout) :: hole_mask(:,:)
    integer, intent(inout) :: q_i(:), q_j(:), tail
    real(wp), intent(inout) :: volume_sum

    if (i < 1 .or. i > size(active_mask, 1) .or. j < 1 .or. j > size(active_mask, 2)) return
    if (active_mask(i, j) .or. hole_mask(i, j)) return
    tail = tail + 1
    q_i(tail) = i
    q_j(tail) = j
    hole_mask(i, j) = .true.
    volume_sum = volume_sum + cell_volume(i, j)
  end subroutine try_push_hole

  pure function trapezoid(x, y) result(intval)
    real(wp), intent(in) :: x(:), y(:)
    real(wp) :: intval
    integer :: i

    if (size(x) /= size(y)) error stop "trapezoid: size mismatch"
    if (size(x) < 2) then
      intval = 0.0_wp
      return
    end if

    intval = 0.0_wp
    do i = 1, size(x) - 1
      intval = intval + 0.5_wp * (x(i+1) - x(i)) * (y(i+1) + y(i))
    end do
  end function trapezoid

  pure subroutine sort_increasing(a)
    real(wp), intent(inout) :: a(:)
    real(wp) :: key
    integer :: i, j

    do i = 2, size(a)
      key = a(i)
      j = i - 1
      do while (j >= 1 .and. a(j) > key)
        a(j+1) = a(j)
        j = j - 1
      end do
      a(j+1) = key
    end do
  end subroutine sort_increasing

end module donu_mod
