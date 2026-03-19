module anderson_optimized
  use precision_mod, only: wp
  use para_mod, only: SDIV, MDIV
  implicit none
  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(lda, *), b(ldb, *)
      double precision, intent(inout) :: c(ldc, *)
    end subroutine dgemm
    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(lda, *), x(*)
      double precision, intent(inout) :: y(*)
    end subroutine dgemv
  end interface
contains

  subroutine anderson_accel_optimized(current_field, target_field, history_f, iter, m_hist, dif)
    real(wp), dimension(SDIV,MDIV), intent(inout) :: current_field
    real(wp), dimension(SDIV,MDIV), intent(in)    :: target_field
    real(wp), dimension(SDIV,MDIV,m_hist), intent(inout) :: history_f
    integer, intent(in) :: iter, m_hist
    real(wp), intent(in) :: dif
    
    real(wp), dimension(m_hist) :: gamma, svals
    real(wp), dimension(m_hist, m_hist) :: F_mat
    real(wp) :: residual(SDIV,MDIV), accel_field(SDIV,MDIV)
    real(wp) :: work_ls(10*m_hist)
    integer :: k, i, info, idx_curr, idx, N, rank_out
    real(wp), parameter :: blend = 0.3e0_wp, damping = 0.3e0_wp
    integer, parameter :: conservative_steps = 3
    real(wp), allocatable, save :: delta(:,:)

    !if ( dif < 1.e-6_wp) then
    !  current_field = damping * current_field + (1.0e0_wp - damping) * target_field
    !  return
    !end if

    N = SDIV*MDIV

    ! Allocate/resize delta buffer (persists across calls, size N x m_hist)
    if (.not. allocated(delta)) then
      allocate(delta(N, m_hist))
    else if (size(delta,1) /= N .or. size(delta,2) /= m_hist) then
      deallocate(delta)
      allocate(delta(N, m_hist))
    end if

    ! Compute residual
    residual = current_field - target_field
    idx_curr = modulo(iter - 1, m_hist) + 1

    ! Early iterations: simple damping, but still seed valid history.
    if (iter < conservative_steps) then
      history_f(:,:,idx_curr) = residual  ! store pre-update residual
      current_field = damping * current_field + (1.0e0_wp - damping) * target_field
      return
    end if

    ! Store current residual
    history_f(:,:,idx_curr) = residual

    ! Determine history depth
    k = min(iter - conservative_steps + 1, m_hist)

    ! Precompute delta(:,i) = residual - history_f(:,:,idx(i))
    ! Each history slice is read exactly once (vs k*(k+1)/2 times in the scalar loop)
    do i = 1, k
      idx = modulo(iter - i - 1, m_hist) + 1
      delta(:,i) = reshape(residual, [N]) - reshape(history_f(:,:,idx), [N])
    end do

    ! F_mat(1:k,1:k) = delta(:,1:k)^T @ delta(:,1:k)  — all k^2 dot products in one DGEMM
    call dgemm('T', 'N', k, k, N, 1.e0_wp, delta, N, delta, N, 0.e0_wp, F_mat, m_hist)

    ! gamma(1:k) = delta(:,1:k)^T @ residual
    call dgemv('T', N, k, 1.e0_wp, delta, N, residual(1,1), 1, 0.e0_wp, gamma, 1)

    call DPOSV('U', k, 1, F_mat, m_hist, gamma, m_hist, info)

    associate (unused => dif); end associate
    if (info == 0) then
      ! accel_field = target_field - delta(:,1:k) @ gamma(1:k)
      accel_field = target_field
      call dgemv('N', N, k, -1.e0_wp, delta, N, gamma, 1, 1.e0_wp, accel_field(1,1), 1)
      !current_field = damping * accel_field + (1.0e0_wp - damping) * current_field
      current_field = blend * current_field + (1.0e0_wp - blend) * accel_field
    else
      !write(*,*) "Fallback to damped Picard"
      current_field = blend * current_field + (1.0e0_wp - blend) * target_field
    end if
    
  end subroutine anderson_accel_optimized

end module anderson_optimized

module JFNK_mod
  use precision_mod, only: wp
  use para_mod,      only: SDIV, MDIV
  implicit none
  private

  public :: jfnk_options
  public :: jfnk_solve
  public :: pack_relaxation_fields
  public :: unpack_relaxation_fields
  public :: relaxation_vector_size

  ! ------------------------------------------------------------------
  type :: jfnk_options
    integer  :: gmres_restart   = 12
    integer  :: max_newton      = 8
    integer  :: max_line_search = 8
    real(wp) :: eta             = 5.e-2_wp
    real(wp) :: atol            = 1.e-10_wp
    real(wp) :: rtol            = 1.e-8_wp
    real(wp) :: fd_eps          = 1.e-7_wp
    real(wp) :: armijo_c1       = 1.e-4_wp
    real(wp) :: line_reduce     = 5.e-1_wp
    real(wp) :: min_alpha       = 1.e-4_wp
  end type jfnk_options

  abstract interface
    subroutine residual_operator(x, fx)
      import :: wp
      real(wp), intent(in)  :: x(:)
      real(wp), intent(out) :: fx(:)
    end subroutine residual_operator
  end interface

contains

  ! ------------------------------------------------------------------
  pure integer function relaxation_vector_size(include_scalar) result(n)
    logical, intent(in), optional :: include_scalar
    logical :: use_scalar
    use_scalar = .true.
    if (present(include_scalar)) use_scalar = include_scalar
    n = merge(4, 3, use_scalar) * SDIV * MDIV
  end function relaxation_vector_size

  ! ------------------------------------------------------------------
  subroutine pack_relaxation_fields(rho, gama, ww, sphi, x, include_scalar)
    real(wp), intent(in),  contiguous :: rho(:,:), gama(:,:), ww(:,:), sphi(:,:)
    real(wp), intent(out), contiguous :: x(:)
    logical,  intent(in),  optional   :: include_scalar

    logical :: use_scalar
    integer :: nb   ! points per field block

    use_scalar = .true.
    if (present(include_scalar)) use_scalar = include_scalar
    nb = SDIV * MDIV

    if (size(x) /= relaxation_vector_size(use_scalar)) &
      error stop "pack_relaxation_fields: size mismatch"

    associate( xr => x(     1:  nb), &
               xg => x(  nb+1:2*nb), &
               xw => x(2*nb+1:3*nb)  )
      xr = reshape(rho,  [nb])
      xg = reshape(gama, [nb])
      xw = reshape(ww,   [nb])
    end associate

    if (use_scalar) x(3*nb+1:4*nb) = reshape(sphi, [nb])
  end subroutine pack_relaxation_fields

  ! ------------------------------------------------------------------
  subroutine unpack_relaxation_fields(x, rho, gama, ww, sphi, include_scalar)
    real(wp), intent(in),  contiguous :: x(:)
    real(wp), intent(out), contiguous :: rho(:,:), gama(:,:), ww(:,:), sphi(:,:)
    logical,  intent(in),  optional   :: include_scalar

    logical :: use_scalar
    integer :: nb

    use_scalar = .true.
    if (present(include_scalar)) use_scalar = include_scalar
    nb = SDIV * MDIV

    if (size(x) /= relaxation_vector_size(use_scalar)) &
      error stop "unpack_relaxation_fields: size mismatch"

    associate( xr => x(     1:  nb), &
               xg => x(  nb+1:2*nb), &
               xw => x(2*nb+1:3*nb)  )
      rho  = reshape(xr, [SDIV, MDIV])
      gama = reshape(xg, [SDIV, MDIV])
      ww   = reshape(xw, [SDIV, MDIV])
    end associate

    if (use_scalar) then
      sphi = reshape(x(3*nb+1:4*nb), [SDIV, MDIV])
    else
      sphi = 0.0_wp
    end if
  end subroutine unpack_relaxation_fields

  ! ------------------------------------------------------------------
  subroutine jfnk_solve(x, evaluate_residual, opts, info, final_res_norm, n_newton)
    real(wp),             intent(inout) :: x(:)
    procedure(residual_operator)        :: evaluate_residual
    type(jfnk_options),  intent(in), optional :: opts
    integer,             intent(out)          :: info
    real(wp),            intent(out), optional :: final_res_norm
    integer,             intent(out), optional :: n_newton

    type(jfnk_options)  :: cfg
    real(wp), allocatable :: fx(:), step(:), x_trial(:), f_trial(:)
    real(wp) :: res0, res_norm, target_tol, alpha
    integer  :: iter, gmres_info
    logical  :: accepted

    cfg = jfnk_options()
    if (present(opts)) cfg = opts

    allocate( fx(size(x)), step(size(x)), x_trial(size(x)), f_trial(size(x)) )

    call evaluate_residual(x, fx)
    res0       = norm2(fx)
    res_norm   = res0
    target_tol = max(cfg%atol, cfg%rtol * max(1.0_wp, res0))

    info = 1  ! not converged unless proven otherwise
    do iter = 1, cfg%max_newton
      if (res_norm <= target_tol) then
        info = 0
        exit
      end if

      call gmres_jfnk_step(x, fx, evaluate_residual, cfg, step, gmres_info)
      if (gmres_info /= 0) then
        info = 2; exit
      end if

      call backtracking_update(x, fx, step, evaluate_residual, cfg, &
                               x_trial, f_trial, alpha, accepted)
      if (.not. accepted) then
        info = 3; exit
      end if

      x        = x_trial
      fx       = f_trial
      res_norm = norm2(fx)
    end do

    ! Final convergence check after last iteration
    if (info == 1 .and. res_norm <= target_tol) info = 0

    if (present(final_res_norm)) final_res_norm = res_norm
    if (present(n_newton))       n_newton = iter - 1   ! completed steps
  end subroutine jfnk_solve

  ! ------------------------------------------------------------------
  subroutine gmres_jfnk_step(x_base, f_base, evaluate_residual, cfg, step, info)
    real(wp),           intent(in),  contiguous :: x_base(:), f_base(:)
    type(jfnk_options), intent(in)              :: cfg
    real(wp),           intent(out), contiguous :: step(:)
    integer,            intent(out)             :: info
    procedure(residual_operator)                :: evaluate_residual

    integer :: n, m, j, i, k
    real(wp) :: beta, lin_tol, resid, h_ij, h_jp1j
    real(wp), allocatable :: v(:,:), h(:,:), g(:), cs(:), sn(:), w(:), y(:), dscale(:)

    n = size(x_base)
    beta = norm2(f_base)

    if (beta <= cfg%atol) then
      step = 0.0_wp;  info = 0;  return
    end if

    m       = min(cfg%gmres_restart, n)
    lin_tol = max(cfg%atol, cfg%eta * beta)

    allocate( v(n,m+1), source=0.0_wp )
    allocate( h(m+1,m), source=0.0_wp )
    allocate( g(m+1),   source=0.0_wp )
    allocate( cs(m),    source=0.0_wp )
    allocate( sn(m),    source=0.0_wp )
    allocate( w(n), y(m), dscale(n) )

    call build_block_variable_scaling(x_base, dscale)

    v(:,1) = -f_base / beta
    g(1)   = beta
    k      = 0

    krylov: do j = 1, m
      ! Apply scaling in the FD direction only (right-preconditioning)
      call apply_jacobian_fd(x_base, f_base, dscale*v(:,j), &
                             evaluate_residual, cfg%fd_eps, w)

      ! Modified Gram-Schmidt
      do i = 1, j
        h(i,j) = dot_product(w, v(:,i))
        w      = w - h(i,j) * v(:,i)
      end do

      h(j+1,j) = norm2(w)
      if (h(j+1,j) > epsilon(h(j+1,j))) v(:,j+1) = w / h(j+1,j)

      ! Apply previous Givens rotations
      do i = 1, j-1
        call apply_givens(cs(i), sn(i), h(i,j), h(i+1,j))
      end do

      ! Build and apply new Givens rotation
      call build_givens(h(j,j), h(j+1,j), cs(j), sn(j))
      call apply_givens(cs(j), sn(j), h(j,j), h(j+1,j))
      call apply_givens(cs(j), sn(j), g(j),   g(j+1))

      k     = j
      resid = abs(g(j+1))
      if (resid <= lin_tol) exit krylov
    end do krylov

    if (k == 0) then
      step = 0.0_wp;  info = 1;  return
    end if

    y(1:k) = g(1:k)
    call solve_upper_triangular(h(1:k,1:k), y(1:k), info)
    if (info /= 0) then
      step = 0.0_wp;  return
    end if

    ! Recover step in original (unscaled) space
    step = dscale * matmul(v(:,1:k), y(1:k))
    info = 0
  end subroutine gmres_jfnk_step

  ! ------------------------------------------------------------------
  subroutine build_block_variable_scaling(x_base, dscale)
    real(wp), intent(in),  contiguous :: x_base(:)
    real(wp), intent(out), contiguous :: dscale(:)

    integer  :: n, nb
    real(wp) :: s

    n  = size(x_base)
    nb = SDIV * MDIV

    if (n /= 3*nb .and. n /= 4*nb) &
      error stop "build_block_variable_scaling: unexpected state size"

    associate( xr => x_base(     1:  nb), &
               xg => x_base(  nb+1:2*nb), &
               xw => x_base(2*nb+1:3*nb)  )
      dscale(     1:  nb) = 1.0_wp / max(1.0_wp, maxval(abs(xr)))
      dscale(  nb+1:2*nb) = 1.0_wp / max(1.0_wp, maxval(abs(xg)))
      dscale(2*nb+1:3*nb) = 1.0_wp / max(1.0_wp, maxval(abs(xw)))
    end associate

    if (n == 4*nb) then
      s = 1.0_wp / max(1.0_wp, maxval(abs(x_base(3*nb+1:4*nb))))
      dscale(3*nb+1:4*nb) = s
    end if
  end subroutine build_block_variable_scaling

  ! ------------------------------------------------------------------
  subroutine backtracking_update(x, fx, step, evaluate_residual, cfg, &
                                 x_trial, f_trial, alpha, accepted)
    real(wp),           intent(in),  contiguous :: x(:), fx(:), step(:)
    type(jfnk_options), intent(in)              :: cfg
    real(wp),           intent(out), contiguous :: x_trial(:), f_trial(:)
    real(wp),           intent(out)             :: alpha
    logical,            intent(out)             :: accepted
    procedure(residual_operator)                :: evaluate_residual

    real(wp) :: phi0, phi_trial, slope0
    real(wp), allocatable :: jstep(:)
    integer  :: i

    phi0 = 0.5_wp * dot_product(fx, fx)

    allocate(jstep(size(step)))
    call apply_jacobian_fd(x, fx, step, evaluate_residual, cfg%fd_eps, jstep)
    slope0 = dot_product(fx, jstep)
    if (slope0 >= 0.0_wp) slope0 = -2.0_wp * phi0   ! ensure descent

    alpha    = 1.0_wp
    accepted = .false.

    do i = 1, cfg%max_line_search
      x_trial = x + alpha * step
      call evaluate_residual(x_trial, f_trial)
      phi_trial = 0.5_wp * dot_product(f_trial, f_trial)

      if (phi_trial <= phi0 + cfg%armijo_c1 * alpha * slope0) then
        accepted = .true.;  return
      end if

      alpha = alpha * cfg%line_reduce
      if (alpha < cfg%min_alpha) exit
    end do
  end subroutine backtracking_update

  ! ------------------------------------------------------------------
  subroutine apply_jacobian_fd(x_base, f_base, v, evaluate_residual, fd_eps, jv)
    real(wp), intent(in),  contiguous :: x_base(:), f_base(:), v(:)
    real(wp), intent(in)              :: fd_eps
    real(wp), intent(out), contiguous :: jv(:)
    procedure(residual_operator)      :: evaluate_residual

    real(wp), allocatable :: x_pert(:), f_pert(:)
    real(wp) :: eps_fd, v_norm, x_norm

    v_norm = norm2(v)
    if (v_norm <= epsilon(v_norm)) then
      jv = 0.0_wp;  return
    end if

    x_norm = norm2(x_base)
    ! Floor x_norm at 1 so the perturbation is sensible near the origin
    eps_fd = fd_eps * max(1.0_wp, x_norm) / v_norm

    allocate(x_pert(size(x_base)), f_pert(size(x_base)))
    x_pert = x_base + eps_fd * v
    call evaluate_residual(x_pert, f_pert)
    jv = (f_pert - f_base) / eps_fd
  end subroutine apply_jacobian_fd

  ! ------------------------------------------------------------------
  subroutine solve_upper_triangular(rmat, rhs, info)
    real(wp), intent(in)    :: rmat(:,:)
    real(wp), intent(inout) :: rhs(:)
    integer,  intent(out)   :: info

    integer :: i, n
    n    = size(rhs)
    info = 0

    do i = n, 1, -1
      if (abs(rmat(i,i)) <= epsilon(rmat(i,i))) then
        info = 1;  return
      end if
      if (i < n) rhs(i) = rhs(i) - dot_product(rmat(i, i+1:n), rhs(i+1:n))
      rhs(i) = rhs(i) / rmat(i,i)
    end do
  end subroutine solve_upper_triangular

  ! ------------------------------------------------------------------
  subroutine build_givens(a, b, c, s)
    real(wp), intent(in)  :: a, b
    real(wp), intent(out) :: c, s
    real(wp) :: scale

    ! Use tiny() not epsilon(): epsilon tests relative magnitude,
    ! tiny() tests absolute underflow — correct zero guard before division
    if (abs(b) <= tiny(b)) then
      c = 1.0_wp;  s = 0.0_wp
    else if (abs(a) <= tiny(a)) then
      c = 0.0_wp;  s = 1.0_wp
    else
      scale = sqrt(a*a + b*b)
      c = a / scale
      s = b / scale
    end if
  end subroutine build_givens

  ! ------------------------------------------------------------------
  pure subroutine apply_givens(c, s, a, b)
    real(wp), intent(in)    :: c, s
    real(wp), intent(inout) :: a, b
    real(wp) :: tmp
    tmp = c*a + s*b
    b   = c*b - s*a     ! note: equivalent to -s*a + c*b
    a   = tmp
  end subroutine apply_givens

end module JFNK_mod

module aitken_mod
  use precision_mod, only: wp
  use para_mod, only: SDIV, MDIV
  implicit none
  private
  public :: aitken_delta2

  real(wp), allocatable, save :: aitk1_rho(:,:), aitk1_gama(:,:), aitk1_ww(:,:), aitk1_sphi(:,:)
  real(wp), allocatable, save :: aitk2_rho(:,:), aitk2_gama(:,:), aitk2_ww(:,:), aitk2_sphi(:,:)
  integer, save :: n_aitk_hist = 0

contains

  subroutine aitken_delta2(rho, gama, ww, sphi, prev_rho, prev_gama, prev_ww, prev_sphi, &
                           n_rho_locked, N_RHO_LOCK, has_scalar, fired)
    real(wp), intent(inout) :: rho(SDIV,MDIV), gama(SDIV,MDIV), ww(SDIV,MDIV), sphi(SDIV,MDIV)
    real(wp), intent(inout) :: prev_rho(SDIV,MDIV), prev_gama(SDIV,MDIV)
    real(wp), intent(inout) :: prev_ww(SDIV,MDIV), prev_sphi(SDIV,MDIV)
    integer,  intent(in)    :: n_rho_locked, N_RHO_LOCK
    logical,  intent(in)    :: has_scalar
    logical,  intent(out)   :: fired

    real(wp) :: x_k_rho(SDIV,MDIV), x_k_gama(SDIV,MDIV), x_k_ww(SDIV,MDIV), x_k_sphi(SDIV,MDIV)
    real(wp) :: denom, step_n
    real(wp), parameter :: MAX_AITKEN_JUMP = 10.e0_wp
    integer  :: s, m

    ! Allocate on first call or grid change
    if (.not. allocated(aitk1_rho) .or. size(aitk1_rho,1) /= SDIV .or. size(aitk1_rho,2) /= MDIV) then
      if (allocated(aitk1_rho)) deallocate(aitk1_rho, aitk1_gama, aitk1_ww, aitk1_sphi)
      if (allocated(aitk2_rho)) deallocate(aitk2_rho, aitk2_gama, aitk2_ww, aitk2_sphi)
      allocate(aitk1_rho(SDIV,MDIV),  source=rho)
      allocate(aitk1_gama(SDIV,MDIV), source=gama)
      allocate(aitk1_ww(SDIV,MDIV),   source=ww)
      allocate(aitk1_sphi(SDIV,MDIV), source=sphi)
      allocate(aitk2_rho(SDIV,MDIV),  source=rho)
      allocate(aitk2_gama(SDIV,MDIV), source=gama)
      allocate(aitk2_ww(SDIV,MDIV),   source=ww)
      allocate(aitk2_sphi(SDIV,MDIV), source=sphi)
      n_aitk_hist = 0
    end if

    fired = .false.
    x_k_rho = rho;  x_k_gama = gama;  x_k_ww = ww;  x_k_sphi = sphi

    if (n_rho_locked >= N_RHO_LOCK .and. n_aitk_hist >= 2) then
      fired = .true.

      do m = 1, MDIV
        do s = 1, SDIV
          step_n = rho(s,m) - aitk1_rho(s,m)
          denom  = rho(s,m) - 2.e0_wp*aitk1_rho(s,m) + aitk2_rho(s,m)
          if (abs(denom) > 1.e-14_wp .and. abs(step_n) <= MAX_AITKEN_JUMP * abs(denom)) &
            rho(s,m) = rho(s,m) - step_n**2 / denom

          step_n = gama(s,m) - aitk1_gama(s,m)
          denom  = gama(s,m) - 2.e0_wp*aitk1_gama(s,m) + aitk2_gama(s,m)
          if (abs(denom) > 1.e-14_wp .and. abs(step_n) <= MAX_AITKEN_JUMP * abs(denom)) &
            gama(s,m) = gama(s,m) - step_n**2 / denom

          step_n = ww(s,m) - aitk1_ww(s,m)
          denom  = ww(s,m) - 2.e0_wp*aitk1_ww(s,m) + aitk2_ww(s,m)
          if (abs(denom) > 1.e-14_wp .and. abs(step_n) <= MAX_AITKEN_JUMP * abs(denom)) &
            ww(s,m) = ww(s,m) - step_n**2 / denom
        end do
      end do

      if (has_scalar) then
        do m = 1, MDIV
          do s = 1, SDIV
            step_n = sphi(s,m) - aitk1_sphi(s,m)
            denom  = sphi(s,m) - 2.e0_wp*aitk1_sphi(s,m) + aitk2_sphi(s,m)
            if (abs(denom) > 1.e-14_wp .and. abs(step_n) <= MAX_AITKEN_JUMP * abs(denom)) &
              sphi(s,m) = sphi(s,m) - step_n**2 / denom
          end do
        end do
      end if

      where (abs(rho)  > 100.e0_wp) rho  = x_k_rho
      where (abs(gama) > 300.e0_wp) gama = x_k_gama
      where (abs(ww)   > 100.e0_wp) ww   = x_k_ww
      if (has_scalar) where (abs(sphi) > 10.e0_wp) sphi = x_k_sphi

      prev_rho = rho;  prev_gama = gama;  prev_ww = ww
      if (has_scalar) prev_sphi = sphi

      n_aitk_hist = 0
    end if

    ! Shift history: aitk2 <- aitk1 <- entry values (pre-modification)
    aitk2_rho  = aitk1_rho;   aitk2_gama = aitk1_gama;  aitk2_ww  = aitk1_ww;   aitk2_sphi = aitk1_sphi
    aitk1_rho  = x_k_rho;     aitk1_gama = x_k_gama;    aitk1_ww  = x_k_ww;     aitk1_sphi = x_k_sphi
    if (.not. fired) n_aitk_hist = n_aitk_hist + 1

  end subroutine aitken_delta2

end module aitken_mod