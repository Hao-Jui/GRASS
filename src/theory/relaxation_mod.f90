module anderson_optimized
  use precision_mod, only: wp
  use para_mod, only: SDIV, MDIV
  implicit none
  
contains

  subroutine anderson_accel_optimized(current_field, target_field, history_f, iter, m_hist)
    real(wp), dimension(SDIV,MDIV), intent(inout) :: current_field
    real(wp), dimension(SDIV,MDIV), intent(in)    :: target_field
    real(wp), dimension(SDIV,MDIV,m_hist), intent(inout) :: history_f
    integer, intent(in) :: iter, m_hist
    
    real(wp), dimension(m_hist) :: gamma
    real(wp), dimension(m_hist, m_hist) :: F_mat
    real(wp) :: residual(SDIV,MDIV), accel_field(SDIV,MDIV)
    integer :: k, i, info, idx_curr, idx, N
    real(wp), parameter :: blend = 0.1e0_wp, damping = 0.8e0_wp
    integer, parameter :: conservative_steps = 3
    real(wp), allocatable, save :: delta(:,:)

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
      history_f(:,:,idx_curr) = residual            ! store pre-update residual
      current_field = 0.5e0_wp * current_field + 0.5e0_wp * target_field
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

    if (info == 0) then
      ! accel_field = target_field - delta(:,1:k) @ gamma(1:k)
      accel_field = target_field
      call dgemv('N', N, k, -1.e0_wp, delta, N, gamma, 1, 1.e0_wp, accel_field(1,1), 1)
      current_field = damping * accel_field + (1.0e0_wp - damping) * current_field
    else
      ! Fallback to damped Picard
      current_field = blend * current_field + (1.0e0_wp - blend) * target_field
    end if
    
  end subroutine anderson_accel_optimized

end module anderson_optimized
