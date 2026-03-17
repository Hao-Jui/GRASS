subroutine loadEos
    use para_mod, only: p_at_PT, eos_file, num_tab, &
                        log_e, log_p, log_h, log_n0, &
                        enthalpy_min, C, KSCALE
    implicit none
    integer :: i

    p_at_PT = 1
    open(16,file="./eos/"//trim(adjustl(eos_file))//".dat",status="unknown",form='formatted')
    read(16,*) num_tab
    allocate(log_e(num_tab),log_p(num_tab),log_h(num_tab),log_n0(num_tab))
    do i = 1, num_tab
      read(16,*) log_e(i), log_p(i), log_h(i), log_n0(i)
      if ( i == 1 ) then
        enthalpy_min =   log( log_h(i) )
      else
        if ( log_p(i)/log_p(i-1) < 1.d0 + 1.d-15 ) then
          p_at_PT = i
          !write(*,*) log( log_p(i-1)*KSCALE ), log( log_p(i)*KSCALE ) 
          !write(*,*) log( log_e(i-1)*C*C*KSCALE ), log( log_e(i)*C*C*KSCALE ) 
        endif
      endif
    enddo
    close(16)
    
    log_e(:) = log( log_e(:)*C*C*KSCALE)
    log_p(:) = log( log_p(:)*KSCALE)
    log_h(:) = log( log( log_h(:) ) )
    log_n0(:)= log( log_n0(:))

    write(*,*) ' '
    write(*,*) '# EOS: ', eos_file
    write(*,*) 'EOS data is in with log(h_min) =',enthalpy_min
    
    if ( any(isnan(log_e)) .or. any(isnan(log_p)) .or. any(isnan(log_h)) ) stop "wrong table"

  end subroutine loadEos


  pure elemental real(8) function e_at_p(pp)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_e, num_tab, phase_transition
    implicit none
    real(8), intent(in) :: pp
    real(8) :: pwr
    if (phase_transition) then
      call interp_pt(log_p, log_e, num_tab, log(pp), pwr)
    else
      call interp(log_p, log_e, num_tab, log(pp), pwr)
    end if
    e_at_p = exp(pwr)
  end function e_at_p

  pure elemental real(8) function p_at_e(ee)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_e, num_tab, phase_transition
    implicit none
    real(8), intent(in) :: ee
    real(8) :: pwr
    if (phase_transition) then
      call interp_pt(log_e, log_p, num_tab, log(ee), pwr)
    else
      call interp(log_e, log_p, num_tab, log(ee), pwr)
    end if
    p_at_e = exp(pwr)
  end function p_at_e

  pure elemental function p_at_e_dual(ee) result(res)
    use toolkit_mod, only: interp_dual
    use ad_mod, only: dual, log, exp
    use para_mod, only : log_p, log_e, num_tab
    implicit none
    type(dual), intent(in) :: ee
    type(dual) :: pwr
    type(dual) :: res

    call interp_dual(log_e, log_p, num_tab, log(ee), pwr)
    res = exp(pwr)
  end function p_at_e_dual

  pure elemental real(8) function n0_at_e(ee)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_n0, log_e, num_tab, phase_transition
    implicit none
    real(8), intent(in) :: ee
    real(8) :: pwr
    if (phase_transition) then
      call interp_pt(log_e, log_n0, num_tab, log(ee), pwr)
    else
      call interp(log_e, log_n0, num_tab, log(ee), pwr)
    end if
    n0_at_e = exp(pwr)
  end function n0_at_e

  pure elemental real(8) function n0_at_h(hh)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_n0, log_h, num_tab, phase_transition
    implicit none
    real(8), intent(in) :: hh
    real(8) :: pwr
    if (phase_transition) then
      call interp_pt(log_h, log_n0, num_tab, log(hh), pwr)
    else
      call interp(log_h, log_n0, num_tab, log(hh), pwr)
    end if
    n0_at_h = exp(pwr)
  end function n0_at_h

  pure elemental real(8) function e_at_h(hh)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_e, log_h, num_tab, phase_transition
    implicit none
    real(8), intent(in) :: hh
    real(8) :: pwr
    if (phase_transition) then
      call interp_pt(log_h, log_e, num_tab, log(hh), pwr)
    else
      call interp(log_h, log_e, num_tab, log(hh), pwr)
    end if
    e_at_h = exp(pwr)
  end function e_at_h

  pure elemental real(8) function p_at_h(hh)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_h, num_tab
    implicit none
    real(8), intent(in) :: hh
    real(8) :: pwr
    call interp(log_h, log_p, num_tab, log(hh), pwr)
    p_at_h = exp(pwr)
  end function p_at_h

  pure elemental real(8) function h_at_p(pp)
    use toolkit_mod, only: interp, interp_pt
    use para_mod, only : log_p, log_h, num_tab
    implicit none
    real(8), intent(in) :: pp
    real(8) :: pwr
    call interp(log_p, log_h, num_tab, log(pp), pwr)
    h_at_p = exp(pwr)
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
    real(8), intent(in) :: ee
    integer, intent(in) :: n
    real(8), intent(out) :: derivative
    integer, intent(out), optional :: status

    integer :: info, idx, half_width, left, right, n_points, i
    real(8) :: x0, min_e, max_e
    real(8), allocatable :: nodes(:), values(:), coeffs(:, :)
    real(8) :: ee_clamped
    logical :: fatal_error
    real(8), parameter :: clamp_tol = 1.d-12

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
      real(8), intent(in) :: x0_local
      real(8), intent(in) :: x(m)
      real(8), intent(out) :: c(m, n_deriv + 1)

      integer :: i, j, k, max_k
      real(8) :: c1, c2, c3, c4, c5

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
