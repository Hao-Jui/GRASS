module nag_compat_mod
  use precision_mod, only: wp
  implicit none

  ! to include: F06QWF, D02NBF
  abstract interface
    subroutine nag_ode_rhs(t, y, dy)
      import wp
      real(wp), intent(in) :: t
      real(wp), intent(in) :: y(:)
      real(wp), intent(out) :: dy(:)
    end subroutine nag_ode_rhs
  end interface
contains

  subroutine d01gaf(x, y, n, ans, er, ifail)
    integer, intent(in)  :: n
    real(wp), intent(in)  :: x(n), y(n)
    real(wp), intent(out) :: ans, er
    integer, intent(out) :: ifail

    integer :: i, last
    real(wp) :: c, d1, d2, d3, h1, h2, h3, h4
    real(wp) :: r1, r2, r3, r4, s
    real(wp) :: first_step

    ifail = 0
    ans   = 0.0_wp
    er    = 0.0_wp

    if (n < 4) then
      ifail = 1
      return
    end if

    if (size(x) /= size(y)) then
      ifail = 2
      return
    end if

    first_step = x(2) - x(1)
    if (abs(first_step) < epsilon(first_step)) then
      ifail = 3
      return
    end if

    do i = 3, n
      h3 = x(i) - x(i-1)
      if (h3 * first_step <= 0.0_wp) then
        ifail = 3
        return
      end if
    end do

    ! piecewise cubic / higher-order Newton–Cotes adaptive-like interpolatory method
    h2 = first_step
    d3 = (y(2) - y(1)) / h2
    h3 = x(3) - x(2)
    d1 = (y(3) - y(2)) / h3
    h1 = h2 + h3
    d2 = (d1 - d3) / h1
    h4 = x(4) - x(3)
    r1 = (y(4) - y(3)) / h4
    r2 = (r1 - d1) / (h4 + h3)
    h1 = h1 + h4
    r3 = (r2 - d2) / h1

    ans = h2 * ( y(1) + h2 * ( d3/2.0_wp - h2 * ( d2/6.0_wp - (h2 + 2.0_wp*h3) * r3 / 12.0_wp ) ) )
    s   = -h2**3 * ( h2*(3.0_wp*h2 + 5.0_wp*h4) + 10.0_wp*h3*h1 ) / 60.0_wp
    r4  = 0.0_wp

    last = n - 1
    do i = 3, last
      ans = ans + h3 * ( (y(i) + y(i-1))/2.0_wp - h3*h3*(d2 + r2 + (h2 - h4)*r3)/12.0_wp )
      c   = h3**3 * ( 2.0_wp*h3*h3 + 5.0_wp*( h3*(h4 + h2) + 2.0_wp*h4*h2 ) ) / 120.0_wp
      er  = er + (c + s) * r4

      if (i == 3) then
        s = s + 2.0_wp * c
      else
        s = c
      end if

      if (i < last) then
        h1 = h2
        h2 = h3
        h3 = h4
        d1 = r1
        d2 = r2
        d3 = r3
        h4 = x(i+2) - x(i+1)
        r1 = (y(i+2) - y(i+1)) / h4
        r4 = h4 + h3
        r2 = (r1 - d1) / r4
        r4 = r4 + h2
        r3 = (r2 - d2) / r4
        r4 = r4 + h1
        r4 = (r3 - d3) / r4
      else
        ans = ans + h4 * ( y(n) - h4 * ( r1/2.0_wp + h4 * ( r2/6.0_wp + (2.0_wp*h3 + h4) * r3 / 12.0_wp ) ) )
        er  = er - h4**3 * r4 * ( h4*(3.0_wp*h4 + 5.0_wp*h2) + 10.0_wp*h3*(h2 + h3 + h4) ) / 60.0_wp + s * r4
        ans = ans + er
      end if
    end do
  end subroutine d01gaf

  subroutine d02pcf(f, neqn, y, yp, t, tout, relerr, abserr, flag, step_count, out, max_step)
    !! Runge-Kutta-Fehlberg 4(5) adaptive integrator.
    !! Integrates dy/dt = f(t,y) from t to tout.
    !! flag: in=1 (first call) or 2 (continuation); out=2 (success), 4 (eval limit), 6 (h<hmin), 8 (bad input).
    use ieee_arithmetic, only: ieee_is_nan
    procedure(nag_ode_rhs) :: f
    integer,  intent(in)    :: neqn
    real(wp), intent(inout) :: y(neqn), yp(neqn), t, relerr, abserr
    real(wp), intent(in)    :: tout
    integer,  intent(inout) :: flag
    integer,  intent(out)   :: step_count
    logical,  intent(in)    :: out
    real(wp), intent(in), optional :: max_step

    integer,  parameter :: MAX_STEPS = 200000
    real(wp), parameter :: HMIN = 1.0e-12_wp
    real(wp), parameter :: FAC_MAX = 5.0_wp, FAC_MIN = 0.1_wp, SAFETY = 0.9_wp
    real(wp), parameter :: PI_BETA = 0.4_wp / 5, PI_ALPHA = 0.7_wp / 5  ! PI-controller exponents
    real(wp) :: dir, h, hmax, err, err_prev, fac, fac_ceil, scale, dist, err_i
    real(wp) :: inv_sqrt_neqn
    real(wp) :: y5(neqn), k(neqn,5), w(neqn)
    integer :: i
    logical :: rejected

    ! --- Butcher tableau (RKF45) ---
    real(wp), parameter :: &
      a(5)   = [0.25_wp, 3.0_wp/8, 12.0_wp/13, 1.0_wp, 0.5_wp], &
      b21    = 0.25_wp, &
      b3(2)  = [3.0_wp/32, 9.0_wp/32], &
      b4(3)  = [1932.0_wp/2197, -7200.0_wp/2197, 7296.0_wp/2197], &
      b5(4)  = [439.0_wp/216, -8.0_wp, 3680.0_wp/513, -845.0_wp/4104], &
      b6(5)  = [-8.0_wp/27, 2.0_wp, -3544.0_wp/2565, 1859.0_wp/4104, -11.0_wp/40], &
      c5th(6)= [16.0_wp/135, 0.0_wp, 6656.0_wp/12825, 28561.0_wp/56430, -9.0_wp/50, 2.0_wp/55], &
      errc(6)= [1.0_wp/360, 0.0_wp, -128.0_wp/4275, -2197.0_wp/75240, 1.0_wp/50, 2.0_wp/55]

    ! --- Input validation ---
    if (neqn <= 0 .or. flag == 0 .or. abs(flag) > 2) then
      flag = 8; return
    end if
    relerr = max(relerr, 1.0e-12_wp)
    abserr = max(abserr, 1.0e-18_wp)

    dist = tout - t
    if (abs(dist) < epsilon(dist)) then
      call f(t, y, yp); flag = 2; return
    end if

    ! --- Initial step size ---
    dir  = sign(1.0_wp, dist)
    hmax = abs(dist)
    if (present(max_step)) then
      if (max_step > 0.0_wp) hmax = min(hmax, max_step)
    end if
    h = dir * max(1.0e-6_wp, min(abs(dist) * 0.1_wp, hmax))
    step_count = 0
    err_prev = 1.0e-4_wp
    err = 0.0_wp
    inv_sqrt_neqn = 1.0_wp / sqrt(real(neqn, wp))
    rejected = .false.
    call f(t, y, yp)

    ! --- Main integration loop ---
    do while (dir * (tout - t) > 0.0_wp)
      if (step_count >= MAX_STEPS) then
        if (out) call write_step_failure("d02pcf: maximum step count reached", t, h, err, step_count, y, yp)
        flag = 4; return
      end if
      h = dir * min(abs(h), abs(tout - t), hmax)
      if (abs(h) < HMIN) then
        if (out) call write_step_failure("d02pcf: step size fell below HMIN", t, h, err, step_count, y, yp)
        flag = 6; return
      end if

      ! RKF45 stages: reuse yp as k1 and keep only k2..k6 in workspace
      do i = 1, neqn
        w(i) = y(i) + h * b21 * yp(i)
      end do
      call f(t + a(1)*h, w, k(:,1))

      do i = 1, neqn
        w(i) = y(i) + h * (b3(1) * yp(i) + b3(2) * k(i,1))
      end do
      call f(t + a(2)*h, w, k(:,2))

      do i = 1, neqn
        w(i) = y(i) + h * (b4(1) * yp(i) + b4(2) * k(i,1) + b4(3) * k(i,2))
      end do
      call f(t + a(3)*h, w, k(:,3))

      do i = 1, neqn
        w(i) = y(i) + h * (b5(1) * yp(i) + b5(2) * k(i,1) + b5(3) * k(i,2) + b5(4) * k(i,3))
      end do
      call f(t + a(4)*h, w, k(:,4))

      do i = 1, neqn
        w(i) = y(i) + h * (b6(1) * yp(i) + b6(2) * k(i,1) + b6(3) * k(i,2) + b6(4) * k(i,3) + b6(5) * k(i,4))
      end do
      call f(t + a(5)*h, w, k(:,5))

      ! Error estimate (NaN-aware)
      err = 0.0_wp
      do i = 1, neqn
        y5(i) = y(i) + h * (c5th(1) * yp(i) + c5th(3) * k(i,2) + c5th(4) * k(i,3) + c5th(5) * k(i,4) + c5th(6) * k(i,5))
        err_i = h * (errc(1) * yp(i) + errc(3) * k(i,2) + errc(4) * k(i,3) + errc(5) * k(i,4) + errc(6) * k(i,5))
        scale = abserr + relerr * max(abs(y(i)), abs(y5(i)))
        err = max(err, abs(err_i) / scale)
      end do
      err = err * inv_sqrt_neqn
      if (ieee_is_nan(err)) then
        ! NaN detected: shrink step and retry
        if (out) call write_step_failure("d02pcf: NaN error estimate, retrying with smaller step", t, h, err, step_count, y, yp)
        h = dir * abs(h) * FAC_MIN
        rejected = .true.
        cycle
      end if

      if (err <= 1.0_wp) then
        ! Accept step — PI controller with post-rejection ceiling
        t = t + h
        y = y5
        step_count = step_count + 1
        call f(t, y, yp)
        if (out) write(*,"(99es27.17e3)") t, y, yp
        fac_ceil = merge(1.0_wp, FAC_MAX, rejected)
        fac = min(fac_ceil, SAFETY * err**(-PI_ALPHA) * err_prev**PI_BETA)
        err_prev = max(err, 1.0e-8_wp)
        rejected = .false.
      else
        ! Reject step
        fac = max(FAC_MIN, SAFETY * err**(-0.25_wp))
        rejected = .true.
      end if
      h = dir * min(abs(h) * fac, hmax)
    end do

    flag = 2
  contains
    subroutine write_step_failure(msg, t_now, h_now, err_now, nstep, y_now, yp_now)
      character(len=*), intent(in) :: msg
      real(wp), intent(in) :: t_now, h_now, err_now
      integer, intent(in) :: nstep
      real(wp), intent(in) :: y_now(:), yp_now(:)

      write(*,'(A)') trim(msg)
      write(*,'(A,1X,ES22.14)') "  t         :", t_now
      write(*,'(A,1X,ES22.14)') "  h         :", h_now
      write(*,'(A,1X,ES22.14)') "  err       :", err_now
      write(*,'(A,1X,I0)')      "  step_count:", nstep
      write(*,'(A,1X,*(ES22.14,1X))') "  y         :", y_now
      write(*,'(A,1X,*(ES22.14,1X))') "  yp        :", yp_now
    end subroutine write_step_failure
  end subroutine d02pcf

end module nag_compat_mod
