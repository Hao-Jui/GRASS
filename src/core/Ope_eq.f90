module ope_eq_mod
  use precision_mod, only: wp
  use para_mod, only: SDIV, MDIV, DS, DM, r_e, s_gp, mu
  implicit none
  private

  public :: laplacian_operator, gradient_vector

  type :: gradient_vector
    real(wp), allocatable :: r(:,:), t(:,:)
  end type gradient_vector

  type :: laplacian_operator
  contains
    procedure :: deriv_r
    procedure :: deriv_rr
    procedure :: deriv_t
    procedure :: deriv_tt
    procedure :: divr
    procedure :: grad
    procedure :: scal
    procedure :: lap2
    procedure :: laplacian
  end type laplacian_operator

contains

  function deriv_r(self, f) result(df_dr)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dr(size(f))
    real(wp) :: df_ds(size(f))
    real(wp) :: ds_dr(size(f))

    if (size(f) /= SDIV) stop "deriv_r: size mismatch"

    df_ds = deriv_s_1d(f)
    ds_dr = ((1.e0_wp - s_gp)**2) / max(r_e, 1.e-30_wp)
    df_dr = df_ds * ds_dr
  end function deriv_r

  function deriv_rr(self, f) result(df_drr)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_drr(size(f))
    real(wp) :: df_ds(size(f)), d2f_ds2(size(f))
    real(wp) :: ds_dr(size(f)), d2s_dr2(size(f))

    if (size(f) /= SDIV) stop "deriv_rr: size mismatch"

    df_ds = deriv_s_1d(f)
    d2f_ds2 = deriv_ss_1d(f)
    ds_dr = ((1.e0_wp - s_gp)**2) / max(r_e, 1.e-30_wp)
    d2s_dr2 = -2.e0_wp * (1.e0_wp - s_gp)**3 / max(r_e**2, 1.e-30_wp)
    df_drr = d2f_ds2 * ds_dr**2 + df_ds * d2s_dr2
  end function deriv_rr

  function deriv_t(self, f) result(df_dt)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dt(size(f))
    real(wp) :: df_dmu(size(f))

    if (size(f) /= MDIV) stop "deriv_t: size mismatch"

    df_dmu = deriv_mu_1d(f)
    df_dt = -sqrt(max(0.e0_wp, 1.e0_wp - mu**2)) * df_dmu
  end function deriv_t

  function deriv_tt(self, f) result(df_dtt)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dtt(size(f))
    real(wp) :: df_dmu(size(f)), d2f_dmu2(size(f))

    if (size(f) /= MDIV) stop "deriv_tt: size mismatch"

    df_dmu = deriv_mu_1d(f)
    d2f_dmu2 = deriv_mumu_1d(f)
    df_dtt = (1.e0_wp - mu**2) * d2f_dmu2 - mu * df_dmu
  end function deriv_tt

  function divr(self, f) result(df_over_r)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:,:)
    real(wp) :: df_over_r(size(f,1), size(f,2))
    real(wp) :: r_phys(SDIV)

    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "divr: size mismatch"

    r_phys = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)
    where (spread(r_phys, dim=2, ncopies=MDIV) < 1.e-30_wp)
      df_over_r = 0.e0_wp
    elsewhere
      df_over_r = f / spread(r_phys, dim=2, ncopies=MDIV)
    end where
  end function divr

  function grad(self, f) result(g)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:,:)
    type(gradient_vector) :: g
    integer :: s, m

    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "grad: size mismatch"

    allocate(g%r(SDIV,MDIV), g%t(SDIV,MDIV))
    do m = 1, MDIV
      g%r(:,m) = self%deriv_r(f(:,m))
    end do
    do s = 1, SDIV
      g%t(s,:) = self%deriv_t(f(s,:))
    end do
  end function grad

  function scal(self, x, y) result(prod)
    class(laplacian_operator), intent(in) :: self
    type(gradient_vector), intent(in) :: x, y
    real(wp) :: prod(size(x%r,1), size(x%r,2))
    real(wp) :: r_phys(SDIV)

    if (.not. allocated(x%r) .or. .not. allocated(x%t)) stop "scal: x gradient not allocated"
    if (.not. allocated(y%r) .or. .not. allocated(y%t)) stop "scal: y gradient not allocated"
    if (size(x%r,1) /= SDIV .or. size(x%r,2) /= MDIV) stop "scal: x%r size mismatch"
    if (size(x%t,1) /= SDIV .or. size(x%t,2) /= MDIV) stop "scal: x%t size mismatch"
    if (size(y%r,1) /= SDIV .or. size(y%r,2) /= MDIV) stop "scal: y%r size mismatch"
    if (size(y%t,1) /= SDIV .or. size(y%t,2) /= MDIV) stop "scal: y%t size mismatch"

    r_phys = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)

    prod = x%r * y%r
    where (spread(r_phys**2, dim=2, ncopies=MDIV) >= 1.e-30_wp)
      prod = prod + x%t * y%t / spread(r_phys**2, dim=2, ncopies=MDIV)
    end where
  end function scal

  function lap2(self, f) result(lap)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:,:)
    real(wp) :: lap(size(f,1), size(f,2))
    real(wp) :: r_phys(SDIV)
    real(wp) :: radial_r(SDIV,MDIV), radial_rr(SDIV,MDIV), theta_tt(SDIV,MDIV)
    integer :: s, m

    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "lap2: size mismatch"

    r_phys = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)

    do m = 1, MDIV
      radial_r(:,m) = self%deriv_r(f(:,m))
      radial_rr(:,m) = self%deriv_rr(f(:,m))
    end do
    do s = 1, SDIV
      theta_tt(s,:) = self%deriv_tt(f(s,:))
    end do

    lap = radial_rr
    where (spread(r_phys, dim=2, ncopies=MDIV) >= 1.e-30_wp)
      lap = lap + radial_r / spread(r_phys, dim=2, ncopies=MDIV)
    end where
    where (spread(r_phys**2, dim=2, ncopies=MDIV) >= 1.e-30_wp)
      lap = lap + theta_tt / spread(r_phys**2, dim=2, ncopies=MDIV)
    end where
  end function lap2

  function laplacian(self, f) result(lap)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:,:)
    real(wp) :: lap(size(f,1), size(f,2))
    real(wp) :: r_phys(SDIV)
    real(wp) :: radial_r(SDIV,MDIV), theta_term(SDIV,MDIV)
    integer :: s, m

    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "laplacian: size mismatch"

    r_phys = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)
    lap = self%lap2(f)

    ! lap2 contributes (1/r)f_r; add another (1/r)f_r to reach (2/r)f_r for 3D flat Laplacian
    do m = 1, MDIV
      radial_r(:,m) = self%deriv_r(f(:,m))
    end do
    where (spread(r_phys, dim=2, ncopies=MDIV) >= 1.e-30_wp)
      lap = lap + radial_r / spread(r_phys, dim=2, ncopies=MDIV)
    end where

    ! cot(theta)/r^2 * f_theta = -mu * f_mu / r^2
    do s = 1, SDIV
      theta_term(s,:) = -mu * deriv_mu_1d(f(s,:))
    end do
    where (spread(r_phys**2, dim=2, ncopies=MDIV) >= 1.e-30_wp)
      lap = lap + theta_term / spread(r_phys**2, dim=2, ncopies=MDIV)
    end where
  end function laplacian

  ! All four derivative helpers use 6th-order stencils (7-point).
  ! Interior:  centered 7-point,  O(h^6) for both f' and f''.
  ! Boundary:  one-sided 7-point, O(h^6) for f', O(h^5) for f''.
  ! Fallback to 4th-order when n < 7, and to 1st/2nd-order when n < 5.

  pure function deriv_s_1d(f) result(df_ds)
    real(wp), intent(in) :: f(:)
    real(wp) :: df_ds(size(f))
    integer :: n

    n = size(f)
    if (n <= 1) then
      df_ds = 0.e0_wp
      return
    end if

    if (n < 5) then
      df_ds(1) = (f(2) - f(1)) / DS
      if (n > 2) df_ds(2:n-1) = (f(3:n) - f(1:n-2)) / (2.e0_wp * DS)
      df_ds(n) = (f(n) - f(n-1)) / DS
      return
    end if

    if (n < 7) then
      df_ds(1) = (-25.e0_wp*f(1) + 48.e0_wp*f(2) - 36.e0_wp*f(3) + 16.e0_wp*f(4) - 3.e0_wp*f(5)) / (12.e0_wp*DS)
      df_ds(2) = (-3.e0_wp*f(1) - 10.e0_wp*f(2) + 18.e0_wp*f(3) - 6.e0_wp*f(4) + f(5)) / (12.e0_wp*DS)
      df_ds(3:n-2) = (-f(5:n) + 8.e0_wp*f(4:n-1) - 8.e0_wp*f(2:n-3) + f(1:n-4)) / (12.e0_wp*DS)
      df_ds(n-1) = (3.e0_wp*f(n) + 10.e0_wp*f(n-1) - 18.e0_wp*f(n-2) + 6.e0_wp*f(n-3) - f(n-4)) / (12.e0_wp*DS)
      df_ds(n) = (25.e0_wp*f(n) - 48.e0_wp*f(n-1) + 36.e0_wp*f(n-2) - 16.e0_wp*f(n-3) + 3.e0_wp*f(n-4)) / (12.e0_wp*DS)
      return
    end if

    ! 6th-order one-sided at left boundary (7-point forward)
    df_ds(1) = (-147.e0_wp*f(1) + 360.e0_wp*f(2) - 450.e0_wp*f(3) + 400.e0_wp*f(4) &
                -225.e0_wp*f(5) +  72.e0_wp*f(6) -  10.e0_wp*f(7)) / (60.e0_wp*DS)
    df_ds(2) = ( -10.e0_wp*f(1) -  77.e0_wp*f(2) + 150.e0_wp*f(3) - 100.e0_wp*f(4) &
                +  50.e0_wp*f(5) -  15.e0_wp*f(6) +   2.e0_wp*f(7)) / (60.e0_wp*DS)
    df_ds(3) = (   2.e0_wp*f(1) -  24.e0_wp*f(2) -  35.e0_wp*f(3) +  80.e0_wp*f(4) &
                -  30.e0_wp*f(5) +   8.e0_wp*f(6) -        f(7)) / (60.e0_wp*DS)

    ! 6th-order centered interior
    df_ds(4:n-3) = (  -f(1:n-6) + 9.e0_wp*f(2:n-5) - 45.e0_wp*f(3:n-4) &
                   +45.e0_wp*f(5:n-2) - 9.e0_wp*f(6:n-1) +  f(7:n)) / (60.e0_wp*DS)

    ! 6th-order one-sided at right boundary (7-point backward)
    df_ds(n-2) = (        f(n-6) -   8.e0_wp*f(n-5) +  30.e0_wp*f(n-4) - 80.e0_wp*f(n-3) &
                 + 35.e0_wp*f(n-2) +  24.e0_wp*f(n-1) -   2.e0_wp*f(n)) / (60.e0_wp*DS)
    df_ds(n-1) = (  -2.e0_wp*f(n-6) +  15.e0_wp*f(n-5) -  50.e0_wp*f(n-4) + 100.e0_wp*f(n-3) &
                 -150.e0_wp*f(n-2) +  77.e0_wp*f(n-1) +  10.e0_wp*f(n)) / (60.e0_wp*DS)
    df_ds(n)   = (  10.e0_wp*f(n-6) -  72.e0_wp*f(n-5) + 225.e0_wp*f(n-4) - 400.e0_wp*f(n-3) &
                 +450.e0_wp*f(n-2) - 360.e0_wp*f(n-1) + 147.e0_wp*f(n)) / (60.e0_wp*DS)
  end function deriv_s_1d

  pure function deriv_ss_1d(f) result(d2f_ds2)
    real(wp), intent(in) :: f(:)
    real(wp) :: d2f_ds2(size(f))
    integer :: n

    n = size(f)
    if (n <= 2) then
      d2f_ds2 = 0.e0_wp
      return
    end if

    if (n < 5) then
      d2f_ds2(1) = (f(3) - 2.e0_wp*f(2) + f(1)) / (DS**2)
      if (n > 3) d2f_ds2(2:n-1) = (f(3:n) - 2.e0_wp*f(2:n-1) + f(1:n-2)) / (DS**2)
      d2f_ds2(n) = (f(n) - 2.e0_wp*f(n-1) + f(n-2)) / (DS**2)
      return
    end if

    if (n < 7) then
      d2f_ds2(1) = (35.e0_wp*f(1) - 104.e0_wp*f(2) + 114.e0_wp*f(3) - 56.e0_wp*f(4) + 11.e0_wp*f(5)) / (12.e0_wp*DS**2)
      d2f_ds2(2) = (11.e0_wp*f(1) - 20.e0_wp*f(2) + 6.e0_wp*f(3) + 4.e0_wp*f(4) - f(5)) / (12.e0_wp*DS**2)
      d2f_ds2(3:n-2) = (-f(5:n) + 16.e0_wp*f(4:n-1) - 30.e0_wp*f(3:n-2) + 16.e0_wp*f(2:n-3) - f(1:n-4)) / (12.e0_wp*DS**2)
      d2f_ds2(n-1) = (11.e0_wp*f(n) - 20.e0_wp*f(n-1) + 6.e0_wp*f(n-2) + 4.e0_wp*f(n-3) - f(n-4)) / (12.e0_wp*DS**2)
      d2f_ds2(n) = (35.e0_wp*f(n) - 104.e0_wp*f(n-1) + 114.e0_wp*f(n-2) - 56.e0_wp*f(n-3) + 11.e0_wp*f(n-4)) / (12.e0_wp*DS**2)
      return
    end if

    ! 5th-order one-sided at left boundary (7-point forward)
    d2f_ds2(1) = ( 812.e0_wp*f(1) - 3132.e0_wp*f(2) + 5265.e0_wp*f(3) - 5080.e0_wp*f(4) &
                 +2970.e0_wp*f(5) -  972.e0_wp*f(6) +  137.e0_wp*f(7)) / (180.e0_wp*DS**2)
    d2f_ds2(2) = ( 137.e0_wp*f(1) -  147.e0_wp*f(2) -  255.e0_wp*f(3) +  470.e0_wp*f(4) &
                 - 285.e0_wp*f(5) +   93.e0_wp*f(6) -   13.e0_wp*f(7)) / (180.e0_wp*DS**2)
    d2f_ds2(3) = ( -13.e0_wp*f(1) +  228.e0_wp*f(2) -  420.e0_wp*f(3) +  200.e0_wp*f(4) &
                 +  15.e0_wp*f(5) -   12.e0_wp*f(6) +    2.e0_wp*f(7)) / (180.e0_wp*DS**2)

    ! 6th-order centered interior
    d2f_ds2(4:n-3) = (  2.e0_wp*f(1:n-6) - 27.e0_wp*f(2:n-5) + 270.e0_wp*f(3:n-4) - 490.e0_wp*f(4:n-3) &
                     +270.e0_wp*f(5:n-2) - 27.e0_wp*f(6:n-1) +   2.e0_wp*f(7:n)) / (180.e0_wp*DS**2)

    ! 5th-order one-sided at right boundary (7-point backward, even derivative: same coefficients reversed)
    d2f_ds2(n-2) = (   2.e0_wp*f(n-6) -   12.e0_wp*f(n-5) +   15.e0_wp*f(n-4) + 200.e0_wp*f(n-3) &
                   - 420.e0_wp*f(n-2) +  228.e0_wp*f(n-1) -   13.e0_wp*f(n)) / (180.e0_wp*DS**2)
    d2f_ds2(n-1) = ( -13.e0_wp*f(n-6) +   93.e0_wp*f(n-5) -  285.e0_wp*f(n-4) + 470.e0_wp*f(n-3) &
                   - 255.e0_wp*f(n-2) -  147.e0_wp*f(n-1) +  137.e0_wp*f(n)) / (180.e0_wp*DS**2)
    d2f_ds2(n)   = ( 137.e0_wp*f(n-6) -  972.e0_wp*f(n-5) + 2970.e0_wp*f(n-4) - 5080.e0_wp*f(n-3) &
                   +5265.e0_wp*f(n-2) - 3132.e0_wp*f(n-1) +  812.e0_wp*f(n)) / (180.e0_wp*DS**2)
  end function deriv_ss_1d

  pure function deriv_mu_1d(f) result(df_dmu)
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dmu(size(f))
    integer :: n

    n = size(f)
    if (n <= 1) then
      df_dmu = 0.e0_wp
      return
    end if

    if (n < 5) then
      df_dmu(1) = (f(2) - f(1)) / DM
      if (n > 2) df_dmu(2:n-1) = (f(3:n) - f(1:n-2)) / (2.e0_wp * DM)
      df_dmu(n) = (f(n) - f(n-1)) / DM
      return
    end if

    if (n < 7) then
      df_dmu(1) = (-25.e0_wp*f(1) + 48.e0_wp*f(2) - 36.e0_wp*f(3) + 16.e0_wp*f(4) - 3.e0_wp*f(5)) / (12.e0_wp*DM)
      df_dmu(2) = (-3.e0_wp*f(1) - 10.e0_wp*f(2) + 18.e0_wp*f(3) - 6.e0_wp*f(4) + f(5)) / (12.e0_wp*DM)
      df_dmu(3:n-2) = (-f(5:n) + 8.e0_wp*f(4:n-1) - 8.e0_wp*f(2:n-3) + f(1:n-4)) / (12.e0_wp*DM)
      df_dmu(n-1) = (3.e0_wp*f(n) + 10.e0_wp*f(n-1) - 18.e0_wp*f(n-2) + 6.e0_wp*f(n-3) - f(n-4)) / (12.e0_wp*DM)
      df_dmu(n) = (25.e0_wp*f(n) - 48.e0_wp*f(n-1) + 36.e0_wp*f(n-2) - 16.e0_wp*f(n-3) + 3.e0_wp*f(n-4)) / (12.e0_wp*DM)
      return
    end if

    df_dmu(1) = (-147.e0_wp*f(1) + 360.e0_wp*f(2) - 450.e0_wp*f(3) + 400.e0_wp*f(4) &
                 -225.e0_wp*f(5) +  72.e0_wp*f(6) -  10.e0_wp*f(7)) / (60.e0_wp*DM)
    df_dmu(2) = ( -10.e0_wp*f(1) -  77.e0_wp*f(2) + 150.e0_wp*f(3) - 100.e0_wp*f(4) &
                 +  50.e0_wp*f(5) -  15.e0_wp*f(6) +   2.e0_wp*f(7)) / (60.e0_wp*DM)
    df_dmu(3) = (   2.e0_wp*f(1) -  24.e0_wp*f(2) -  35.e0_wp*f(3) +  80.e0_wp*f(4) &
                 -  30.e0_wp*f(5) +   8.e0_wp*f(6) -        f(7)) / (60.e0_wp*DM)

    df_dmu(4:n-3) = (  -f(1:n-6) + 9.e0_wp*f(2:n-5) - 45.e0_wp*f(3:n-4) &
                    +45.e0_wp*f(5:n-2) - 9.e0_wp*f(6:n-1) +  f(7:n)) / (60.e0_wp*DM)

    df_dmu(n-2) = (        f(n-6) -   8.e0_wp*f(n-5) +  30.e0_wp*f(n-4) - 80.e0_wp*f(n-3) &
                  + 35.e0_wp*f(n-2) +  24.e0_wp*f(n-1) -   2.e0_wp*f(n)) / (60.e0_wp*DM)
    df_dmu(n-1) = (  -2.e0_wp*f(n-6) +  15.e0_wp*f(n-5) -  50.e0_wp*f(n-4) + 100.e0_wp*f(n-3) &
                  -150.e0_wp*f(n-2) +  77.e0_wp*f(n-1) +  10.e0_wp*f(n)) / (60.e0_wp*DM)
    df_dmu(n)   = (  10.e0_wp*f(n-6) -  72.e0_wp*f(n-5) + 225.e0_wp*f(n-4) - 400.e0_wp*f(n-3) &
                  +450.e0_wp*f(n-2) - 360.e0_wp*f(n-1) + 147.e0_wp*f(n)) / (60.e0_wp*DM)
  end function deriv_mu_1d

  pure function deriv_mumu_1d(f) result(d2f_dmu2)
    real(wp), intent(in) :: f(:)
    real(wp) :: d2f_dmu2(size(f))
    integer :: n

    n = size(f)
    if (n <= 2) then
      d2f_dmu2 = 0.e0_wp
      return
    end if

    if (n < 5) then
      d2f_dmu2(1) = (f(3) - 2.e0_wp*f(2) + f(1)) / (DM**2)
      if (n > 3) d2f_dmu2(2:n-1) = (f(3:n) - 2.e0_wp*f(2:n-1) + f(1:n-2)) / (DM**2)
      d2f_dmu2(n) = (f(n) - 2.e0_wp*f(n-1) + f(n-2)) / (DM**2)
      return
    end if

    if (n < 7) then
      d2f_dmu2(1) = (35.e0_wp*f(1) - 104.e0_wp*f(2) + 114.e0_wp*f(3) - 56.e0_wp*f(4) + 11.e0_wp*f(5)) / (12.e0_wp*DM**2)
      d2f_dmu2(2) = (11.e0_wp*f(1) - 20.e0_wp*f(2) + 6.e0_wp*f(3) + 4.e0_wp*f(4) - f(5)) / (12.e0_wp*DM**2)
      d2f_dmu2(3:n-2) = (-f(5:n) + 16.e0_wp*f(4:n-1) - 30.e0_wp*f(3:n-2) + 16.e0_wp*f(2:n-3) - f(1:n-4)) / (12.e0_wp*DM**2)
      d2f_dmu2(n-1) = (11.e0_wp*f(n) - 20.e0_wp*f(n-1) + 6.e0_wp*f(n-2) + 4.e0_wp*f(n-3) - f(n-4)) / (12.e0_wp*DM**2)
      d2f_dmu2(n) = (35.e0_wp*f(n) - 104.e0_wp*f(n-1) + 114.e0_wp*f(n-2) - 56.e0_wp*f(n-3) + 11.e0_wp*f(n-4)) / (12.e0_wp*DM**2)
      return
    end if

    ! 5th-order one-sided at left boundary (7-point forward)
    d2f_dmu2(1) = ( 812.e0_wp*f(1) - 3132.e0_wp*f(2) + 5265.e0_wp*f(3) - 5080.e0_wp*f(4) &
                 +2970.e0_wp*f(5) -  972.e0_wp*f(6) +  137.e0_wp*f(7)) / (180.e0_wp*DM**2)
    d2f_dmu2(2) = ( 137.e0_wp*f(1) -  147.e0_wp*f(2) -  255.e0_wp*f(3) +  470.e0_wp*f(4) &
                 - 285.e0_wp*f(5) +   93.e0_wp*f(6) -   13.e0_wp*f(7)) / (180.e0_wp*DM**2)
    d2f_dmu2(3) = ( -13.e0_wp*f(1) +  228.e0_wp*f(2) -  420.e0_wp*f(3) +  200.e0_wp*f(4) &
                 +  15.e0_wp*f(5) -   12.e0_wp*f(6) +    2.e0_wp*f(7)) / (180.e0_wp*DM**2)

    ! 6th-order centered interior
    d2f_dmu2(4:n-3) = (  2.e0_wp*f(1:n-6) - 27.e0_wp*f(2:n-5) + 270.e0_wp*f(3:n-4) - 490.e0_wp*f(4:n-3) &
                     +270.e0_wp*f(5:n-2) - 27.e0_wp*f(6:n-1) +   2.e0_wp*f(7:n)) / (180.e0_wp*DM**2)

    ! 5th-order one-sided at right boundary (7-point backward, even derivative: same coefficients reversed)
    d2f_dmu2(n-2) = (   2.e0_wp*f(n-6) -   12.e0_wp*f(n-5) +   15.e0_wp*f(n-4) + 200.e0_wp*f(n-3) &
                   - 420.e0_wp*f(n-2) +  228.e0_wp*f(n-1) -   13.e0_wp*f(n)) / (180.e0_wp*DM**2)
    d2f_dmu2(n-1) = ( -13.e0_wp*f(n-6) +   93.e0_wp*f(n-5) -  285.e0_wp*f(n-4) + 470.e0_wp*f(n-3) &
                   - 255.e0_wp*f(n-2) -  147.e0_wp*f(n-1) +  137.e0_wp*f(n)) / (180.e0_wp*DM**2)
    d2f_dmu2(n)   = ( 137.e0_wp*f(n-6) -  972.e0_wp*f(n-5) + 2970.e0_wp*f(n-4) - 5080.e0_wp*f(n-3) &
                   +5265.e0_wp*f(n-2) - 3132.e0_wp*f(n-1) +  812.e0_wp*f(n)) / (180.e0_wp*DM**2)
  end function deriv_mumu_1d

end module ope_eq_mod
