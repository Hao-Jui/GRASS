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

  ! -----------------------------------------------------------------------
  ! Public type-bound procedures
  ! -----------------------------------------------------------------------

  function deriv_r(self, f) result(df_dr)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dr(size(f))
    if (size(f) /= SDIV) stop "deriv_r: size mismatch"
    df_dr = deriv_s_1d(f) * ((1.e0_wp - s_gp)**2) / max(r_e, 1.e-30_wp)
  end function deriv_r

  function deriv_rr(self, f) result(df_drr)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_drr(size(f))
    real(wp) :: df_ds(SDIV), d2f_ds2(SDIV)
    real(wp) :: ds_dr(SDIV), d2s_dr2(SDIV)
    if (size(f) /= SDIV) stop "deriv_rr: size mismatch"
    call deriv_s_and_ss_1d(f, df_ds, d2f_ds2)
    ds_dr   = (1.e0_wp - s_gp)**2  / max(r_e,    1.e-30_wp)
    d2s_dr2 = -2.e0_wp * (1.e0_wp - s_gp)**3 / max(r_e**2, 1.e-30_wp)
    df_drr  = d2f_ds2 * ds_dr**2 + df_ds * d2s_dr2
  end function deriv_rr

  function deriv_t(self, f) result(df_dt)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dt(size(f))
    if (size(f) /= MDIV) stop "deriv_t: size mismatch"
    df_dt = -sqrt(max(0.e0_wp, 1.e0_wp - mu**2)) * deriv_mu_1d(f)
  end function deriv_t

  function deriv_tt(self, f) result(df_dtt)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dtt(size(f))
    real(wp) :: df_dmu(MDIV), d2f_dmu2(MDIV)
    if (size(f) /= MDIV) stop "deriv_tt: size mismatch"
    call deriv_mu_and_mumu_1d(f, df_dmu, d2f_dmu2)
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
    real(wp) :: f_T(MDIV,SDIV)
    real(wp) :: ds_dr(SDIV)
    integer :: s, m
    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "grad: size mismatch"
    allocate(g%r(SDIV,MDIV), g%t(SDIV,MDIV))
    ds_dr = (1.e0_wp - s_gp)**2 / max(r_e, 1.e-30_wp)
    do m = 1, MDIV
      g%r(:,m) = deriv_s_1d(f(:,m)) * ds_dr
    end do
    f_T = transpose(f)
    do s = 1, SDIV
      g%t(s,:) = -sqrt(max(0.e0_wp, 1.e0_wp - mu**2)) * deriv_mu_1d(f_T(:,s))
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

  ! lap2 = f_rr + (1/r) f_r + (1/r²) f_θθ   (KEH operator)
  ! Optimisations vs. original:
  !   - one combined radial pass (no separate radial_r / radial_rr temporaries)
  !   - transpose(f) for cache-friendly angular access
  !   - precomputed inv_r, inv_r2 (no repeated spread)
  function lap2(self, f) result(lap)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:,:)
    real(wp) :: lap(SDIV,MDIV)
    real(wp) :: r_phys(SDIV), ds_dr(SDIV), d2s_dr2(SDIV), inv_r(SDIV), inv_r2(SDIV)
    real(wp) :: df_ds(SDIV), d2f_ds2(SDIV)
    real(wp) :: df_dmu(MDIV), d2f_dmu2(MDIV)
    real(wp) :: f_T(MDIV,SDIV)
    integer :: s, m

    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "lap2: size mismatch"

    r_phys  = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)
    ds_dr   = (1.e0_wp - s_gp)**2 / max(r_e, 1.e-30_wp)
    d2s_dr2 = -2.e0_wp * (1.e0_wp - s_gp)**3 / max(r_e**2, 1.e-30_wp)
    where (r_phys >= 1.e-30_wp)
      inv_r = 1.e0_wp / r_phys
    elsewhere
      inv_r = 0.e0_wp
    end where
    inv_r2 = inv_r**2

    ! Radial: f_rr + (1/r) f_r  — one combined pass, no large intermediates
    do m = 1, MDIV
      call deriv_s_and_ss_1d(f(:,m), df_ds, d2f_ds2)
      lap(:,m) = d2f_ds2 * ds_dr**2 + df_ds * (d2s_dr2 + inv_r * ds_dr)
    end do

    ! Angular: (1/r²) f_θθ = (1/r²) [(1-μ²) f_μμ - μ f_μ]
    ! Transpose for contiguous row access (f_T(:,s) is stride-1)
    f_T = transpose(f)
    do s = 1, SDIV
      call deriv_mu_and_mumu_1d(f_T(:,s), df_dmu, d2f_dmu2)
      lap(s,:) = lap(s,:) + inv_r2(s) * ((1.e0_wp - mu**2) * d2f_dmu2 - mu * df_dmu)
    end do
  end function lap2

  ! laplacian = f_rr + (2/r) f_r + (1/r²) [(1-μ²) f_μμ - 2μ f_μ]   (full 3D flat)
  ! Does NOT call lap2 — computes all terms in one radial pass and one angular pass.
  function laplacian(self, f) result(lap)
    class(laplacian_operator), intent(in) :: self
    real(wp), intent(in) :: f(:,:)
    real(wp) :: lap(SDIV,MDIV)
    real(wp) :: r_phys(SDIV), ds_dr(SDIV), d2s_dr2(SDIV), inv_r(SDIV), inv_r2(SDIV)
    real(wp) :: df_ds(SDIV), d2f_ds2(SDIV)
    real(wp) :: df_dmu(MDIV), d2f_dmu2(MDIV)
    real(wp) :: f_T(MDIV,SDIV)
    integer :: s, m

    if (size(f,1) /= SDIV .or. size(f,2) /= MDIV) stop "laplacian: size mismatch"

    r_phys  = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)
    ds_dr   = (1.e0_wp - s_gp)**2 / max(r_e, 1.e-30_wp)
    d2s_dr2 = -2.e0_wp * (1.e0_wp - s_gp)**3 / max(r_e**2, 1.e-30_wp)
    where (r_phys >= 1.e-30_wp)
      inv_r = 1.e0_wp / r_phys
    elsewhere
      inv_r = 0.e0_wp
    end where
    inv_r2 = inv_r**2

    ! Radial: f_rr + (2/r) f_r  — fused into single stencil evaluation
    do m = 1, MDIV
      call deriv_s_and_ss_1d(f(:,m), df_ds, d2f_ds2)
      lap(:,m) = d2f_ds2 * ds_dr**2 + df_ds * (d2s_dr2 + 2.e0_wp * inv_r * ds_dr)
    end do

    ! Angular: (1/r²) [(1-μ²) f_μμ - 2μ f_μ]  (full 3D: f_θθ + cot θ f_θ) / r²
    f_T = transpose(f)
    do s = 1, SDIV
      call deriv_mu_and_mumu_1d(f_T(:,s), df_dmu, d2f_dmu2)
      lap(s,:) = lap(s,:) + inv_r2(s) * ((1.e0_wp - mu**2) * d2f_dmu2 - 2.e0_wp * mu * df_dmu)
    end do
  end function laplacian

  ! -----------------------------------------------------------------------
  ! Private combined-derivative subroutines (fused interior loop)
  ! -----------------------------------------------------------------------

  ! All four derivative helpers use 6th-order stencils (7-point).
  ! Interior:  centered 7-point,  O(h^6) for both f' and f''.
  ! Boundary:  one-sided 7-point, O(h^6) for f', O(h^5) for f''.
  ! Fallback to 4th-order when n < 7, and to 1st/2nd-order when n < 5.

  pure subroutine deriv_s_and_ss_1d(f, df_ds, d2f_ds2)
    real(wp), intent(in)  :: f(:)
    real(wp), intent(out) :: df_ds(:), d2f_ds2(:)
    integer  :: n, i
    real(wp) :: inv60DS, inv180DS2

    n = size(f)

    if (n <= 1) then
      df_ds = 0.e0_wp;  d2f_ds2 = 0.e0_wp
      return
    end if
    if (n <= 2) then
      df_ds(1) = (f(2) - f(1)) / DS;  df_ds(2) = df_ds(1)
      d2f_ds2  = 0.e0_wp
      return
    end if

    if (n < 5) then
      df_ds(1) = (f(2) - f(1)) / DS
      if (n > 2) df_ds(2:n-1) = (f(3:n) - f(1:n-2)) / (2.e0_wp * DS)
      df_ds(n) = (f(n) - f(n-1)) / DS
      d2f_ds2(1) = (f(3) - 2.e0_wp*f(2) + f(1)) / DS**2
      if (n > 3) d2f_ds2(2:n-1) = (f(3:n) - 2.e0_wp*f(2:n-1) + f(1:n-2)) / DS**2
      d2f_ds2(n) = (f(n) - 2.e0_wp*f(n-1) + f(n-2)) / DS**2
      return
    end if

    if (n < 7) then
      df_ds(1)   = (-25.e0_wp*f(1)+48.e0_wp*f(2)-36.e0_wp*f(3)+16.e0_wp*f(4)-3.e0_wp*f(5)) / (12.e0_wp*DS)
      df_ds(2)   = (-3.e0_wp*f(1)-10.e0_wp*f(2)+18.e0_wp*f(3)-6.e0_wp*f(4)+f(5)) / (12.e0_wp*DS)
      df_ds(3:n-2) = (-f(5:n)+8.e0_wp*f(4:n-1)-8.e0_wp*f(2:n-3)+f(1:n-4)) / (12.e0_wp*DS)
      df_ds(n-1) = (3.e0_wp*f(n)+10.e0_wp*f(n-1)-18.e0_wp*f(n-2)+6.e0_wp*f(n-3)-f(n-4)) / (12.e0_wp*DS)
      df_ds(n)   = (25.e0_wp*f(n)-48.e0_wp*f(n-1)+36.e0_wp*f(n-2)-16.e0_wp*f(n-3)+3.e0_wp*f(n-4)) / (12.e0_wp*DS)
      d2f_ds2(1)   = (35.e0_wp*f(1)-104.e0_wp*f(2)+114.e0_wp*f(3)-56.e0_wp*f(4)+11.e0_wp*f(5)) / (12.e0_wp*DS**2)
      d2f_ds2(2)   = (11.e0_wp*f(1)-20.e0_wp*f(2)+6.e0_wp*f(3)+4.e0_wp*f(4)-f(5)) / (12.e0_wp*DS**2)
      d2f_ds2(3:n-2) = (-f(5:n)+16.e0_wp*f(4:n-1)-30.e0_wp*f(3:n-2)+16.e0_wp*f(2:n-3)-f(1:n-4)) / (12.e0_wp*DS**2)
      d2f_ds2(n-1) = (11.e0_wp*f(n)-20.e0_wp*f(n-1)+6.e0_wp*f(n-2)+4.e0_wp*f(n-3)-f(n-4)) / (12.e0_wp*DS**2)
      d2f_ds2(n)   = (35.e0_wp*f(n)-104.e0_wp*f(n-1)+114.e0_wp*f(n-2)-56.e0_wp*f(n-3)+11.e0_wp*f(n-4)) / (12.e0_wp*DS**2)
      return
    end if

    inv60DS   = 1.e0_wp / (60.e0_wp * DS)
    inv180DS2 = 1.e0_wp / (180.e0_wp * DS**2)

    ! 6th-order one-sided boundaries (f')
    df_ds(1) = (-147.e0_wp*f(1)+360.e0_wp*f(2)-450.e0_wp*f(3)+400.e0_wp*f(4) &
                -225.e0_wp*f(5)+ 72.e0_wp*f(6)- 10.e0_wp*f(7)) * inv60DS
    df_ds(2) = ( -10.e0_wp*f(1)- 77.e0_wp*f(2)+150.e0_wp*f(3)-100.e0_wp*f(4) &
                +  50.e0_wp*f(5)- 15.e0_wp*f(6)+  2.e0_wp*f(7)) * inv60DS
    df_ds(3) = (   2.e0_wp*f(1)- 24.e0_wp*f(2)- 35.e0_wp*f(3)+ 80.e0_wp*f(4) &
                -  30.e0_wp*f(5)+  8.e0_wp*f(6)-       f(7)) * inv60DS
    df_ds(n-2) = (       f(n-6)-  8.e0_wp*f(n-5)+ 30.e0_wp*f(n-4)- 80.e0_wp*f(n-3) &
                 + 35.e0_wp*f(n-2)+ 24.e0_wp*f(n-1)-  2.e0_wp*f(n)) * inv60DS
    df_ds(n-1) = (  -2.e0_wp*f(n-6)+ 15.e0_wp*f(n-5)- 50.e0_wp*f(n-4)+100.e0_wp*f(n-3) &
                 -150.e0_wp*f(n-2)+ 77.e0_wp*f(n-1)+ 10.e0_wp*f(n)) * inv60DS
    df_ds(n)   = (  10.e0_wp*f(n-6)- 72.e0_wp*f(n-5)+225.e0_wp*f(n-4)-400.e0_wp*f(n-3) &
                 + 450.e0_wp*f(n-2)-360.e0_wp*f(n-1)+147.e0_wp*f(n)) * inv60DS

    ! 5th-order one-sided boundaries (f'')
    d2f_ds2(1) = ( 812.e0_wp*f(1)-3132.e0_wp*f(2)+5265.e0_wp*f(3)-5080.e0_wp*f(4) &
                 +2970.e0_wp*f(5)- 972.e0_wp*f(6)+ 137.e0_wp*f(7)) * inv180DS2
    d2f_ds2(2) = ( 137.e0_wp*f(1)- 147.e0_wp*f(2)- 255.e0_wp*f(3)+ 470.e0_wp*f(4) &
                 - 285.e0_wp*f(5)+  93.e0_wp*f(6)-  13.e0_wp*f(7)) * inv180DS2
    d2f_ds2(3) = ( -13.e0_wp*f(1)+ 228.e0_wp*f(2)- 420.e0_wp*f(3)+ 200.e0_wp*f(4) &
                 +  15.e0_wp*f(5)-  12.e0_wp*f(6)+   2.e0_wp*f(7)) * inv180DS2
    d2f_ds2(n-2) = (   2.e0_wp*f(n-6)-  12.e0_wp*f(n-5)+  15.e0_wp*f(n-4)+200.e0_wp*f(n-3) &
                   - 420.e0_wp*f(n-2)+ 228.e0_wp*f(n-1)-  13.e0_wp*f(n)) * inv180DS2
    d2f_ds2(n-1) = ( -13.e0_wp*f(n-6)+  93.e0_wp*f(n-5)- 285.e0_wp*f(n-4)+470.e0_wp*f(n-3) &
                   - 255.e0_wp*f(n-2)- 147.e0_wp*f(n-1)+ 137.e0_wp*f(n)) * inv180DS2
    d2f_ds2(n)   = ( 137.e0_wp*f(n-6)- 972.e0_wp*f(n-5)+2970.e0_wp*f(n-4)-5080.e0_wp*f(n-3) &
                   +5265.e0_wp*f(n-2)-3132.e0_wp*f(n-1)+ 812.e0_wp*f(n)) * inv180DS2

    ! Fused 6th-order centered interior: f(i-3:i+3) loaded once for both stencils
    do i = 4, n-3
      df_ds(i)   = (   -f(i-3) + 9.e0_wp*f(i-2) - 45.e0_wp*f(i-1) &
                   + 45.e0_wp*f(i+1) -  9.e0_wp*f(i+2) +       f(i+3)) * inv60DS
      d2f_ds2(i) = (2.e0_wp*f(i-3) - 27.e0_wp*f(i-2) + 270.e0_wp*f(i-1) - 490.e0_wp*f(i) &
                   +270.e0_wp*f(i+1) - 27.e0_wp*f(i+2) +  2.e0_wp*f(i+3)) * inv180DS2
    end do
  end subroutine deriv_s_and_ss_1d

  pure subroutine deriv_mu_and_mumu_1d(f, df_dmu, d2f_dmu2)
    real(wp), intent(in)  :: f(:)
    real(wp), intent(out) :: df_dmu(:), d2f_dmu2(:)
    integer  :: n, i
    real(wp) :: inv60DM, inv180DM2

    n = size(f)

    if (n <= 1) then
      df_dmu = 0.e0_wp;  d2f_dmu2 = 0.e0_wp
      return
    end if
    if (n <= 2) then
      df_dmu(1) = (f(2) - f(1)) / DM;  df_dmu(2) = df_dmu(1)
      d2f_dmu2  = 0.e0_wp
      return
    end if

    if (n < 5) then
      df_dmu(1) = (f(2) - f(1)) / DM
      if (n > 2) df_dmu(2:n-1) = (f(3:n) - f(1:n-2)) / (2.e0_wp * DM)
      df_dmu(n) = (f(n) - f(n-1)) / DM
      d2f_dmu2(1) = (f(3) - 2.e0_wp*f(2) + f(1)) / DM**2
      if (n > 3) d2f_dmu2(2:n-1) = (f(3:n) - 2.e0_wp*f(2:n-1) + f(1:n-2)) / DM**2
      d2f_dmu2(n) = (f(n) - 2.e0_wp*f(n-1) + f(n-2)) / DM**2
      return
    end if

    if (n < 7) then
      df_dmu(1)   = (-25.e0_wp*f(1)+48.e0_wp*f(2)-36.e0_wp*f(3)+16.e0_wp*f(4)-3.e0_wp*f(5)) / (12.e0_wp*DM)
      df_dmu(2)   = (-3.e0_wp*f(1)-10.e0_wp*f(2)+18.e0_wp*f(3)-6.e0_wp*f(4)+f(5)) / (12.e0_wp*DM)
      df_dmu(3:n-2) = (-f(5:n)+8.e0_wp*f(4:n-1)-8.e0_wp*f(2:n-3)+f(1:n-4)) / (12.e0_wp*DM)
      df_dmu(n-1) = (3.e0_wp*f(n)+10.e0_wp*f(n-1)-18.e0_wp*f(n-2)+6.e0_wp*f(n-3)-f(n-4)) / (12.e0_wp*DM)
      df_dmu(n)   = (25.e0_wp*f(n)-48.e0_wp*f(n-1)+36.e0_wp*f(n-2)-16.e0_wp*f(n-3)+3.e0_wp*f(n-4)) / (12.e0_wp*DM)
      d2f_dmu2(1)   = (35.e0_wp*f(1)-104.e0_wp*f(2)+114.e0_wp*f(3)-56.e0_wp*f(4)+11.e0_wp*f(5)) / (12.e0_wp*DM**2)
      d2f_dmu2(2)   = (11.e0_wp*f(1)-20.e0_wp*f(2)+6.e0_wp*f(3)+4.e0_wp*f(4)-f(5)) / (12.e0_wp*DM**2)
      d2f_dmu2(3:n-2) = (-f(5:n)+16.e0_wp*f(4:n-1)-30.e0_wp*f(3:n-2)+16.e0_wp*f(2:n-3)-f(1:n-4)) / (12.e0_wp*DM**2)
      d2f_dmu2(n-1) = (11.e0_wp*f(n)-20.e0_wp*f(n-1)+6.e0_wp*f(n-2)+4.e0_wp*f(n-3)-f(n-4)) / (12.e0_wp*DM**2)
      d2f_dmu2(n)   = (35.e0_wp*f(n)-104.e0_wp*f(n-1)+114.e0_wp*f(n-2)-56.e0_wp*f(n-3)+11.e0_wp*f(n-4)) / (12.e0_wp*DM**2)
      return
    end if

    inv60DM   = 1.e0_wp / (60.e0_wp * DM)
    inv180DM2 = 1.e0_wp / (180.e0_wp * DM**2)

    ! 6th-order one-sided boundaries (f')
    df_dmu(1) = (-147.e0_wp*f(1)+360.e0_wp*f(2)-450.e0_wp*f(3)+400.e0_wp*f(4) &
                 -225.e0_wp*f(5)+ 72.e0_wp*f(6)- 10.e0_wp*f(7)) * inv60DM
    df_dmu(2) = ( -10.e0_wp*f(1)- 77.e0_wp*f(2)+150.e0_wp*f(3)-100.e0_wp*f(4) &
                 +  50.e0_wp*f(5)- 15.e0_wp*f(6)+  2.e0_wp*f(7)) * inv60DM
    df_dmu(3) = (   2.e0_wp*f(1)- 24.e0_wp*f(2)- 35.e0_wp*f(3)+ 80.e0_wp*f(4) &
                 -  30.e0_wp*f(5)+  8.e0_wp*f(6)-       f(7)) * inv60DM
    df_dmu(n-2) = (       f(n-6)-  8.e0_wp*f(n-5)+ 30.e0_wp*f(n-4)- 80.e0_wp*f(n-3) &
                   + 35.e0_wp*f(n-2)+ 24.e0_wp*f(n-1)-  2.e0_wp*f(n)) * inv60DM
    df_dmu(n-1) = (  -2.e0_wp*f(n-6)+ 15.e0_wp*f(n-5)- 50.e0_wp*f(n-4)+100.e0_wp*f(n-3) &
                   -150.e0_wp*f(n-2)+ 77.e0_wp*f(n-1)+ 10.e0_wp*f(n)) * inv60DM
    df_dmu(n)   = (  10.e0_wp*f(n-6)- 72.e0_wp*f(n-5)+225.e0_wp*f(n-4)-400.e0_wp*f(n-3) &
                   + 450.e0_wp*f(n-2)-360.e0_wp*f(n-1)+147.e0_wp*f(n)) * inv60DM

    ! 5th-order one-sided boundaries (f'')
    d2f_dmu2(1) = ( 812.e0_wp*f(1)-3132.e0_wp*f(2)+5265.e0_wp*f(3)-5080.e0_wp*f(4) &
                  +2970.e0_wp*f(5)- 972.e0_wp*f(6)+ 137.e0_wp*f(7)) * inv180DM2
    d2f_dmu2(2) = ( 137.e0_wp*f(1)- 147.e0_wp*f(2)- 255.e0_wp*f(3)+ 470.e0_wp*f(4) &
                  - 285.e0_wp*f(5)+  93.e0_wp*f(6)-  13.e0_wp*f(7)) * inv180DM2
    d2f_dmu2(3) = ( -13.e0_wp*f(1)+ 228.e0_wp*f(2)- 420.e0_wp*f(3)+ 200.e0_wp*f(4) &
                  +  15.e0_wp*f(5)-  12.e0_wp*f(6)+   2.e0_wp*f(7)) * inv180DM2
    d2f_dmu2(n-2) = (   2.e0_wp*f(n-6)-  12.e0_wp*f(n-5)+  15.e0_wp*f(n-4)+200.e0_wp*f(n-3) &
                     - 420.e0_wp*f(n-2)+ 228.e0_wp*f(n-1)-  13.e0_wp*f(n)) * inv180DM2
    d2f_dmu2(n-1) = ( -13.e0_wp*f(n-6)+  93.e0_wp*f(n-5)- 285.e0_wp*f(n-4)+470.e0_wp*f(n-3) &
                     - 255.e0_wp*f(n-2)- 147.e0_wp*f(n-1)+ 137.e0_wp*f(n)) * inv180DM2
    d2f_dmu2(n)   = ( 137.e0_wp*f(n-6)- 972.e0_wp*f(n-5)+2970.e0_wp*f(n-4)-5080.e0_wp*f(n-3) &
                     +5265.e0_wp*f(n-2)-3132.e0_wp*f(n-1)+ 812.e0_wp*f(n)) * inv180DM2

    ! Fused 6th-order centered interior: f(i-3:i+3) loaded once for both stencils
    do i = 4, n-3
      df_dmu(i)   = (   -f(i-3) + 9.e0_wp*f(i-2) - 45.e0_wp*f(i-1) &
                    + 45.e0_wp*f(i+1) -  9.e0_wp*f(i+2) +       f(i+3)) * inv60DM
      d2f_dmu2(i) = (2.e0_wp*f(i-3) - 27.e0_wp*f(i-2) + 270.e0_wp*f(i-1) - 490.e0_wp*f(i) &
                    +270.e0_wp*f(i+1) - 27.e0_wp*f(i+2) +  2.e0_wp*f(i+3)) * inv180DM2
    end do
  end subroutine deriv_mu_and_mumu_1d

  ! -----------------------------------------------------------------------
  ! Scalar wrappers kept for backward compatibility (deriv_r, deriv_t, etc.)
  ! -----------------------------------------------------------------------

  pure function deriv_s_1d(f) result(df_ds)
    real(wp), intent(in) :: f(:)
    real(wp) :: df_ds(size(f))
    real(wp) :: d2f_ds2(size(f))
    call deriv_s_and_ss_1d(f, df_ds, d2f_ds2)
  end function deriv_s_1d

  pure function deriv_ss_1d(f) result(d2f_ds2)
    real(wp), intent(in) :: f(:)
    real(wp) :: d2f_ds2(size(f))
    real(wp) :: df_ds(size(f))
    call deriv_s_and_ss_1d(f, df_ds, d2f_ds2)
  end function deriv_ss_1d

  pure function deriv_mu_1d(f) result(df_dmu)
    real(wp), intent(in) :: f(:)
    real(wp) :: df_dmu(size(f))
    real(wp) :: d2f_dmu2(size(f))
    call deriv_mu_and_mumu_1d(f, df_dmu, d2f_dmu2)
  end function deriv_mu_1d

  pure function deriv_mumu_1d(f) result(d2f_dmu2)
    real(wp), intent(in) :: f(:)
    real(wp) :: d2f_dmu2(size(f))
    real(wp) :: df_dmu(size(f))
    call deriv_mu_and_mumu_1d(f, df_dmu, d2f_dmu2)
  end function deriv_mumu_1d

end module ope_eq_mod
