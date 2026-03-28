module spin_workspace
  use para_mod, only: wp, SDIV, MDIV, LMAX, s_pwr, &
                      s_gp, mu, sin_theta, &
                      P_2n, P1_2n_1, sin_2n_1_theta, &
                      angular_collocation, COLLOCATION_UNI, w_mu
  use nag_compat_mod, only: d01gaf
  implicit none

  real(wp) :: mphi_tran = 1.e-11_wp ! threshold to swift to massless solver
  real(wp) :: dif
  real(wp), allocatable, target :: e_gsm_cache(:,:), e_rsm_cache(:,:), e2alpha_r2_cache(:,:)
  real(wp), allocatable, target :: Acoup4_cache(:,:)
  real(wp), allocatable, target :: dg_s_cache(:,:), dg_m_cache(:,:), d2g_ss_cache(:,:), d2g_mm_cache(:,:)
  real(wp), allocatable, target :: dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
  real(wp), allocatable, target :: ds_s_cache(:,:), ds_m_cache(:,:)
  real(wp), allocatable, target :: mr_cache(:), besseli_cache(:,:), besselk_cache(:,:), wfac_cache(:)
  real(wp), allocatable :: s1_geom(:), s1_sq_geom(:), s1_one_minus_s_geom(:), s2_geom(:), sgp4_geom(:), m1_geom(:)
  real(wp), allocatable :: rad_inv_s1(:), rad_left_rho(:), rad_right_gama(:), rad_ratio_s(:), rad_ratio_g(:)
  real(wp), allocatable :: sin_theta_inv(:)
  real(wp), allocatable :: sgp_term_2d_cache(:,:), sin_theta_2d_cache(:,:), sgp_2d_cache(:,:)
  real(wp), allocatable :: radial_quad_weights(:)
  real(wp), allocatable :: angular_quad_weights(:)
  real(wp), allocatable :: weighted_even_basis(:,:), weighted_gama_basis(:,:), weighted_omega_basis(:,:)
  real(wp), allocatable :: recon_even_massive_basis(:,:), recon_gama_basis(:,:), recon_omega_basis(:,:)
  real(wp), allocatable, target :: scratch_exp_mhalf_gsm(:,:), scratch_exp_rsm_mhalf_gsm(:,:)
  real(wp), allocatable, target :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
  real(wp), allocatable, target :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
  real(wp), allocatable, target :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
  real(wp), allocatable, target :: target_rho(:,:), target_gama(:,:), target_ww(:,:), target_sphi(:,:)
  character(len=10) :: metric_method = 'Picard'
  character(len=10) :: scalar_method = 'Picard'

contains

  subroutine allocate_workspace
    real(wp), allocatable :: sgp_term_1d(:)
    integer :: n
    allocate(e_gsm_cache(SDIV,MDIV), e_rsm_cache(SDIV,MDIV), e2alpha_r2_cache(SDIV,MDIV))
    allocate(Acoup4_cache(SDIV,MDIV))
    allocate(dg_s_cache(SDIV,MDIV), dg_m_cache(SDIV,MDIV), d2g_ss_cache(SDIV,MDIV), d2g_mm_cache(SDIV,MDIV))
    allocate(dr_s_cache(SDIV,MDIV), dr_m_cache(SDIV,MDIV), dww_s_cache(SDIV,MDIV), dww_m_cache(SDIV,MDIV))
    allocate(ds_s_cache(SDIV,MDIV), ds_m_cache(SDIV,MDIV))
    allocate(mr_cache(SDIV), besseli_cache(LMAX+1,SDIV), besselk_cache(LMAX+1,SDIV), wfac_cache(SDIV))
    allocate(s1_geom(SDIV), s1_sq_geom(SDIV), s1_one_minus_s_geom(SDIV), s2_geom(SDIV), sgp4_geom(SDIV), m1_geom(MDIV))
    allocate(rad_inv_s1(SDIV), rad_left_rho(SDIV), rad_right_gama(SDIV), rad_ratio_s(SDIV), rad_ratio_g(SDIV))
    allocate(sin_theta_inv(MDIV))
    allocate(sgp_term_2d_cache(SDIV,MDIV), sin_theta_2d_cache(SDIV,MDIV), sgp_2d_cache(SDIV,MDIV))
    allocate(radial_quad_weights(SDIV), angular_quad_weights(MDIV))
    allocate(weighted_even_basis(MDIV,LMAX+1))
    allocate(weighted_gama_basis(MDIV,LMAX))
    allocate(weighted_omega_basis(MDIV,LMAX))
    allocate(recon_even_massive_basis(MDIV,LMAX+1))
    allocate(recon_gama_basis(MDIV,LMAX))
    allocate(recon_omega_basis(MDIV,LMAX))
    allocate(scratch_exp_mhalf_gsm(SDIV,MDIV), scratch_exp_rsm_mhalf_gsm(SDIV,MDIV))
    allocate(target_rho(SDIV,MDIV), target_gama(SDIV,MDIV), target_ww(SDIV,MDIV), target_sphi(SDIV,MDIV))
    allocate(S_metric_rho(SDIV,MDIV), S_metric_gama(SDIV,MDIV), S_metric_omega(SDIV,MDIV), S_metric_sphi(SDIV,MDIV))
    allocate(D1_metric_rho(LMAX+1,SDIV), D1_metric_gama(LMAX+1,SDIV), D1_metric_omega(LMAX+1,SDIV), D1_metric_sphi(LMAX+1,SDIV))
    allocate(D2_metric_rho(SDIV,LMAX+1), D2_metric_gama(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1), D2_metric_sphi(SDIV,LMAX+1))

    call compute_effective_d01gaf_weights(s_gp, radial_quad_weights)
    if (angular_collocation /= COLLOCATION_UNI) then
      angular_quad_weights = w_mu
    else
      call compute_effective_d01gaf_weights(mu, angular_quad_weights)
    end if
    wfac_cache = 1.e0_wp / (1.e0_wp - s_gp)**2
    s1_geom = s_gp * (1.e0_wp - s_gp)
    s1_sq_geom = s1_geom**2
    s1_one_minus_s_geom = s1_geom * (1.e0_wp - s_gp)
    s2_geom = (s_gp / (1.e0_wp - s_gp))**2
    sgp4_geom = s_gp**4
    m1_geom = 1.e0_wp - mu**2
    rad_inv_s1(2:) = dble(s_pwr) / (s_gp(2:) * (1.e0_wp - s_gp(2:)))
    rad_left_rho(2:) = dble(s_pwr) * s_gp(2:)**(s_pwr - 1) / (1.e0_wp - s_gp(2:))**(s_pwr + 1)
    rad_right_gama(2:) = dble(s_pwr) * s_gp(2:)**(2 * s_pwr - 1) / (1.e0_wp - s_gp(2:))**(2 * s_pwr + 1)
    rad_ratio_s(2:) = ((1.e0_wp - s_gp(2:)) / s_gp(2:))**s_pwr
    rad_ratio_g(2:) = ((1.e0_wp - s_gp(2:)) / s_gp(2:))**(2 * s_pwr)
    sin_theta_inv = 0.e0_wp
    sin_theta_inv(2:) = 1.e0_wp / sin_theta(2:)
    allocate(sgp_term_1d(SDIV))
    sgp_term_1d = s_gp / (1.e0_wp - s_gp)
    sgp_term_2d_cache = spread(sgp_term_1d, 2, MDIV)
    sin_theta_2d_cache = spread(sin_theta, 1, SDIV)
    sgp_2d_cache = spread(s_gp, 2, MDIV)
    deallocate(sgp_term_1d)
    weighted_even_basis = spread(angular_quad_weights, 2, LMAX+1) * P_2n
    do n = 0, LMAX
      recon_even_massive_basis(:,n+1) = real(2*n+1, wp) * P_2n(:,n+1)
    end do
    if (LMAX > 0) then
      weighted_gama_basis = spread(angular_quad_weights, 2, LMAX) * sin_2n_1_theta
      weighted_omega_basis = spread(angular_quad_weights * sin_theta, 2, LMAX) * P1_2n_1(:,2:LMAX+1)
      recon_gama_basis = 0.e0_wp
      recon_gama_basis(:,1) = 1.e0_wp
      do n = 2, LMAX
        recon_gama_basis(1:MDIV-1,n) = sin_2n_1_theta(1:MDIV-1,n) * sin_theta_inv(1:MDIV-1) / real(2*n-1, wp)
        recon_gama_basis(MDIV,n) = 1.e0_wp
      end do
      do n = 1, LMAX
        recon_omega_basis(1:MDIV-1,n) = -P1_2n_1(1:MDIV-1,n+1) * sin_theta_inv(1:MDIV-1) / real(2*n*(2*n-1), wp)
        recon_omega_basis(MDIV,n) = 0.5_wp
      end do
    end if
  end subroutine allocate_workspace

  subroutine deallocate_workspace
    deallocate(e_gsm_cache, e_rsm_cache, e2alpha_r2_cache)
    deallocate(Acoup4_cache)
    deallocate(dg_s_cache, dg_m_cache, d2g_ss_cache, d2g_mm_cache)
    deallocate(dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache)
    deallocate(ds_s_cache, ds_m_cache)
    deallocate(mr_cache, besseli_cache, besselk_cache, wfac_cache)
    deallocate(s1_geom, s1_sq_geom, s1_one_minus_s_geom, s2_geom, sgp4_geom, m1_geom)
    deallocate(rad_inv_s1, rad_left_rho, rad_right_gama, rad_ratio_s, rad_ratio_g)
    deallocate(sin_theta_inv)
    deallocate(sgp_term_2d_cache, sin_theta_2d_cache, sgp_2d_cache)
    deallocate(weighted_even_basis, weighted_gama_basis, weighted_omega_basis)
    deallocate(recon_even_massive_basis, recon_gama_basis, recon_omega_basis)
    deallocate(scratch_exp_mhalf_gsm, scratch_exp_rsm_mhalf_gsm)
    deallocate(target_rho, target_gama, target_ww, target_sphi)
    deallocate(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi)
    deallocate(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    deallocate(D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    if (allocated(radial_quad_weights)) deallocate(radial_quad_weights)
    if (allocated(angular_quad_weights)) deallocate(angular_quad_weights)
  end subroutine deallocate_workspace

  subroutine compute_effective_d01gaf_weights(x, w)
    real(wp), intent(in) :: x(:)
    real(wp), intent(out) :: w(:)
    real(wp), allocatable :: basis(:)
    real(wp) :: ans, er
    integer :: i, ifail, n

    n = size(x)
    if (size(w) /= n) stop "compute_effective_d01gaf_weights: size mismatch"

    allocate(basis(n))
    basis = 0.e0_wp

    do i = 1, n
      basis(i) = 1.e0_wp
      call d01gaf(x, basis, n, ans, er, ifail)
      if (ifail /= 0) stop "compute_effective_d01gaf_weights: d01gaf failed"
      w(i) = ans
      basis(i) = 0.e0_wp
    end do

    deallocate(basis)
  end subroutine compute_effective_d01gaf_weights

end module spin_workspace
