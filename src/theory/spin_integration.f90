module spin_integration
  use para_mod, only: wp, SDIV, MDIV, LMAX, &
                      s_gp, mu, sin_theta, s_pwr, DM, &
                      rho, gama, alpha, ww, omg, sphi, &
                      energy, pressure, enthalpy, velocity_sq, &
                      P_2n, P1_2n_1, &
                      r_ratio, r_e, enthalpy_min, &
                      Omega_c, Omega_e, M2, S3, M4, sphi_m, &
                      B_coup, mphi_r, &
                      l_uni, KAPPA, C, G, MSUN, MB, pi, KSCALE, &
                      mass, mass_0, ang_mom, &
                      A_diff, lambda1, lambda2, &
                      solver_type, output, timing, eos_file
  use eos_mod, only: n0_at_e
  use toolkit_mod, only: besseli, besselk
  use spin_derivatives, only: deriv_s_sub, deriv_m_sub
  use spin_workspace
  implicit none
  private
  public :: get_all_targets, update_alpha_potential, output_helper

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(lda,*), b(ldb,*)
      double precision, intent(inout) :: c(ldc,*)
    end subroutine dgemm
    subroutine dpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, kd, nrhs, ldab, ldb
      double precision, intent(inout) :: ab(ldab,*), b(ldb,*)
      integer, intent(out) :: info
    end subroutine dpbsv
  end interface

contains

  subroutine precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, besseli_cache, besselk_cache, &
                                                dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                                e2alpha_r2_cache, Acoup4_cache)
    real(wp), intent(in)  :: r_e_new
    real(wp), intent(in)  :: root_mphi_re
    real(wp), intent(out) :: mr_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    real(wp), intent(out) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(out) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(out) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(wp), intent(out) :: e2alpha_r2_cache(:,:), Acoup4_cache(:,:)
    integer :: s, m, n
    logical :: use_massive_scalar
    real(wp) :: s1(SDIV), m1(MDIV)
    real(wp) :: one_minus_2s(SDIV)

    use_massive_scalar = mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp
    if (use_massive_scalar) then
      mr_cache(:) = root_mphi_re * s_gp(:) / (1.e0_wp - s_gp(:))
      do n = 0, LMAX
        do s = 1, SDIV
          besseli_cache(n+1,s) = besseli(2*n, mr_cache(s))
          besselk_cache(n+1,s) = besselk(2*n, mr_cache(s))
        end do
      end do
    end if

    call deriv_s_sub(gama, dg_s_cache)
    call deriv_s_sub(rho, dr_s_cache)
    call deriv_s_sub(ww, dww_s_cache)
    call deriv_s_sub(sphi, ds_s_cache)
    if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
      dg_m_cache = 0.e0_wp
      dr_m_cache = 0.e0_wp
      dww_m_cache = 0.e0_wp
      ds_m_cache = 0.e0_wp
    else
      call deriv_m_sub(gama, dg_m_cache)
      call deriv_m_sub(rho, dr_m_cache)
      call deriv_m_sub(ww, dww_m_cache)
      call deriv_m_sub(sphi, ds_m_cache)
    end if

    s1 = s_gp * (1.e0_wp - s_gp)
    m1 = 1.e0_wp - mu**2
    one_minus_2s = 1.e0_wp - 2.e0_wp * s_gp
    call deriv_s_sub(dg_s_cache, d2g_ss_cache)
    do m = 1, MDIV
      d2g_ss_cache(:,m) = s1 * d2g_ss_cache(:,m) + one_minus_2s * dg_s_cache(:,m)
    end do
    if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
      d2g_mm_cache = 0.e0_wp
    else
      call deriv_m_sub(dg_m_cache, d2g_mm_cache)
      do m = 1, MDIV
        d2g_mm_cache(:,m) = m1(m) * d2g_mm_cache(:,m) - 2.e0_wp * mu(m) * dg_m_cache(:,m)
      end do
    end if
    e_gsm_cache      = exp(0.5e0_wp * gama)
    e_rsm_cache      = exp(-rho)
    e2alpha_r2_cache = exp(2.e0_wp * alpha) * r_e_new**2
    Acoup4_cache     = exp(-sphi**2 * B_coup)
  end subroutine precompute_derivatives_and_bessels

  subroutine build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, &
                                ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, &
                                e2alpha_r2_cache, Acoup4_cache)
    real(wp), intent(in)  :: r_e_new
    real(wp), intent(out) :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(wp), intent(in)  :: dr_s_cache(:,:), dr_m_cache(:,:), dg_s_cache(:,:), dg_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(in)  :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(in)  :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_gsm_cache(:,:), e_rsm_cache(:,:)
    real(wp), intent(in)  :: e2alpha_r2_cache(:,:), Acoup4_cache(:,:)
    integer :: m
    real(wp) :: mum, m1
    real(wp), dimension(SDIV) :: esm_col, psm_col, vphi_col, scal_p_col
    real(wp), dimension(SDIV) :: vel_fac_col, matter_sum_col, matter_trace_col
    real(wp), dimension(SDIV) :: e2alpha_s2_col, source_common_col, e_rsm2_col
    real(wp), dimension(SDIV) :: dg_s_scaled_col, dg_m_scaled_col, rho_bracket_col
    real(wp), dimension(SDIV) :: omega_matter_col, omega_bracket_col, sphi_source_col
    real(wp), dimension(SDIV) :: vsq_col, one_plus_vsq_col

    do m = 1, MDIV
      mum = mu(m)
      m1 = m1_geom(m)

      associate( gs => dg_s_cache(:,m), gm => dg_m_cache(:,m), &
                 rs => dr_s_cache(:,m), rm => dr_m_cache(:,m), &
                 wws => dww_s_cache(:,m), wwm => dww_m_cache(:,m), &
                 ss => ds_s_cache(:,m), sm => ds_m_cache(:,m), &
                 gss => d2g_ss_cache(:,m), gmm => d2g_mm_cache(:,m), &
                 egsm => e_gsm_cache(:,m), ersm => e_rsm_cache(:,m), &
                 e2ar2 => e2alpha_r2_cache(:,m), Ac4 => Acoup4_cache(:,m), &
                 ww_col => ww(:,m), omg_col => omg(:,m), &
                 rho_col => rho(:,m), gama_col => gama(:,m), sphi_col => sphi(:,m) )

      esm_col = energy(:,m) * Ac4
      psm_col = pressure(:,m) * Ac4
      vphi_col = sphi_col**2 * mphi_r * 0.5e0_wp * e2ar2
      scal_p_col = sphi_col * egsm
      vsq_col = velocity_sq(:,m)
      one_plus_vsq_col = 1.e0_wp + vsq_col
      vel_fac_col = 1.e0_wp / (1.e0_wp - vsq_col)
      matter_sum_col = esm_col + psm_col
      matter_trace_col = esm_col - 3.e0_wp * psm_col
      e_rsm2_col = ersm**2
      e2alpha_s2_col = e2ar2 * s2_geom
      source_common_col = 16.e0_wp * pi * e2alpha_s2_col * psm_col - 4.e0_wp * vphi_col * s2_geom
      dg_s_scaled_col = s1_geom * gs
      dg_m_scaled_col = m1 * gm

      rho_bracket_col = source_common_col &
        - dg_s_scaled_col * (0.5e0_wp * dg_s_scaled_col + 1.e0_wp) &
        - gm * (0.5e0_wp * dg_m_scaled_col - mum)

      S_metric_rho(:,m) = egsm * ( &
          8.e0_wp * pi * e2alpha_s2_col * matter_sum_col * one_plus_vsq_col * vel_fac_col &
        + s2_geom * m1 * e_rsm2_col * ((s1_geom * wws)**2 + m1 * wwm**2) &
        + dg_s_scaled_col - mum * gm &
        + rho_col * 0.5e0_wp * rho_bracket_col )

      S_metric_gama(:,m) = egsm * (source_common_col &
        + gama_col * 0.5e0_wp * (source_common_col - 0.5e0_wp * dg_s_scaled_col**2 &
        - 0.5e0_wp * dg_m_scaled_col * gm) )

      if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
        S_metric_omega(:,m) = 0.e0_wp
      else
        omega_matter_col = (one_plus_vsq_col * esm_col + 2.e0_wp * vsq_col * psm_col) * vel_fac_col
        omega_bracket_col = -8.e0_wp * pi * e2alpha_s2_col * omega_matter_col &
          - s1_geom * (2.e0_wp * rs + 0.5e0_wp * gs) &
          + mum * (2.e0_wp * rm + 0.5e0_wp * gm) &
          + 0.25e0_wp * s1_sq_geom * (4.e0_wp * rs**2 - gs**2) &
          + 0.25e0_wp * m1 * (4.e0_wp * rm**2 - gm**2) &
          - m1 * e_rsm2_col * (sgp4_geom * wws**2 + s2_geom * m1 * wwm**2) &
          - 2.e0_wp * vphi_col * s2_geom

        S_metric_omega(:,m) = egsm * ersm * ( &
          -16.e0_wp * pi * e2alpha_s2_col * (omg_col - ww_col) * matter_sum_col * vel_fac_col &
          + ww_col * omega_bracket_col )
      endif

      if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp ) then
        sphi_source_col = -2.e0_wp * pi * B_coup * matter_trace_col + mphi_r
        S_metric_sphi(:,m) = -r_e_new**2 * s2_geom * scal_p_col * mphi_r &
          + e2alpha_s2_col * scal_p_col * sphi_source_col &
          + scal_p_col * (s1_one_minus_s_geom * gs + s1_sq_geom * (0.5e0_wp * gss + 0.25e0_wp * gs**2) &
          + m1 * (0.5e0_wp * gmm + 0.25e0_wp * gm**2) - mum * gm)
      else
        sphi_source_col = -2.e0_wp * pi * B_coup * matter_trace_col
        S_metric_sphi(:,m) = -s1_sq_geom * gs * ss - dg_m_scaled_col * sm &
          + sphi_col * sphi_source_col * e2alpha_s2_col
      endif
      end associate
    end do
  end subroutine build_source_terms

  subroutine angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
                                D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    real(wp), intent(in)  :: S_metric_rho(:,:), S_metric_gama(:,:), S_metric_omega(:,:), S_metric_sphi(:,:)
    real(wp), intent(out) :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    call dgemm('T', 'T', LMAX+1, SDIV, MDIV, 1.e0_wp, weighted_even_basis, MDIV, S_metric_rho, SDIV, 0.e0_wp, D1_metric_rho, LMAX+1)
    call dgemm('T', 'T', LMAX+1, SDIV, MDIV, 1.e0_wp, weighted_even_basis, MDIV, S_metric_sphi, SDIV, 0.e0_wp, D1_metric_sphi, LMAX+1)
    D1_metric_gama = 0.e0_wp
    D1_metric_omega = 0.e0_wp
    if (LMAX > 0) then
      block
        real(wp) :: tmp_gama(LMAX, SDIV), tmp_omega(LMAX, SDIV)
        call dgemm('T', 'T', LMAX, SDIV, MDIV, 1.e0_wp, weighted_gama_basis, MDIV, S_metric_gama, SDIV, 0.e0_wp, tmp_gama, LMAX)
        call dgemm('T', 'T', LMAX, SDIV, MDIV, 1.e0_wp, weighted_omega_basis, MDIV, S_metric_omega, SDIV, 0.e0_wp, tmp_omega, LMAX)
        D1_metric_gama(2:LMAX+1,:) = tmp_gama
        D1_metric_omega(2:LMAX+1,:) = tmp_omega
      end block
    end if
  end subroutine angular_integration

  subroutine radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
                                D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, &
                                root_mphi_re, wfac_cache, besseli_cache, besselk_cache)
    real(wp), intent(in)  :: D1_metric_rho(:,:), D1_metric_gama(:,:), D1_metric_omega(:,:), D1_metric_sphi(:,:)
    real(wp), intent(out) :: D2_metric_rho(:,:), D2_metric_gama(:,:), D2_metric_omega(:,:), D2_metric_sphi(:,:)
    real(wp), intent(in)  :: root_mphi_re
    real(wp), intent(in)  :: wfac_cache(:), besseli_cache(:,:), besselk_cache(:,:)
    integer :: n
    real(wp) :: weighted_source(SDIV)

    D2_metric_gama (:,1) = 0.e0_wp
    D2_metric_omega(:,1) = 0.e0_wp

    do n = 0, LMAX
      weighted_source = radial_quad_weights * D1_metric_rho(n+1,:)
      call integrate_rho_like(n, weighted_source, D2_metric_rho(:,n+1))

      if (mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp) then
        weighted_source = radial_quad_weights * wfac_cache * D1_metric_sphi(n+1,:) * root_mphi_re
        call integrate_massive_sphi(n+1, weighted_source, D2_metric_sphi(:,n+1))
      else
        weighted_source = radial_quad_weights * D1_metric_sphi(n+1,:)
        call integrate_rho_like(n, weighted_source, D2_metric_sphi(:,n+1))
      end if
    end do

    do n = 1, LMAX
      weighted_source = radial_quad_weights * D1_metric_gama(n+1,:)
      call integrate_gama(n, weighted_source, D2_metric_gama(:,n+1))

      weighted_source = radial_quad_weights * D1_metric_omega(n+1,:)
      call integrate_omega(n, weighted_source, D2_metric_omega(:,n+1))
    end do

  contains
    subroutine integrate_rho_like(n_phys, source_weights, out_values)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: j, k
      real(wp) :: prefactor, ratio_j
      real(wp) :: f2n_vals(SDIV), left_terms(SDIV), right_terms(SDIV)
      real(wp) :: left_prefix(SDIV), right_suffix(SDIV)

      prefactor = dble(s_pwr)
      f2n_vals = 1.e0_wp
      left_terms = 0.e0_wp
      right_terms = 0.e0_wp
      left_prefix = 0.e0_wp
      right_suffix = 0.e0_wp

      if (n_phys == 0) then
        do k = 2, SDIV
          left_terms(k) = prefactor * s_gp(k)**(s_pwr - 1) * source_weights(k) / (1.e0_wp - s_gp(k))**(s_pwr + 1)
          right_terms(k) = prefactor * source_weights(k) / (s_gp(k) * (1.e0_wp - s_gp(k)))
        end do
      else
        do k = 2, SDIV
          f2n_vals(k) = ((1.e0_wp - s_gp(k)) / s_gp(k))**(2 * s_pwr * n_phys)
          left_terms(k) = prefactor * s_gp(k)**(s_pwr - 1) * source_weights(k) &
                        / (f2n_vals(k) * (1.e0_wp - s_gp(k))**(s_pwr + 1))
          right_terms(k) = prefactor * f2n_vals(k) * source_weights(k) / (s_gp(k) * (1.e0_wp - s_gp(k)))
        end do
      end if

      do k = 2, SDIV
        left_prefix(k) = left_prefix(k-1) + left_terms(k)
      end do
      right_suffix(SDIV) = right_terms(SDIV)
      do k = SDIV-1, 1, -1
        right_suffix(k) = right_suffix(k+1) + right_terms(k)
      end do

      out_values(1) = merge(right_suffix(1), 0.e0_wp, n_phys == 0)
      do j = 2, SDIV
        ratio_j = ((1.e0_wp - s_gp(j)) / s_gp(j))**s_pwr
        if (n_phys == 0) then
          out_values(j) = ratio_j * left_prefix(j-1) + right_suffix(j)
        else
          out_values(j) = f2n_vals(j) * ratio_j * left_prefix(j-1) + right_suffix(j) / f2n_vals(j)
        end if
      end do
    end subroutine integrate_rho_like

    subroutine integrate_gama(n_phys, source_weights, out_values)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: j, k
      real(wp) :: prefactor, ratio_j
      real(wp) :: f2n_vals(SDIV), left_terms(SDIV), right_terms(SDIV)
      real(wp) :: left_prefix(SDIV), right_suffix(SDIV), boundary_sum

      prefactor = dble(s_pwr)
      f2n_vals = 1.e0_wp
      left_terms = 0.e0_wp
      right_terms = 0.e0_wp
      left_prefix = 0.e0_wp
      right_suffix = 0.e0_wp
      boundary_sum = 0.e0_wp

      do k = 2, SDIV
        f2n_vals(k) = ((1.e0_wp - s_gp(k)) / s_gp(k))**(2 * s_pwr * n_phys)
        left_terms(k) = prefactor * source_weights(k) / (f2n_vals(k) * s_gp(k) * (1.e0_wp - s_gp(k)))
        right_terms(k) = prefactor * f2n_vals(k) * s_gp(k)**(2 * s_pwr - 1) * source_weights(k) &
                      / (1.e0_wp - s_gp(k))**(2 * s_pwr + 1)
        boundary_sum = boundary_sum + prefactor * source_weights(k) / (s_gp(k) * (1.e0_wp - s_gp(k)))
      end do

      do k = 2, SDIV
        left_prefix(k) = left_prefix(k-1) + left_terms(k)
      end do
      right_suffix(SDIV) = right_terms(SDIV)
      do k = SDIV-1, 1, -1
        right_suffix(k) = right_suffix(k+1) + right_terms(k)
      end do

      out_values(1) = merge(boundary_sum, 0.e0_wp, n_phys == 1)
      do j = 2, SDIV
        ratio_j = ((1.e0_wp - s_gp(j)) / s_gp(j))**(2 * s_pwr)
        out_values(j) = prefactor * f2n_vals(j) * left_prefix(j-1) + ratio_j * right_suffix(j) / f2n_vals(j)
      end do
    end subroutine integrate_gama

    subroutine integrate_omega(n_phys, source_weights, out_values)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: j, k
      real(wp) :: prefactor, ratio_rho_j, ratio_gama_j
      real(wp) :: f2n_vals(SDIV), left_terms(SDIV), right_terms(SDIV)
      real(wp) :: left_prefix(SDIV), right_suffix(SDIV), boundary_sum

      prefactor = dble(s_pwr)
      f2n_vals = 1.e0_wp
      left_terms = 0.e0_wp
      right_terms = 0.e0_wp
      left_prefix = 0.e0_wp
      right_suffix = 0.e0_wp
      boundary_sum = 0.e0_wp

      do k = 2, SDIV
        f2n_vals(k) = ((1.e0_wp - s_gp(k)) / s_gp(k))**(2 * s_pwr * n_phys)
        left_terms(k) = prefactor * s_gp(k)**(s_pwr - 1) * source_weights(k) &
                      / (f2n_vals(k) * (1.e0_wp - s_gp(k))**(s_pwr + 1))
        right_terms(k) = prefactor * f2n_vals(k) * s_gp(k)**(2 * s_pwr - 1) * source_weights(k) &
                      / (1.e0_wp - s_gp(k))**(2 * s_pwr + 1)
        boundary_sum = boundary_sum + prefactor * source_weights(k) / (s_gp(k) * (1.e0_wp - s_gp(k)))
      end do

      do k = 2, SDIV
        left_prefix(k) = left_prefix(k-1) + left_terms(k)
      end do
      right_suffix(SDIV) = right_terms(SDIV)
      do k = SDIV-1, 1, -1
        right_suffix(k) = right_suffix(k+1) + right_terms(k)
      end do

      out_values(1) = merge(boundary_sum, 0.e0_wp, n_phys == 1)
      do j = 2, SDIV
        ratio_rho_j = ((1.e0_wp - s_gp(j)) / s_gp(j))**s_pwr
        ratio_gama_j = ((1.e0_wp - s_gp(j)) / s_gp(j))**(2 * s_pwr)
        out_values(j) = f2n_vals(j) * ratio_rho_j * left_prefix(j-1) + ratio_gama_j * right_suffix(j) / f2n_vals(j)
      end do
    end subroutine integrate_omega

    subroutine integrate_massive_sphi(n_idx, source_weights, out_values)
      integer, intent(in) :: n_idx
      real(wp), intent(in) :: source_weights(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      real(wp) :: left_prefix(SDIV), right_suffix(SDIV)
      real(wp) :: source_left(SDIV), source_right(SDIV)
      integer :: s

      source_left = source_weights * besseli_cache(n_idx,:)
      source_right = source_weights * besselk_cache(n_idx,:)

      left_prefix(1) = source_left(1)
      do s = 2, SDIV
        left_prefix(s) = left_prefix(s-1) + source_left(s)
      end do

      right_suffix(SDIV) = source_right(SDIV)
      do s = SDIV-1, 1, -1
        right_suffix(s) = right_suffix(s+1) + source_right(s)
      end do

      out_values(1) = besseli_cache(n_idx,1) * right_suffix(1)
      do s = 2, SDIV
        out_values(s) = besselk_cache(n_idx,s) * left_prefix(s-1) + &
                        besseli_cache(n_idx,s) * right_suffix(s)
      end do
    end subroutine integrate_massive_sphi
  end subroutine radial_integration

  subroutine sum_coefficients_and_get_targets(out_target_rho, out_target_gama, out_target_ww, out_target_sphi, &
                                              D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    real(wp), intent(out) :: out_target_rho(:,:), out_target_gama(:,:), out_target_ww(:,:), out_target_sphi(:,:)
    real(wp), intent(in)  :: D2_metric_rho(SDIV,LMAX+1), D2_metric_gama(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1), D2_metric_sphi(SDIV,LMAX+1)
    integer :: n, m
    real(wp), dimension(SDIV,MDIV) :: exp_mhalf_gsm, exp_rsm_mhalf_gsm
    real(wp), dimension(SDIV,MDIV) :: sum_rho, sum_sphi, sum_gama, sum_omega
    real(wp), dimension(MDIV) :: sin_theta_inv

    exp_mhalf_gsm = 1.0e0_wp / e_gsm_cache
    exp_rsm_mhalf_gsm = exp_mhalf_gsm / e_rsm_cache
    sin_theta_inv = 0.e0_wp; where (sin_theta > 1.e-12_wp) sin_theta_inv = 1.e0_wp / sin_theta

    call dgemm('N','T', SDIV, MDIV, LMAX+1, 1.e0_wp, D2_metric_rho, SDIV, P_2n, MDIV, 0.e0_wp, sum_rho, SDIV)
    sum_rho = -exp_mhalf_gsm * sum_rho

    sum_gama = 0.e0_wp
    sum_gama(:,1:MDIV-1) = -(2.e0_wp/pi) * exp_mhalf_gsm(:,1:MDIV-1) * spread(D2_metric_gama(:,2), 2, MDIV-1)
    sum_gama(:,MDIV) = -(2.e0_wp/pi) * exp_mhalf_gsm(:, MDIV) * D2_metric_gama(:, 2)
    sum_omega = 0.e0_wp

    if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp ) then
      block
        real(wp) :: D2_sphi_scaled(SDIV, LMAX+1)
        integer :: nn
        do nn = 0, LMAX
          D2_sphi_scaled(:,nn+1) = real(2*nn+1, wp) * D2_metric_sphi(:,nn+1)
        end do
        call dgemm('N','T', SDIV, MDIV, LMAX+1, 1.e0_wp, D2_sphi_scaled, SDIV, P_2n, MDIV, 0.e0_wp, sum_sphi, SDIV)
      end block
      sum_sphi = -exp_mhalf_gsm * sum_sphi
    else
      call dgemm('N','T', SDIV, MDIV, LMAX+1, -1.e0_wp, D2_metric_sphi, SDIV, P_2n, MDIV, 0.e0_wp, sum_sphi, SDIV)
    endif

    if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
      out_target_rho  = sum_rho
      out_target_gama = sum_gama
      out_target_ww   = sum_omega
      out_target_sphi = sum_sphi
      return
    endif
    do n = 1, LMAX
      do m = 1, MDIV-1
        sum_omega(:,m) = sum_omega(:,m) - exp_rsm_mhalf_gsm(:,m) * D2_metric_omega(:,n+1) * &
          (P1_2n_1(m,n+1) * sin_theta_inv(m) / (2.e0_wp*n*(2.e0_wp*n-1.e0_wp)))
      end do
      sum_omega(:,MDIV) = sum_omega(:,MDIV) + exp_rsm_mhalf_gsm(:,MDIV) * D2_metric_omega(:,n+1) / 2.e0_wp
    end do

    do n = 2, LMAX
      do m = 1, MDIV-1
        sum_gama(:,m) = sum_gama(:,m) - (2.e0_wp/pi) * exp_mhalf_gsm(:,m) * D2_metric_gama(:,n+1) * &
          (sin_2n_1_theta(m,n) * sin_theta_inv(m) / (2.e0_wp*n-1.e0_wp))
      end do
      sum_gama(:,MDIV) = sum_gama(:,MDIV) - (2.e0_wp/pi) * exp_mhalf_gsm(:,MDIV) * D2_metric_gama(:,n+1)
    end do

    out_target_rho  = sum_rho
    out_target_gama = sum_gama
    out_target_ww   = sum_omega
    out_target_sphi = sum_sphi
  end subroutine sum_coefficients_and_get_targets

  subroutine update_alpha_potential(r_e_new, dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                    ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_rsm_cache)
    real(wp), intent(in) :: r_e_new
    real(wp), intent(in) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(in) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(in) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_rsm_cache(:,:)
    integer :: m
    real(wp) :: m1, mu_m
    real(wp), dimension(SDIV,MDIV) :: da_dm, d_gama_sm_all
    real(wp), dimension(SDIV) :: sgp_ratio_cache
    real(wp), dimension(SDIV) :: temp1_col, temp2_col, temp3_col, temp4_col, temp5_col, temp6_col, temp7_col, temp8_col, temp9_col
    real(wp), dimension(SDIV) :: numer_m, one_plus_s1dgs, da_col
    real(wp) :: adj_const(SDIV)

    alpha(:,:) = 0.e0_wp
    if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
      return
    else
      sgp_ratio_cache = s_gp / (1.e0_wp-s_gp)

      da_dm(1,:) = 0.0e0_wp
      call deriv_m_sub(dg_s_cache, d_gama_sm_all)
      do m = 1, MDIV
        mu_m = mu(m)
        m1   = 1.e0_wp - mu_m**2
        associate( gs => dg_s_cache(:,m), gm => dg_m_cache(:,m), &
             rs => dr_s_cache(:,m), rm => dr_m_cache(:,m), &
             ss => ds_s_cache(:,m), sm => ds_m_cache(:,m), &
             gsm => d_gama_sm_all(:,m), wws => dww_s_cache(:,m), &
             wwm => dww_m_cache(:,m), gss => d2g_ss_cache(:,m), &
             gmm => d2g_mm_cache(:,m), e_cache => e_rsm_cache(:,m) )

        numer_m        = -mu_m + m1 * gm
        one_plus_s1dgs =  1.e0_wp + s1_geom * gs

        temp1_col = 2.e0_wp * s_gp**2 * sgp_ratio_cache * m1 * wws * wwm * one_plus_s1dgs &
          - ( (s_gp**2 * wws)**2 - (s_gp * wwm * sgp_ratio_cache)**2 * m1 ) * numer_m
        temp2_col = 1.e0_wp / ( m1 * one_plus_s1dgs**2 + numer_m**2 )
        temp3_col = s1_geom * gss + (s1_geom * gs)**2
        temp4_col = gm * numer_m
        temp5_col = ( (s1_geom * (rs + gs))**2 - m1 * (rm + gm)**2 ) * numer_m
        temp6_col = s1_geom * m1 * (  (rs + gs) * (rm + gm) / 2.e0_wp + gsm + gs * gm  ) * one_plus_s1dgs
        temp7_col = s1_geom * mu_m * gs * one_plus_s1dgs
        temp8_col = m1 * e_cache * e_cache
        temp9_col = -temp2_col * numer_m * ( (s1_geom * ss)**2 - m1 * sm**2 ) &
              - m1 * s1_geom * one_plus_s1dgs * 2.e0_wp * sm * ss

        da_col = - (rm + gm) / 2.e0_wp &
          - temp2_col * ( (temp3_col - gmm - temp4_col) * numer_m / 2.e0_wp &
          + temp5_col / 4.e0_wp - temp6_col  + temp7_col + temp8_col * temp1_col / 4.e0_wp ) + temp9_col
        da_dm(2:SDIV,m) = da_col(2:SDIV)
        end associate
      end do

      do m = 1, MDIV-1
        alpha(:,m+1) = alpha(:,m) + dm * ( da_dm(:,m+1) + da_dm(:,m) ) * 0.5e0_wp
      enddo

      alpha(SDIV,:) = 0.e0_wp
      adj_const = alpha(:,MDIV) - ( gama(:,MDIV) - rho(:,MDIV) )/2.e0_wp
      alpha = alpha - spread(adj_const, DIM=2, NCOPIES=MDIV)
    end if
    if (any(alpha .ge. 300.0)) then
      write(*,*) "Error: Alpha fails in at least one row."
      stop "alpha fails"
    end if
  end subroutine update_alpha_potential

  subroutine get_all_targets(r_e_new, root_mphi_re, &
                            out_target_rho, out_target_gama, out_target_ww, out_target_sphi)
    real(wp), intent(in) :: r_e_new
    real(wp), intent(in) :: root_mphi_re
    real(wp), intent(out) :: out_target_rho(SDIV,MDIV), out_target_gama(SDIV,MDIV), out_target_ww(SDIV,MDIV), out_target_sphi(SDIV,MDIV)
    real(wp) :: t0, t1, dt_precompute, dt_build, dt_angular, dt_radial, dt_sum
    integer, parameter :: timing_calls = 5
    integer, save :: target_call_count = 0
    real(wp), save :: sum_dt_precompute = 0.e0_wp, sum_dt_build = 0.e0_wp, sum_dt_angular = 0.e0_wp
    real(wp), save :: sum_dt_radial = 0.e0_wp, sum_dt_sum = 0.e0_wp

    if (timing) then
      if (target_call_count == 0) then
        sum_dt_precompute = 0.e0_wp
        sum_dt_build = 0.e0_wp
        sum_dt_angular = 0.e0_wp
        sum_dt_radial = 0.e0_wp
        sum_dt_sum = 0.e0_wp
      end if
      target_call_count = target_call_count + 1
      dt_precompute = 0.e0_wp
      dt_build = 0.e0_wp
      dt_angular = 0.e0_wp
      dt_radial = 0.e0_wp
      dt_sum = 0.e0_wp
      call cpu_time(t0)
    end if
    call precompute_derivatives_and_bessels(r_e_new, root_mphi_re, mr_cache, besseli_cache, besselk_cache, &
        dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
        d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup4_cache)
    if (timing) then
      call cpu_time(t1); dt_precompute = t1 - t0; call cpu_time(t0)
    end if

    call build_source_terms(r_e_new, S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
        dr_s_cache, dr_m_cache, dg_s_cache, dg_m_cache, dww_s_cache, dww_m_cache, ds_s_cache, ds_m_cache, &
        d2g_ss_cache, d2g_mm_cache, e_gsm_cache, e_rsm_cache, e2alpha_r2_cache, Acoup4_cache)
    if (timing) then
      call cpu_time(t1); dt_build = t1 - t0; call cpu_time(t0)
    end if

    call angular_integration(S_metric_rho, S_metric_gama, S_metric_omega, S_metric_sphi, &
        D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi)
    if (timing) then
      call cpu_time(t1); dt_angular = t1 - t0; call cpu_time(t0)
    end if

    call radial_integration(D1_metric_rho, D1_metric_gama, D1_metric_omega, D1_metric_sphi, &
        D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi, root_mphi_re, wfac_cache, &
        besseli_cache, besselk_cache)
    if (timing) then
      call cpu_time(t1); dt_radial = t1 - t0; call cpu_time(t0)
    end if

    call sum_coefficients_and_get_targets(out_target_rho, out_target_gama, out_target_ww, out_target_sphi, &
        D2_metric_rho, D2_metric_gama, D2_metric_omega, D2_metric_sphi)
    if (timing) then
      call cpu_time(t1); dt_sum = t1 - t0; call cpu_time(t0)
    end if

    if (timing) then
      sum_dt_precompute = sum_dt_precompute + dt_precompute
      sum_dt_build = sum_dt_build + dt_build
      sum_dt_angular = sum_dt_angular + dt_angular
      sum_dt_radial = sum_dt_radial + dt_radial
      sum_dt_sum = sum_dt_sum + dt_sum

      write(*,'(A,I0,A,7(1X,ES12.5))') 'get_all_targets call ', target_call_count, ':', &
        dt_precompute, dt_build, dt_angular, dt_radial, dt_sum, &
        dt_precompute + dt_build + dt_angular + dt_radial + dt_sum
      if (target_call_count >= timing_calls) then
        write(*,'(A,I0,A,7(1X,ES12.5))') 'get_all_targets avg over ', timing_calls, ':', &
          sum_dt_precompute / timing_calls, sum_dt_build / timing_calls, sum_dt_angular / timing_calls, &
          sum_dt_radial / timing_calls, sum_dt_sum / timing_calls, &
          (sum_dt_precompute + sum_dt_build + sum_dt_angular + sum_dt_radial + sum_dt_sum) / timing_calls
        stop "Finish profiling."
      end if
    end if
  end subroutine get_all_targets

  subroutine output_helper(D2_metric_rho, D2_metric_omega)
    real(wp), intent(in) :: D2_metric_rho(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1)
    real(wp) :: r_inf, rho_0, t0, t1
    character(32) :: fil1, fil2, fil3, fil4, fil5, fil6, fil7, fil8
    character(64) :: tail_fmt
    integer :: s, unit, ios, sig_digits, field_width

    r_inf = r_e * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1.e0_wp - s_gp(SDIV - 1) ))**s_pwr
    M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2.e0_wp * r_inf**3 * ( C**2 / G / Mass )**3
    S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2.e0_wp * r_inf**5 * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2.e0_wp * r_inf**5 * ( C**2 / G / Mass )**5

    if (.not. output) &
      return

    call cpu_time(t0); write(*,*) " "
    write(*,"(A)",advance='no') " Off-loading data ..."

    open(newunit=unit, file="Cont/moment_tail.dat", status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim("Cont/moment_tail.dat")
      return
    end if
    sig_digits = precision(1.0_wp) - 1
    field_width = sig_digits + 10
    write(tail_fmt,'("(30es",i0,".",i0,"e3)")') field_width, sig_digits
    do s = 1, SDIV-1
      write(unit,tail_fmt) s_gp(s), D2_metric_rho(s,:), D2_metric_omega(s,:)
    end do
    close(unit)

    write(fil1,"(f6.2)") ang_mom
    write(fil2,"(f16.3)") mass_0/MSUN
    write(fil3,"(es15.2)") B_coup
    write(fil4,"(es15.2)") sqrt(mphi_r*1.e10_wp/KAPPA)*l_uni
    rho_0 = n0_at_e( energy(1,1) ) * MB
    write(fil5,"(es15.3)") rho_0
    write(fil6,"(f15.3)") sphi_m
    select case(trim(solver_type))
    case("uniform")
      write(fil7,"(A)") ".dat"
    case("const_j")
      write(fil7,"(f5.2)") A_diff
      write(fil7,"(A)") "_A"//trim(adjustl(fil7))//".dat"
    case("uryu")
      write(fil7,"(f5.2)") lambda1
      write(fil8,"(f5.2)") lambda2
      write(fil7,"(A)") "_L1"//trim(adjustl(fil7))//"_L2"//trim(adjustl(fil8))//".dat"
    case default
      stop "Unknown solver type"
    end select

    call write_output_file("./Cont/"//trim(adjustl(eos_file))//&
      "_J"//trim(adjustl(fil1))//&
      "_Mb"//trim(adjustl(fil2))//"_B"//trim(adjustl(fil3))//&
      "_mphi"//trim(adjustl(fil4))//"_rhoc"//trim(adjustl(fil5))//"_sphim"//trim(adjustl(fil6))//trim(adjustl(fil7)) )

    call write_output_file("./Res/res.dat")

    call cpu_time(t1); write(*,"(A, f12.6, A)") "   took ", t1-t0, " [s]"
  end subroutine output_helper

  subroutine write_output_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios, s, m
    integer :: sig_digits, exp_digits, field_width
    real(wp) :: rho_0_val
    character(len=256) :: header_fmt, data_fmt

    sig_digits = precision(1.0_wp) - 1
    exp_digits = 3
    if (range(1.0_wp) > 999) exp_digits = 4
    field_width = sig_digits + exp_digits + 10
    write(header_fmt,'("(3(i0,3x), 5es",i0,".",i0,"e",i0,")")') field_width, sig_digits, exp_digits
    write(data_fmt,  '("(13es",i0,".",i0,"e",i0,")")') field_width, sig_digits, exp_digits

    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim(filename)
      return
    end if

    write(unit, fmt=header_fmt) SDIV, MDIV, s_pwr, r_e*sqrt(KAPPA)/1.e5_wp, &
            energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(kappa)), Omega_c* (C/sqrt(kappa))

    do s = 1, SDIV
      do m = 1, MDIV
        if (enthalpy(s,m) > enthalpy_min) then
          rho_0_val = n0_at_e(energy(s,m)) * MB
        else
          rho_0_val = 0.e0_wp
        end if

        write(unit, fmt=data_fmt) s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), &
          ww(s,m) * (C/sqrt(kappa)), pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), &
          enthalpy(s,m), rho_0_val, velocity_sq(s,m), omg(s,m) * (C/sqrt(kappa)), &
          sphi(s,m) * sqrt(B_coup)
      end do
    end do

    close(unit)
  end subroutine write_output_file

end module spin_integration
