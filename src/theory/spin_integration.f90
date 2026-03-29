module spin_integration
  use para_mod, only: wp, SDIV, MDIV, LMAX, s_gp, mu, sin_theta, s_pwr, &
                      rho, gama, alpha, ww, omg, sphi, &
                      energy, pressure, enthalpy, velocity_sq, &
                      P_2n, P1_2n_1, l_uni, KAPPA, C, G, MSUN, MB, pi, KSCALE, &
                      r_ratio, r_e, enthalpy_min, &
                      Omega_c, Omega_e, M2, S3, M4, sphi_m, &
                      B_coup, mphi_r, mass, mass_0, ang_mom, &
                      A_diff, lambda1, lambda2, solver_type, output, timing, eos_file
  use eos_mod, only: n0_at_e
  use toolkit_mod, only: bessel_even_tables
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
  elemental real(wp) function pade_exp(x) result(y)
    real(wp), intent(in) :: x
    real(wp) :: t, t2
    t = x * 0.00390625_wp
    t2 = t * t
    y = (12.0_wp + 6.0_wp * t + t2) / (12.0_wp - 6.0_wp * t + t2)
    y = y*y; y = y*y; y = y*y; y = y*y
    y = y*y; y = y*y; y = y*y; y = y*y
  end function

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
    real(wp) :: t0, t1, dt_bessel, dt_first_derivs, dt_second_derivs, dt_cache_fill
    integer, parameter :: timing_calls = 5
    integer, save :: precompute_call_count = 0
    real(wp), save :: sum_dt_bessel = 0.0_wp, sum_dt_first_derivs = 0.0_wp
    real(wp), save :: sum_dt_second_derivs = 0.0_wp, sum_dt_cache_fill = 0.0_wp

    if (timing) then
      if (precompute_call_count == 0) then
        sum_dt_bessel = 0.0_wp
        sum_dt_first_derivs = 0.0_wp
        sum_dt_second_derivs = 0.0_wp
        sum_dt_cache_fill = 0.0_wp
      end if
      precompute_call_count = precompute_call_count + 1
      dt_bessel = 0.0_wp
      dt_first_derivs = 0.0_wp
      dt_second_derivs = 0.0_wp
      dt_cache_fill = 0.0_wp
      call cpu_time(t0)
    end if

    use_massive_scalar = mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp
    if (use_massive_scalar) then
      mr_cache(:) = root_mphi_re * s_gp(:) / (1.0_wp - s_gp(:))
      do s = 1, SDIV
        call bessel_even_tables(mr_cache(s), LMAX, besseli_cache(:,s), besselk_cache(:,s))
      end do
    end if

    if (timing) then
      call cpu_time(t1); dt_bessel = t1 - t0; call cpu_time(t0)
    end if

    call deriv_s_sub(gama, dg_s_cache)
    call deriv_s_sub(rho, dr_s_cache)
    call deriv_s_sub(ww, dww_s_cache)
    call deriv_s_sub(sphi, ds_s_cache)
    call deriv_m_sub(gama, dg_m_cache)
    call deriv_m_sub(rho, dr_m_cache)
    call deriv_m_sub(ww, dww_m_cache)
    call deriv_m_sub(sphi, ds_m_cache)
    if (timing) then
      call cpu_time(t1); dt_first_derivs = t1 - t0; call cpu_time(t0)
    end if

    call deriv_s_sub(dg_s_cache, d2g_ss_cache)
    do m = 1, MDIV
      d2g_ss_cache(:,m) = s1_geom * d2g_ss_cache(:,m) + (1.0_wp - 2.0_wp * s_gp) * dg_s_cache(:,m)
    end do
    call deriv_m_sub(dg_m_cache, d2g_mm_cache)
    do m = 1, MDIV
      d2g_mm_cache(:,m) = m1_geom(m) * d2g_mm_cache(:,m) - 2.0_wp * mu(m) * dg_m_cache(:,m)
    end do

    if (timing) then
      call cpu_time(t1); dt_second_derivs = t1 - t0; call cpu_time(t0)
    end if

    block
      integer :: si, mi
      real(wp) :: r2
      r2 = r_e_new * r_e_new
      do mi = 1, MDIV
        do si = 1, SDIV
          e_gsm_cache(si,mi)      = pade_exp(0.5_wp * gama(si,mi))
          e_rsm_cache(si,mi)      = pade_exp(-rho(si,mi))
          e2alpha_r2_cache(si,mi) = pade_exp(2.0_wp * alpha(si,mi)) * r2
          Acoup4_cache(si,mi)     = pade_exp(-sphi(si,mi)**2 * B_coup)
        end do
      end do
    end block

    if (timing) then
      call cpu_time(t1); dt_cache_fill = t1 - t0
      sum_dt_bessel = sum_dt_bessel + dt_bessel
      sum_dt_first_derivs = sum_dt_first_derivs + dt_first_derivs
      sum_dt_second_derivs = sum_dt_second_derivs + dt_second_derivs
      sum_dt_cache_fill = sum_dt_cache_fill + dt_cache_fill

      write(*,'(A,I0,A,5(1X,ES12.5))') 'precompute call ', precompute_call_count, ':', &
        dt_bessel, dt_first_derivs, dt_second_derivs, dt_cache_fill, &
        dt_bessel + dt_first_derivs + dt_second_derivs + dt_cache_fill
      if (precompute_call_count >= timing_calls) then
        write(*,'(A,I0,A)') 'precompute avg over ', timing_calls, ':'
        write(*,'(A,1X,ES12.5)') '  bessel', sum_dt_bessel / timing_calls
        write(*,'(A,1X,ES12.5)') '  first_derivs', sum_dt_first_derivs / timing_calls
        write(*,'(A,1X,ES12.5)') '  second_derivs', sum_dt_second_derivs / timing_calls
        write(*,'(A,1X,ES12.5)') '  cache_fill', sum_dt_cache_fill / timing_calls
        write(*,'(A,1X,ES12.5)') '  total', &
          (sum_dt_bessel + sum_dt_first_derivs + sum_dt_second_derivs + sum_dt_cache_fill) / timing_calls
      end if
    end if
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
      one_plus_vsq_col = 1.0_wp + vsq_col
      vel_fac_col = 1.0_wp / (1.0_wp - vsq_col)
      matter_sum_col = esm_col + psm_col
      matter_trace_col = esm_col - 3.0_wp * psm_col
      e_rsm2_col = ersm**2
      e2alpha_s2_col = e2ar2 * s2_geom
      source_common_col = 16.0_wp * pi * e2alpha_s2_col * psm_col - 4.0_wp * vphi_col * s2_geom
      dg_s_scaled_col = s1_geom * gs
      dg_m_scaled_col = m1 * gm

      rho_bracket_col = source_common_col &
        - dg_s_scaled_col * (0.5e0_wp * dg_s_scaled_col + 1.0_wp) &
        - gm * (0.5e0_wp * dg_m_scaled_col - mum)

      S_metric_rho(:,m) = egsm * ( &
          8.0_wp * pi * e2alpha_s2_col * matter_sum_col * one_plus_vsq_col * vel_fac_col &
        + s2_geom * m1 * e_rsm2_col * (s1_sq_geom * wws**2 + m1 * wwm**2) &
        + dg_s_scaled_col - mum * gm &
        + rho_col * 0.5e0_wp * rho_bracket_col )

      S_metric_gama(:,m) = egsm * (source_common_col &
        + gama_col * 0.5e0_wp * (source_common_col - 0.5e0_wp * dg_s_scaled_col**2 &
        - 0.5e0_wp * dg_m_scaled_col * gm) )

      if (abs(r_ratio - 1.0_wp) < epsilon(r_ratio)) then
        S_metric_omega(:,m) = 0.0_wp
      else
        omega_matter_col = (one_plus_vsq_col * esm_col + 2.0_wp * vsq_col * psm_col) * vel_fac_col
        omega_bracket_col = -8.0_wp * pi * e2alpha_s2_col * omega_matter_col &
          - s1_geom * (2.0_wp * rs + 0.5e0_wp * gs) &
          + mum * (2.0_wp * rm + 0.5e0_wp * gm) &
          + 0.25e0_wp * s1_sq_geom * (4.0_wp * rs**2 - gs**2) &
          + 0.25e0_wp * m1 * (4.0_wp * rm**2 - gm**2) &
          - m1 * e_rsm2_col * (sgp4_geom * wws**2 + s2_geom * m1 * wwm**2) &
          - 2.0_wp * vphi_col * s2_geom

        S_metric_omega(:,m) = egsm * ersm * ( &
          -16.0_wp * pi * e2alpha_s2_col * (omg_col - ww_col) * matter_sum_col * vel_fac_col &
          + ww_col * omega_bracket_col )
      endif

      if ( mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp ) then
        sphi_source_col = -2.0_wp * pi * B_coup * matter_trace_col + mphi_r
        S_metric_sphi(:,m) = -r_e_new**2 * s2_geom * scal_p_col * mphi_r &
          + e2alpha_s2_col * scal_p_col * sphi_source_col &
          + scal_p_col * (s1_one_minus_s_geom * gs + s1_sq_geom * (0.5e0_wp * gss + 0.25e0_wp * gs**2) &
          + m1 * (0.5e0_wp * gmm + 0.25e0_wp * gm**2) - mum * gm)
      else
        sphi_source_col = -2.0_wp * pi * B_coup * matter_trace_col
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
    call dgemm('N', 'N', SDIV, LMAX+1, MDIV, 1.0_wp, S_metric_rho, SDIV, weighted_even_basis, MDIV, 0.0_wp, D1_metric_rho, SDIV)
    call dgemm('N', 'N', SDIV, LMAX+1, MDIV, 1.0_wp, S_metric_sphi, SDIV, weighted_even_basis, MDIV, 0.0_wp, D1_metric_sphi, SDIV)
    D1_metric_gama = 0.0_wp
    D1_metric_omega = 0.0_wp
    if (LMAX > 0) then
      call dgemm('N', 'N', SDIV, LMAX, MDIV, 1.0_wp, S_metric_gama, SDIV, weighted_gama_basis, MDIV, 0.0_wp, D1_metric_gama(:,2:LMAX+1), SDIV)
      call dgemm('N', 'N', SDIV, LMAX, MDIV, 1.0_wp, S_metric_omega, SDIV, weighted_omega_basis, MDIV, 0.0_wp, D1_metric_omega(:,2:LMAX+1), SDIV)
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
    real(wp) :: f2n_vals(SDIV)
    real(wp) :: left_prefix(SDIV), right_suffix(SDIV), massive_sphi_scale(SDIV)
    logical :: use_massive_scalar

    D2_metric_gama (:,1) = 0.0_wp
    D2_metric_omega(:,1) = 0.0_wp
    use_massive_scalar = mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp
    if (use_massive_scalar) massive_sphi_scale = radial_quad_weights * wfac_cache * root_mphi_re

    f2n_vals = 1.0_wp
    if (use_massive_scalar) then
      do n = 0, LMAX
        call integrate_rho_like(n, f2n_vals, D1_metric_rho(:,n+1), radial_quad_weights, D2_metric_rho(:,n+1))
        call integrate_massive_sphi(n+1, D1_metric_sphi(:,n+1), massive_sphi_scale, D2_metric_sphi(:,n+1))
        f2n_vals(2:) = f2n_vals(2:) * rad_ratio_g(2:)
      end do
    else
      do n = 0, LMAX
        call integrate_rho_like(n, f2n_vals, D1_metric_rho(:,n+1), radial_quad_weights, D2_metric_rho(:,n+1))
        call integrate_rho_like(n, f2n_vals, D1_metric_sphi(:,n+1), radial_quad_weights, D2_metric_sphi(:,n+1))
        f2n_vals(2:) = f2n_vals(2:) * rad_ratio_g(2:)
      end do
    end if

    f2n_vals(1) = 1.0_wp
    f2n_vals(2:) = rad_ratio_g(2:)
    do n = 1, LMAX
      call integrate_multipole_like(n, f2n_vals, D1_metric_gama(:,n+1), radial_quad_weights, rad_inv_s1, .false., D2_metric_gama(:,n+1))
      call integrate_multipole_like(n, f2n_vals, D1_metric_omega(:,n+1), radial_quad_weights, rad_left_rho, .true., D2_metric_omega(:,n+1))
      f2n_vals(2:) = f2n_vals(2:) * rad_ratio_g(2:)
    end do

  contains
    subroutine integrate_rho_like(n_phys, f2n_vals, source_values, source_scale, out_values)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: f2n_vals(SDIV)
      real(wp), intent(in) :: source_values(:), source_scale(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: k
      real(wp) :: weighted_value

      left_prefix(1) = 0.0_wp
      if (n_phys == 0) then
        do k = 2, SDIV
          weighted_value = source_scale(k) * source_values(k)
          left_prefix(k) = left_prefix(k-1) + rad_left_rho(k) * weighted_value
        end do
        right_suffix(SDIV) = rad_inv_s1(SDIV) * source_scale(SDIV) * source_values(SDIV)
        do k = SDIV-1, 2, -1
          weighted_value = source_scale(k) * source_values(k)
          right_suffix(k) = right_suffix(k+1) + rad_inv_s1(k) * weighted_value
        end do
        right_suffix(1) = right_suffix(2)
        out_values(1) = right_suffix(1)
        do k = 2, SDIV
          out_values(k) = rad_ratio_s(k) * left_prefix(k-1) + right_suffix(k)
        end do
      else
        do k = 2, SDIV
          weighted_value = source_scale(k) * source_values(k)
          left_prefix(k) = left_prefix(k-1) + rad_left_rho(k) * weighted_value / f2n_vals(k)
        end do
        right_suffix(SDIV) = rad_inv_s1(SDIV) * f2n_vals(SDIV) * source_scale(SDIV) * source_values(SDIV)
        do k = SDIV-1, 2, -1
          weighted_value = source_scale(k) * source_values(k)
          right_suffix(k) = right_suffix(k+1) + rad_inv_s1(k) * f2n_vals(k) * weighted_value
        end do
        right_suffix(1) = right_suffix(2)
        out_values(1) = 0.0_wp
        do k = 2, SDIV
          out_values(k) = f2n_vals(k) * rad_ratio_s(k) * left_prefix(k-1) + right_suffix(k) / f2n_vals(k)
        end do
      end if
    end subroutine integrate_rho_like

    subroutine integrate_multipole_like(n_phys, f2n_vals, source_values, source_scale, left_prefix_coeff, use_rad_ratio_s_output, out_values)
      integer, intent(in) :: n_phys
      real(wp), intent(in) :: f2n_vals(SDIV)
      real(wp), intent(in) :: source_values(:), source_scale(SDIV), left_prefix_coeff(SDIV)
      logical, intent(in) :: use_rad_ratio_s_output
      real(wp), intent(out) :: out_values(SDIV)
      integer :: k
      real(wp) :: boundary_sum, weighted_value

      left_prefix(1) = 0.0_wp
      if (n_phys == 1) then
        boundary_sum = 0.0_wp
        do k = 2, SDIV
          weighted_value = source_scale(k) * source_values(k)
          left_prefix(k) = left_prefix(k-1) + left_prefix_coeff(k) * weighted_value / f2n_vals(k)
          boundary_sum = boundary_sum + rad_inv_s1(k) * weighted_value
        end do
      else
        do k = 2, SDIV
          weighted_value = source_scale(k) * source_values(k)
          left_prefix(k) = left_prefix(k-1) + left_prefix_coeff(k) * weighted_value / f2n_vals(k)
        end do
      end if
      right_suffix(SDIV) = rad_right_gama(SDIV) * f2n_vals(SDIV) * source_scale(SDIV) * source_values(SDIV)
      do k = SDIV-1, 2, -1
        weighted_value = source_scale(k) * source_values(k)
        right_suffix(k) = right_suffix(k+1) + rad_right_gama(k) * f2n_vals(k) * weighted_value
      end do
      right_suffix(1) = right_suffix(2)

      out_values(1) = merge(boundary_sum, 0.0_wp, n_phys == 1)
      do k = 2, SDIV
        if (use_rad_ratio_s_output) then
          out_values(k) = f2n_vals(k) * rad_ratio_s(k) * left_prefix(k-1) + rad_ratio_g(k) * right_suffix(k) / f2n_vals(k)
        else
          out_values(k) = real(s_pwr, wp) * f2n_vals(k) * left_prefix(k-1) + rad_ratio_g(k) * right_suffix(k) / f2n_vals(k)
        end if
      end do
    end subroutine integrate_multipole_like

    subroutine integrate_massive_sphi(n_idx, source_values, source_scale, out_values)
      integer, intent(in) :: n_idx
      real(wp), intent(in) :: source_values(:), source_scale(SDIV)
      real(wp), intent(out) :: out_values(SDIV)
      integer :: s
      real(wp) :: weighted_value

      left_prefix(1) = source_scale(1) * source_values(1) * besseli_cache(n_idx,1)
      do s = 2, SDIV
        weighted_value = source_scale(s) * source_values(s)
        left_prefix(s) = left_prefix(s-1) + weighted_value * besseli_cache(n_idx,s)
      end do

      right_suffix(SDIV) = source_scale(SDIV) * source_values(SDIV) * besselk_cache(n_idx,SDIV)
      do s = SDIV-1, 1, -1
        weighted_value = source_scale(s) * source_values(s)
        right_suffix(s) = right_suffix(s+1) + weighted_value * besselk_cache(n_idx,s)
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
    integer :: si, mi
    real(wp) :: emhalf
    logical :: massive_sphi, spherical

    massive_sphi = mphi_r > (mphi_tran / l_uni)**2 * KAPPA / 1.e10_wp
    spherical = abs(r_ratio - 1.0_wp) < epsilon(r_ratio)

    call dgemm('N','T', SDIV, MDIV, LMAX+1, 1.0_wp, D2_metric_rho, SDIV, P_2n, MDIV, 0.0_wp, out_target_rho, SDIV)

    if (massive_sphi) then
      call dgemm('N','T', SDIV, MDIV, LMAX+1, -1.0_wp, D2_metric_sphi, SDIV, recon_even_massive_basis, MDIV, 0.0_wp, out_target_sphi, SDIV)
    else
      call dgemm('N','T', SDIV, MDIV, LMAX+1, -1.0_wp, D2_metric_sphi, SDIV, P_2n, MDIV, 0.0_wp, out_target_sphi, SDIV)
    endif

    if (spherical) then
      call dgemm('N','T', SDIV, MDIV, 1, -(2.0_wp/pi), D2_metric_gama(:,2:2), SDIV, recon_gama_basis(:,1:1), MDIV, 0.0_wp, out_target_gama, SDIV)
      if (massive_sphi) then
        do mi = 1, MDIV
          do si = 1, SDIV
            emhalf = 1.0_wp / e_gsm_cache(si,mi)
            out_target_rho(si,mi)  = -emhalf * out_target_rho(si,mi)
            out_target_gama(si,mi) = emhalf * out_target_gama(si,mi)
            out_target_sphi(si,mi) = emhalf * out_target_sphi(si,mi)
          end do
        end do
      else
        do mi = 1, MDIV
          do si = 1, SDIV
            out_target_rho(si,mi) = -out_target_rho(si,mi) / e_gsm_cache(si,mi)
            out_target_gama(si,mi) = out_target_gama(si,mi) / e_gsm_cache(si,mi)
          end do
        end do
      end if
      out_target_ww = 0.0_wp
      return
    endif

    call dgemm('N','T', SDIV, MDIV, LMAX, 1.0_wp, D2_metric_omega(:,2:LMAX+1), SDIV, recon_omega_basis, MDIV, 0.0_wp, out_target_ww, SDIV)
    call dgemm('N','T', SDIV, MDIV, LMAX, -(2.0_wp/pi), D2_metric_gama(:,2:LMAX+1), SDIV, recon_gama_basis, MDIV, 0.0_wp, out_target_gama, SDIV)

    if (massive_sphi) then
      do mi = 1, MDIV
        do si = 1, SDIV
          emhalf = 1.0_wp / e_gsm_cache(si,mi)
          out_target_rho(si,mi)  = -emhalf * out_target_rho(si,mi)
          out_target_gama(si,mi) = emhalf * out_target_gama(si,mi)
          out_target_ww(si,mi)   = emhalf / e_rsm_cache(si,mi) * out_target_ww(si,mi)
          out_target_sphi(si,mi) = emhalf * out_target_sphi(si,mi)
        end do
      end do
    else
      do mi = 1, MDIV
        do si = 1, SDIV
          emhalf = 1.0_wp / e_gsm_cache(si,mi)
          out_target_rho(si,mi)  = -emhalf * out_target_rho(si,mi)
          out_target_gama(si,mi) = emhalf * out_target_gama(si,mi)
          out_target_ww(si,mi)   = emhalf / e_rsm_cache(si,mi) * out_target_ww(si,mi)
        end do
      end do
    end if
  end subroutine sum_coefficients_and_get_targets

  subroutine update_alpha_potential(r_e_new, dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
                                    ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_rsm_cache)
    real(wp), intent(in) :: r_e_new
    real(wp), intent(in) :: dg_s_cache(:,:), dg_m_cache(:,:), dr_s_cache(:,:), dr_m_cache(:,:), dww_s_cache(:,:), dww_m_cache(:,:)
    real(wp), intent(in) :: ds_s_cache(:,:), ds_m_cache(:,:)
    real(wp), intent(in) :: d2g_ss_cache(:,:), d2g_mm_cache(:,:), e_rsm_cache(:,:)
    integer :: m, s
    real(wp) :: m1, mu_m
    real(wp) :: t0, t1, dt_deriv_m, dt_column_loop, dt_integrate, dt_adjust
    integer, save :: alpha_call_count = 0
    real(wp), save :: sum_dt_deriv_m = 0.0_wp, sum_dt_column_loop = 0.0_wp
    real(wp), save :: sum_dt_integrate = 0.0_wp, sum_dt_adjust = 0.0_wp
    real(wp), dimension(SDIV,MDIV) :: da_dm, d_gama_sm_all
    real(wp) :: adj_const(SDIV)
    real(wp) :: gs, gm, rs, rm, ss, sm, gsm, wws, wwm, gss, gmm, e_cache
    real(wp) :: sg, s1, sg_ratio, sg4, numer_m, one_plus_s1dgs, inv_denom
    real(wp) :: temp1, temp3, temp4, temp5, temp6, temp7, temp8, temp9

    if (timing) then
      if (alpha_call_count == 0) then
        sum_dt_deriv_m = 0.0_wp
        sum_dt_column_loop = 0.0_wp
        sum_dt_integrate = 0.0_wp
        sum_dt_adjust = 0.0_wp
      end if
      alpha_call_count = alpha_call_count + 1
      dt_deriv_m = 0.0_wp
      dt_column_loop = 0.0_wp
      dt_integrate = 0.0_wp
      dt_adjust = 0.0_wp
      call cpu_time(t0)
    end if

    alpha(:,:) = 0.0_wp
    if (abs(r_ratio - 1.0_wp) < epsilon(r_ratio)) then
      return
    else
      da_dm(1,:) = 0.0e0_wp
      call deriv_m_sub(dg_s_cache, d_gama_sm_all)
      if (timing) then
        call cpu_time(t1); dt_deriv_m = t1 - t0; call cpu_time(t0)
      end if
      do m = 1, MDIV
        mu_m = mu(m)
        m1   = m1_geom(m)
        da_dm(1,m) = 0.0_wp
        do s = 2, SDIV
          gs = dg_s_cache(s,m)
          gm = dg_m_cache(s,m)
          rs = dr_s_cache(s,m)
          rm = dr_m_cache(s,m)
          ss = ds_s_cache(s,m)
          sm = ds_m_cache(s,m)
          gsm = d_gama_sm_all(s,m)
          wws = dww_s_cache(s,m)
          wwm = dww_m_cache(s,m)
          gss = d2g_ss_cache(s,m)
          gmm = d2g_mm_cache(s,m)
          e_cache = e_rsm_cache(s,m)
          sg = s_gp(s)
          s1 = s1_geom(s)
          sg_ratio = sgp_term_2d_cache(s,1)
          sg4 = sgp4_geom(s)

          numer_m = -mu_m + m1 * gm
          one_plus_s1dgs = 1.0_wp + s1 * gs
          inv_denom = 1.0_wp / (m1 * one_plus_s1dgs**2 + numer_m**2)

          temp1 = 2.0_wp * sg**2 * sg_ratio * m1 * wws * wwm * one_plus_s1dgs &
            - (sg4 * wws**2 - (sg * wwm * sg_ratio)**2 * m1) * numer_m
          temp3 = s1 * gss + (s1 * gs)**2
          temp4 = gm * numer_m
          temp5 = ((s1 * (rs + gs))**2 - m1 * (rm + gm)**2) * numer_m
          temp6 = s1 * m1 * ((rs + gs) * (rm + gm) / 2.0_wp + gsm + gs * gm) * one_plus_s1dgs
          temp7 = s1 * mu_m * gs * one_plus_s1dgs
          temp8 = m1 * e_cache * e_cache
          temp9 = -inv_denom * numer_m * ((s1 * ss)**2 - m1 * sm**2) &
            - 2.0_wp * m1 * s1 * one_plus_s1dgs * sm * ss

          da_dm(s,m) = -(rm + gm) / 2.0_wp &
            - inv_denom * ((temp3 - gmm - temp4) * numer_m / 2.0_wp &
            + temp5 / 4.0_wp - temp6 + temp7 + temp8 * temp1 / 4.0_wp) + temp9
        end do
      end do
      if (timing) then
        call cpu_time(t1); dt_column_loop = t1 - t0; call cpu_time(t0)
      end if

      do m = 1, MDIV-1
        alpha(:,m+1) = alpha(:,m) + (mu(m+1) - mu(m)) * ( da_dm(:,m+1) + da_dm(:,m) ) * 0.5e0_wp
      enddo
      if (timing) then
        call cpu_time(t1); dt_integrate = t1 - t0; call cpu_time(t0)
      end if

      alpha(SDIV,:) = 0.0_wp
      adj_const = alpha(:,MDIV) - ( gama(:,MDIV) - rho(:,MDIV) )/2.0_wp
      do m = 1, MDIV
        alpha(:,m) = alpha(:,m) - adj_const
      end do
      if (timing) then
        call cpu_time(t1); dt_adjust = t1 - t0
      end if
    end if
    if (timing) then
      sum_dt_deriv_m = sum_dt_deriv_m + dt_deriv_m
      sum_dt_column_loop = sum_dt_column_loop + dt_column_loop
      sum_dt_integrate = sum_dt_integrate + dt_integrate
      sum_dt_adjust = sum_dt_adjust + dt_adjust
      if (alpha_call_count<4) return
      write(*,'(A,I0,A,5(1X,ES12.5))') 'update_alpha_potential call ', alpha_call_count, ':', &
        dt_deriv_m, dt_column_loop, dt_integrate, dt_adjust, &
        dt_deriv_m + dt_column_loop + dt_integrate + dt_adjust
      write(*,'(A,I0,A)') 'update_alpha_potential running avg over ', alpha_call_count, ':'
      write(*,'(A,1X,ES12.5)') '  deriv_m_sub', sum_dt_deriv_m / alpha_call_count
      write(*,'(A,1X,ES12.5)') '  column_loop', sum_dt_column_loop / alpha_call_count
      write(*,'(A,1X,ES12.5)') '  integrate_mu', sum_dt_integrate / alpha_call_count
      write(*,'(A,1X,ES12.5)') '  adjust_boundary', sum_dt_adjust / alpha_call_count
      write(*,'(A,1X,ES12.5)') '  total', &
        (sum_dt_deriv_m + sum_dt_column_loop + sum_dt_integrate + sum_dt_adjust) / alpha_call_count
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
    real(wp), save :: sum_dt_precompute = 0.0_wp, sum_dt_build = 0.0_wp, sum_dt_angular = 0.0_wp
    real(wp), save :: sum_dt_radial = 0.0_wp, sum_dt_sum = 0.0_wp

    if (timing) then
      if (target_call_count == 0) then
        sum_dt_precompute = 0.0_wp
        sum_dt_build = 0.0_wp
        sum_dt_angular = 0.0_wp
        sum_dt_radial = 0.0_wp
        sum_dt_sum = 0.0_wp
      end if
      target_call_count = target_call_count + 1
      dt_precompute = 0.0_wp
      dt_build = 0.0_wp
      dt_angular = 0.0_wp
      dt_radial = 0.0_wp
      dt_sum = 0.0_wp
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
        write(*,'(A,I0,A)') 'get_all_targets avg over ', timing_calls, ':'
        write(*,'(A,1X,ES12.5)') '  precompute', sum_dt_precompute / timing_calls
        write(*,'(A,1X,ES12.5)') '  build_source_terms', sum_dt_build / timing_calls
        write(*,'(A,1X,ES12.5)') '  angular_integration', sum_dt_angular / timing_calls
        write(*,'(A,1X,ES12.5)') '  radial_integration', sum_dt_radial / timing_calls
        write(*,'(A,1X,ES12.5)') '  sum_coefficients_and_get_targets', sum_dt_sum / timing_calls
        write(*,'(A,1X,ES12.5)') '  total', &
          (sum_dt_precompute + sum_dt_build + sum_dt_angular + sum_dt_radial + sum_dt_sum) / timing_calls
        stop "get_all_targets: profiling window complete"
      end if
    end if
  end subroutine get_all_targets

  subroutine output_helper(D2_metric_rho, D2_metric_omega)
    real(wp), intent(in) :: D2_metric_rho(SDIV,LMAX+1), D2_metric_omega(SDIV,LMAX+1)
    real(wp) :: r_inf, rho_0, t0, t1
    character(512) :: fname
    character(64) :: tail_fmt
    integer :: s, unit, ios, sig_digits, field_width

    r_inf = r_e * sqrt(KAPPA) * (s_gp(SDIV - 1) / ( 1.0_wp - s_gp(SDIV - 1) ))**s_pwr
    M2 = - D2_metric_rho  (SDIV-1,1+1 ) / 2.0_wp * r_inf**3 * ( C**2 / G / Mass )**3
    S3 = - D2_metric_omega(SDIV-1,2+1 ) / 2.0_wp * r_inf**5 * ( C**2 / G / Mass )**4 / sqrt(KAPPA)
    M4 =   D2_metric_rho  (SDIV-1,2+1 ) / 2.0_wp * r_inf**5 * ( C**2 / G / Mass )**5

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

    rho_0 = n0_at_e(energy(1,1)) * MB
    write(fname,"(A, A, A, f0.2, A, f0.3, A, es0.2e2, A,es0.2e2, A, es0.3e2, A, f0.3)") &
      "./Cont/", trim(eos_file), &
      "_J",    ang_mom, &
      "_Mb",   mass_0/MSUN, &
      "_B",    B_coup, &
      "_mphi", sqrt(mphi_r*1.e10_wp/KAPPA)*l_uni, &
      "_rhoc", rho_0, &
      "_sphim",sphi_m
    fname = adjustl(fname)
    select case(trim(solver_type))
    case("uniform")
    case("const_j");  write(fname,"(A,A,f5.2)") trim(fname), "_A",  A_diff
    case("uryu");     write(fname,"(A,A,f5.2,A,f5.2)") trim(fname), "_L1", lambda1, "_L2", lambda2
    case default;     stop "Unknown solver type"
    end select
    fname = trim(fname)//".dat"

    call write_output_file(trim(fname))

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
          rho_0_val = 0.0_wp
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
