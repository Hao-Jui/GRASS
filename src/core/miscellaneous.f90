module miscellaneous_mod
  use para_mod, only: i_isco_m, i_isco_p, res, &
                      V_rr_m, V_rr_p, v_minus, v_plus, &
                      eos_file, chi, pi, C, KAPPA, MSUN, &
                      Omega_c, Omega_e, Omega_K, &
                      s_gp, s_pwr, SDIV, &
                      r_e, mass, Mb_goal, sphi_c
  use nag_compat_mod, only: e02baf, e02bbf
  implicit none

contains
  subroutine write_eq_profile(file_name, s_max, column1, column2, column3, column4, &
                            column5, column6, column7, column8, column9, column10, &
                            column11, column12, column13, column14, column15, column16)
    use para_mod, only : s_gp, s_pwr, SDIV
    implicit none
    character(*), intent(in) :: file_name
    integer, intent(in) :: s_max
    real(8), intent(in) :: column1(:)
    real(8), intent(in), optional :: column2(:), column3(:), column4(:), column5(:), &
                                    column6(:), column7(:), column8(:), column9(:), column10(:), &
                                    column11(:), column12(:), column13(:), column14(:), column15(:), &
                                    column16(:)
    integer :: s, unit, ios, nvals, n_points
    real(8) :: row_values(99)

    if (s_max < 1) return
    if (s_max > SDIV) write(*,*) "write_eq_profile: s_max truncated for file ", trim(file_name)
    n_points = min(s_max, SDIV)

    if (size(column1) < n_points) then
      write(*,*) "write_eq_profile: column1 too short for file ", trim(file_name)
      return
    end if

    if (.not. check_column(column2,  "column2",  n_points)) return
    if (.not. check_column(column3,  "column3",  n_points)) return
    if (.not. check_column(column4,  "column4",  n_points)) return
    if (.not. check_column(column5,  "column5",  n_points)) return
    if (.not. check_column(column6,  "column6",  n_points)) return
    if (.not. check_column(column7,  "column7",  n_points)) return
    if (.not. check_column(column8,  "column8",  n_points)) return
    if (.not. check_column(column9,  "column9",  n_points)) return
    if (.not. check_column(column10, "column10", n_points)) return
    if (.not. check_column(column11, "column11", n_points)) return
    if (.not. check_column(column12, "column12", n_points)) return
    if (.not. check_column(column13, "column13", n_points)) return
    if (.not. check_column(column14, "column14", n_points)) return
    if (.not. check_column(column15, "column15", n_points)) return
    if (.not. check_column(column16, "column16", n_points)) return

    open(newunit=unit, file=file_name, status="replace", action="write", iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim(file_name)
      return
    end if

    do s = 1, n_points
      row_values = 0.d0
      nvals = 2
      row_values(1) = (s_gp(s)/(1.d0 - s_gp(s)))**s_pwr
      row_values(2) = column1(s)
      if (present(column2))  then; nvals=nvals+1; row_values(nvals)=column2(s);  end if
      if (present(column3))  then; nvals=nvals+1; row_values(nvals)=column3(s);  end if
      if (present(column4))  then; nvals=nvals+1; row_values(nvals)=column4(s);  end if
      if (present(column5))  then; nvals=nvals+1; row_values(nvals)=column5(s);  end if
      if (present(column6))  then; nvals=nvals+1; row_values(nvals)=column6(s);  end if
      if (present(column7))  then; nvals=nvals+1; row_values(nvals)=column7(s);  end if
      if (present(column8))  then; nvals=nvals+1; row_values(nvals)=column8(s);  end if
      if (present(column9))  then; nvals=nvals+1; row_values(nvals)=column9(s);  end if
      if (present(column10)) then; nvals=nvals+1; row_values(nvals)=column10(s); end if
      if (present(column11)) then; nvals=nvals+1; row_values(nvals)=column11(s); end if
      if (present(column12)) then; nvals=nvals+1; row_values(nvals)=column12(s); end if
      if (present(column13)) then; nvals=nvals+1; row_values(nvals)=column13(s); end if
      if (present(column14)) then; nvals=nvals+1; row_values(nvals)=column14(s); end if
      if (present(column15)) then; nvals=nvals+1; row_values(nvals)=column15(s); end if
      if (present(column16)) then; nvals=nvals+1; row_values(nvals)=column16(s); end if
      write(unit,"(99es18.9e3)") row_values(1:nvals)
    end do

    close(unit)
  contains
    logical function check_column(col, name, n_points)
      implicit none
      real(8), intent(in), optional :: col(:)
      character(*), intent(in)      :: name
      integer, intent(in)           :: n_points
      ! Default to OK
      check_column = .true.
      if (present(col)) then
        if (size(col) < n_points) then
            write(*,*) "write_eq_profile: ", trim(name), " too short for file."
            check_column = .false.
        end if
      end if
    end function check_column
  end subroutine write_eq_profile

  subroutine print_converged_block(rho0, ee)
  use para_mod, only: has_scalar, rho_uni, MB, C, KSCALE, &
                      r_ratio, Omega_c, pi, KAPPA, Omega_e, Omega_K, &
                      mass, mass_0, ang_mom, chi, B_coup, &
                      mphi_r, l_uni, sphi_c, sphi_m, &
                      I_inertia, M2, S3, M4, T_kin, mass_p, &
                      r_e, r_circ, Fmax_h, h_center
  implicit none
  real(8), intent(in) :: rho0, ee
  integer :: i, prop_unit, out_unit, ios

  if (has_scalar) then
    open(newunit=prop_unit, file="./Cont/properties_st.dat", status="replace", action="write", iostat=ios)
  else
    open(newunit=prop_unit, file="./Cont/properties.dat", status="replace", action="write", iostat=ios)
  end if
  if (ios /= 0) then
    write(*,*) "print_converged_block: failed to open properties file"
    return
  end if

  write(*,*) " ===================================="
  write(*,*) "              Converged              "
  do i = 1, 2
     if (i == 1) then
       out_unit = 6
     else
       out_unit = prop_unit
     end if
     write(out_unit,"(A18,ES18.9,A8,ES18.9)")  &
          "   Central rho =", rho0*MB,"g/cm^3", rho0*MB*rho_uni
     write(out_unit,"(A18,ES18.9,A10)")        &
          "Central energy =", ee/(C*C*KSCALE), "g/cm^3"
     write(out_unit,"(A18,ES18.9)")            "   Axial ratio =", r_ratio
     write(out_unit,"(A18,F18.9,A4)")          " Central Omega =", Omega_c/(2.d0*pi)*(C/sqrt(kappa)), "Hz"
     write(out_unit,"(A18,F18.9,A4)")          " Equator Omega =", Omega_e/(2.d0*pi)*(C/sqrt(kappa)), "Hz"
     if (.not. has_scalar) then
       write(out_unit,"(A18,F18.9,A4)")        "   Kepler Omega =", Omega_K/(2.d0*pi), "Hz"
     end if
     write(out_unit,"(A18,F18.9,A4)")          "      ADM Mass =", mass/MSUN, "M_o"
     write(out_unit,"(A18,F18.9,A16,ES18.9,A2)") &
          " Baryonic Mass =", mass_0/MSUN, "M_o ( binding =", (mass - mass_0)/MSUN, " )"
     write(out_unit,"(A18,F18.9,A10,F7.4,A2)") &
          " Ang. Momentum =", ang_mom, " ( chi =", chi, ")"
     if (has_scalar) then
        write(out_unit,"(A18,ES18.9)")         "    Coupling B =", B_coup
        write(out_unit,"(A18,ES18.9)")         "   Scalar mass =", sqrt(mphi_r*1.d10/KAPPA)*l_uni
        write(out_unit,"(A18,ES18.9)")         "     varphi(0) =", sphi_c
        write(out_unit,"(A18,ES18.9)")         "    varphi_max =", sphi_m
     endif
     write(out_unit,"(A18,F18.9)")             "   Slow rot. I =", I_inertia
     write(out_unit,"(A18,F18.9)")             "        M2/M^3 =", M2
     write(out_unit,"(A18,F18.9)")             "        S3/M^4 =", S3
     write(out_unit,"(A18,F18.9)")             "        M4/M^5 =", M4
     write(out_unit,"(A18,F18.9)")             "           T/W =", T_kin/abs(mass_p - mass + T_kin)
     write(out_unit,"(A18,F18.9,A4)")          "       Coord R =", r_e*sqrt(KAPPA)/1.d5,"km"
     write(out_unit,"(A18,F18.9,A4)")          "       Areal R =", r_circ/1.d5,"km"
  end do

  write(*,*) " "
  write(*,*) "In code unit:"
  write(*,"(A6,ES18.9,2X,A4,ES18.9,2X,A5,ES18.9,2X,A8,ES18.9)")  &
       "h_c", h_center, "r_e", r_e, "Fmax", Fmax_h, "Omega_e", Omega_e*r_e
  write(*,*) " ===================================="
  close(prop_unit)
end subroutine print_converged_block
  
  subroutine initial_data_for_spec(file_name, var1, var2, var3, var4, var5, var6)
    use para_mod, only : s_gp, mu, SDIV,MDIV, l_uni, ang_mom, mass_0, mass
    implicit none
    character(*), intent(in) :: file_name
    real(8), intent(in) :: var1(:,:), var2(:,:), var3(:,:), var4(:,:), var5(:,:), var6(:,:)
    integer :: s, m, unit, ios
    real(8), parameter :: K_km = 218.04217865726338d0
    real(8) :: r_s 
    open(newunit=unit, file=file_name, status="replace", action="write", iostat=ios)
    if (ios /= 0) then
      write(*,*) "write_eq_profile: failed to open file ", trim(file_name)
      return
    end if

    write(unit,"(3(A,e10.4))") "BaryM = ", mass_0 / MSUN * l_uni / sqrt(K_km), &
                            "  AngM = ", ang_mom * l_uni**2 / K_km,         &
                            "   E = ", mass / MSUN * l_uni / sqrt(K_km)
    write(unit,"(2i12)") SDIV, MDIV
    write(unit,"(2es21.12)") r_e * sqrt(KAPPA) / 1.d5 / sqrt(K_km), 1.0d0, &
                             omega_c * ( sqrt(K_km) / sqrt(kappa) ), 0.0d0
    do s = 1, SDIV
      r_s = r_e * sqrt(KAPPA) / 1.d5 / sqrt(K_km) * s_gp(s) / (1.d0 - s_gp(s))
      do m = 1, MDIV
        write(unit,"(4es22.12)") r_s, mu(m), var1(s,m), var2(s,m), &
                                 var3(s,m), var4(s,m), var5(s,m), var6(s,m)
      enddo
    enddo
    close(unit)
  end subroutine initial_data_for_spec

  subroutine log_kepler_sequence()
    integer :: i
    real(8) :: min_Vrr

    min_Vrr = 1.d10
    i_isco_m = 1
    do i = res, int(5.d0 * dble(res) / 3.d0)
      if (min_Vrr > abs(V_rr_m(i))) then
        i_isco_m = i
        min_Vrr = abs(V_rr_m(i))
      endif
    enddo

    min_Vrr = 1.d10
    i_isco_p = 1
    do i = res, int(5.d0 * dble(res) / 3.d0)
      if (min_Vrr > abs(V_rr_p(i))) then
        i_isco_p = i
        min_Vrr = abs(V_rr_p(i))
      endif
    enddo

    open(771,file="./Cont/Kep_"//trim(adjustl(eos_file))//".log",position='append')
    write(771,"(99es18.9)") omega_c / 2.d0 / pi * (C/sqrt(kappa)), &
                            chi, &
                            s_gp(i_isco_m) / (1.d0-s_gp(i_isco_m)), &
                            s_gp(i_isco_p) / (1.d0-s_gp(i_isco_p)), &
                            r_e*sqrt(KAPPA) / 1.d5, &
                            mass / MSUN, &
                            (C/sqrt(kappa)) * v_minus(i_isco_m) / r_e, &
                            (C/sqrt(kappa)) * v_plus(i_isco_p) / r_e, &
                            ( Omega_e * (C/sqrt(kappa)) ) / Omega_K, Mb_goal, &
                            sphi_c
    close(771)
    Mb_goal = Mb_goal + 0.05d0
  end subroutine log_kepler_sequence

  subroutine spectral_tail_fit(field, ell, r_e_current, coeff_leading, coeff_next, success)
    real(8), intent(in) :: field(SDIV)
    integer, intent(in) :: ell
    real(8), intent(in) :: r_e_current
    real(8), intent(out) :: coeff_leading, coeff_next
    logical, intent(out) :: success

    integer, parameter :: tail_points_min = 6
    integer, parameter :: tail_points_max = 24
    integer :: start_idx, end_idx, idx, count, slope_count
    real(8) :: ratio, radius_factor, radius_phys, sqrt_kappa
    real(8) :: weight0, weight1, value
    real(8) :: inv_exponent0, inv_exponent1
    real(8) :: log_slope_sum, avg_slope, denom_log
    real(8) :: radii_samples(tail_points_max), field_samples(tail_points_max)
    real(8) :: weight_samples(tail_points_max), segment_slope(tail_points_max)
    logical :: segment_valid(tail_points_max)
    real(8) :: slope_sum, weight_sum
    real(8) :: fit_slope, slope_residual, fit_intercept, weight
    real(8) :: sum_w, sum_wr, sum_wr2, sum_wf, sum_wrf, denom_fit
    real(8) :: residual, residual_norm, field_norm, sanity_ratio
    real(8) :: sum_num, sum_den, log_r, log_f
    integer :: effective_points
    real(8) :: sum_sign

    coeff_leading = 0.d0
    coeff_next    = 0.d0
    success       = .false.

    if (SDIV <= 2) return

    sqrt_kappa = sqrt(KAPPA)
    end_idx = max(1, SDIV - 5)
    start_idx = max(1, end_idx - tail_points_max + 1)

    count = 0

    inv_exponent0 = -(ell + 1)
    inv_exponent1 = -(ell + 3)

    do idx = start_idx, end_idx
      ratio = s_gp(idx)
      if (ratio <= 0.d0 .or. ratio >= 1.d0) cycle
      radius_factor = ratio / (1.d0 - ratio)
      radius_phys = r_e_current * sqrt_kappa * radius_factor**s_pwr
      if (radius_phys <= 0.d0) cycle

      value = field(idx)
      if (count < tail_points_max) then
        count = count + 1
        radii_samples(count) = radius_phys
        field_samples(count) = value
      else
        cycle
      end if
    end do

    if (count < tail_points_min) return

    sum_sign = 0.d0
    do idx = 1, count
      sum_sign = sum_sign + field_samples(idx)
    end do
    if (sum_sign >= 0.d0) return

    segment_valid(:) = .false.
    segment_slope(:) = 0.d0
    log_slope_sum = 0.d0
    slope_count   = 0
    do idx = 2, count
      if (abs(field_samples(idx-1)) < epsilon(field_samples(idx-1)) .or. &
          abs(field_samples(idx)) < epsilon(field_samples(idx))) cycle
      if (field_samples(idx-1) * field_samples(idx) > 0.d0) then
        denom_log = log(radii_samples(idx)) - log(radii_samples(idx-1))
        if (abs(denom_log) <= tiny(1.d0)) cycle
        segment_slope(idx) = (log(abs(field_samples(idx))) - log(abs(field_samples(idx-1)))) / denom_log
        segment_valid(idx) = .true.
        log_slope_sum = log_slope_sum + segment_slope(idx)
        slope_count   = slope_count + 1
      end if
    end do
    if (slope_count > 0) then
      avg_slope = log_slope_sum / dble(slope_count)
      if (abs(avg_slope + dble(ell + 1)) > 0.4d0) then
        coeff_leading = 0.d0
        coeff_next    = 0.d0
        success = .false.
        return
      end if
    end if

    weight_samples(:) = 1.d0
    effective_points = 0
    do idx = 1, count
      slope_sum = 0.d0
      weight_sum = 0.d0
      if (idx > 1 .and. segment_valid(idx)) then
        slope_sum = slope_sum + segment_slope(idx)
        weight_sum = weight_sum + 1.d0
      end if
      if (idx < count .and. segment_valid(idx + 1)) then
        slope_sum = slope_sum + segment_slope(idx + 1)
        weight_sum = weight_sum + 1.d0
      end if
      if (weight_sum > 0.d0) then
        fit_slope = slope_sum / weight_sum
        slope_residual = fit_slope + dble(ell + 1)
        weight_samples(idx) = 1.d0 / (1.d0 + (slope_residual / 0.2d0)**2)
      else
        weight_samples(idx) = 1.d0
      end if
      if (abs(field_samples(idx)) <= tiny(1.d0)) weight_samples(idx) = 0.d0
      if (weight_samples(idx) > 1.d-6) effective_points = effective_points + 1
    end do

    if (effective_points < tail_points_min) return

    sum_w   = 0.d0
    sum_wr  = 0.d0
    sum_wr2 = 0.d0
    sum_wf  = 0.d0
    sum_wrf = 0.d0
    do idx = 1, count
      if (weight_samples(idx) <= tiny(1.d0)) cycle
      log_r = log(radii_samples(idx))
      log_f = log(abs(field_samples(idx)))
      weight = weight_samples(idx)

      sum_w   = sum_w   + weight
      sum_wr  = sum_wr  + weight * log_r
      sum_wr2 = sum_wr2 + weight * log_r * log_r
      sum_wf  = sum_wf  + weight * log_f
      sum_wrf = sum_wrf + weight * log_r * log_f
    end do

    if (sum_w <= tiny(1.d0)) return

    denom_fit = sum_w * sum_wr2 - sum_wr * sum_wr
    if (abs(denom_fit) <= tiny(1.d0)) return

    fit_slope = (sum_w * sum_wrf - sum_wr * sum_wf) / denom_fit
    slope_residual = fit_slope + dble(ell + 1)
    if (abs(slope_residual) > 0.25d0) return

    fit_intercept = (sum_wf - fit_slope * sum_wr) / sum_w

    coeff_leading = -exp(fit_intercept)
    if (abs(coeff_leading) <= tiny(1.d0)) return

    sum_num = 0.d0
    sum_den = 0.d0
    do idx = 1, count
      if (weight_samples(idx) <= tiny(1.d0)) cycle
      radius_phys = radii_samples(idx)
      value = field_samples(idx)
      weight0 = radius_phys**inv_exponent0
      weight1 = radius_phys**inv_exponent1
      weight = weight_samples(idx)
      residual = value - coeff_leading * weight0

      sum_num = sum_num + weight * weight1 * residual
      sum_den = sum_den + weight * weight1 * weight1
    end do

    if (sum_den > tiny(1.d0)) then
      coeff_next = sum_num / sum_den
    else
      coeff_next = 0.d0
    end if

    residual_norm = 0.d0
    field_norm    = 0.d0
    do idx = 1, count
      if (weight_samples(idx) <= tiny(1.d0)) cycle
      radius_phys = radii_samples(idx)
      value = field_samples(idx)
      weight0 = radius_phys**inv_exponent0
      weight1 = radius_phys**inv_exponent1
      weight = weight_samples(idx)
      residual = value - (coeff_leading * weight0 + coeff_next * weight1)
      residual_norm = residual_norm + weight * residual * residual
      field_norm    = field_norm    + weight * value * value
    end do

    if (field_norm > tiny(1.d0)) then
      sanity_ratio = residual_norm / field_norm
      if (sanity_ratio > 0.1d0) then
        coeff_leading = 0.d0
        coeff_next    = 0.d0
        success = .false.
        return
      end if
    end if

    success = coeff_leading < 0.d0
  end subroutine spectral_tail_fit

  subroutine composite_richardson(field, ell, r_e_current, coeff_leading, coeff_next, success)
    real(8), intent(in) :: field(SDIV)
    integer, intent(in) :: ell
    real(8), intent(in) :: r_e_current
    real(8), intent(out) :: coeff_leading, coeff_next
    logical, intent(out) :: success

    integer, parameter :: sample_points = 24
    integer :: start_idx, end_idx, idx, count, res_count
    real(8) :: ratio, radius_factor, radius_phys, sqrt_kappa
    real(8) :: sum_x, sum_xx, sum_y, sum_xy, denom
    real(8) :: res_sum_x, res_sum_xx, res_sum_y, res_sum_xy, res_denom
    real(8) :: residual, g_value, x_value
    real(8) :: radius_vals(sample_points), field_vals(sample_points), x_vals(sample_points)

    coeff_leading = 0.d0
    coeff_next    = 0.d0
    success       = .false.

    if (SDIV <= 2) return

    sqrt_kappa = sqrt(KAPPA)
    end_idx = max(1, SDIV - sample_points)
    start_idx = max(1, end_idx - sample_points + 1)

    sum_x  = 0.d0
    sum_xx = 0.d0
    sum_y  = 0.d0
    sum_xy = 0.d0
    count  = 0

    do idx = start_idx, end_idx
      ratio = s_gp(idx)
      if (ratio <= 0.d0 .or. ratio >= 1.d0) cycle
      radius_factor = ratio / (1.d0 - ratio)
      radius_phys = r_e_current * sqrt_kappa * radius_factor**s_pwr
      if (radius_phys <= 0.d0) cycle

      g_value = field(idx) * radius_phys**(ell + 1)
      x_value = 1.d0 / max(radius_phys*radius_phys, tiny(1.d0))

      count = count + 1
      if (count > sample_points) exit

      radius_vals(count) = radius_phys
      field_vals(count)  = field(idx)
      x_vals(count)      = x_value

      sum_x  = sum_x  + x_value
      sum_xx = sum_xx + x_value * x_value
      sum_y  = sum_y  + g_value
      sum_xy = sum_xy + x_value * g_value
    end do

    if (count < 3) return

    denom = count * sum_xx - sum_x * sum_x
    if (abs(denom) < 1.d-20) return

    coeff_leading = (sum_y * sum_xx - sum_x * sum_xy) / denom
    if (coeff_leading >= 0.d0) then
      coeff_leading = 0.d0
      coeff_next    = 0.d0
      return
    end if

    res_sum_x  = 0.d0
    res_sum_xx = 0.d0
    res_sum_y  = 0.d0
    res_sum_xy = 0.d0
    res_count  = 0

    do idx = 1, count
      residual = field_vals(idx) - coeff_leading / radius_vals(idx)**(ell + 1)
      g_value  = residual * radius_vals(idx)**(ell + 3)
      if (abs(g_value) <= tiny(1.d0)) cycle

      res_count  = res_count + 1
      res_sum_x  = res_sum_x  + x_vals(idx)
      res_sum_xx = res_sum_xx + x_vals(idx) * x_vals(idx)
      res_sum_y  = res_sum_y  + g_value
      res_sum_xy = res_sum_xy + x_vals(idx) * g_value
    end do

    if (res_count >= 2) then
      res_denom = res_count * res_sum_xx - res_sum_x * res_sum_x
      if (abs(res_denom) > 1.d-20) then
        coeff_next = (res_sum_y * res_sum_xx - res_sum_x * res_sum_xy) / res_denom
      end if
    end if

    success = .true.
  end subroutine composite_richardson

end module miscellaneous_mod
