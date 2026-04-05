module exporter_mod
  use precision_mod, only: wp
  implicit none
  private

  public :: initial_data_for_spec, initial_data_for_sacra_aei

contains

  subroutine initial_data_for_spec(file_name, var1, var2, var3, var4, var5, var6)
    use para_mod, only : s_gp, mu, SDIV, MDIV, l_uni, ang_mom, mass_0, mass, &
                         r_e, KAPPA, Omega_c, MSUN
    implicit none
    character(*), intent(in) :: file_name
    real(wp), intent(in) :: var1(:,:), var2(:,:), var3(:,:), var4(:,:), var5(:,:), var6(:,:)
    integer :: s, m, unit, ios
    real(wp), parameter :: K_KM = 218.04217865726338_wp
    real(wp) :: r_s

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
                            Omega_c * ( sqrt(K_km) / sqrt(KAPPA) ), 0.0d0
    do s = 1, SDIV
      r_s = r_e * sqrt(KAPPA) / 1.d5 / sqrt(K_km) * s_gp(s) / (1.d0 - s_gp(s))
      do m = 1, MDIV
        write(unit,"(4es22.12)") r_s, mu(m), var1(s,m), var2(s,m), &
                                var3(s,m), var4(s,m), var5(s,m), var6(s,m)
      end do
    end do
    close(unit)
  end subroutine initial_data_for_spec

  subroutine initial_data_for_sacra_aei(filename)
    use para_mod, only: SDIV, MDIV, s_pwr, r_e, KAPPA, energy, r_ratio, Omega_e, Omega_c, &
                        C, alpha, gama, rho, ww, pressure, KSCALE, enthalpy, enthalpy_min, &
                        velocity_sq, omg, sphi, B_coup, s_gp, mu, MB
    use eos_mod, only: n0_at_e
    implicit none
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
            energy(1,1)/(C*C*KSCALE), r_ratio, Omega_e* (C/sqrt(KAPPA)), Omega_c* (C/sqrt(KAPPA))

    do s = 1, SDIV
      do m = 1, MDIV
        if (enthalpy(s,m) > enthalpy_min) then
          rho_0_val = n0_at_e(energy(s,m)) * MB
        else
          rho_0_val = 0.0_wp
        end if

        write(unit, fmt=data_fmt) s_gp(s), mu(m), alpha(s,m), gama(s,m), rho(s,m), &
          ww(s,m) * (C/sqrt(KAPPA)), pressure(s,m)/KSCALE, energy(s,m)/(C*C*KSCALE), &
          enthalpy(s,m), rho_0_val, velocity_sq(s,m), omg(s,m) * (C/sqrt(KAPPA)), &
          sphi(s,m) * sqrt(B_coup)
      end do
    end do

    close(unit)
  end subroutine initial_data_for_sacra_aei

end module exporter_mod
