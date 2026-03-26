module constrain_mod
  use precision_mod, only: wp
  use para_mod, only: SDIV, MDIV, DS, DM, r_e, s_gp, mu, rho, gama, alpha, ww, &
                      energy, pressure, sphi, mphi_r, B_coup, pi, has_scalar, velocity_sq
  use ope_eq_mod, only: laplacian_operator, gradient_vector
  implicit none
  private

  public :: hamiltonian

  type(laplacian_operator) :: op   ! allocated once on first hamiltonian call

contains

  subroutine hamiltonian(hamL2)
    real(wp), intent(out) :: hamL2
    real(wp), dimension(SDIV,MDIV) :: ham
    real(wp), dimension(SDIV,MDIV) :: logPsi4, acoup4
    real(wp), dimension(SDIV,MDIV) :: ricci, ricci_lap, ricci_scal, check
    real(wp), dimension(SDIV,MDIV) :: psi4, gutt, gurr, twist, KK, rhoH, dphidphi, Vphi
    real(wp), dimension(SDIV) :: r_phys
    real(wp), dimension(SDIV,MDIV) :: r2_2d, m1_2d
    integer :: unit, ios, s, m
    type(gradient_vector) :: grad_logPsi4, grad_ww, grad_sphi, grad_alp

    r_phys = r_e * s_gp / max(1.e-30_wp, 1.e0_wp - s_gp)
    r2_2d  = spread(r_phys**2, dim=2, ncopies=MDIV)
    m1_2d  = spread(1.e0_wp - mu**2, dim=1, ncopies=SDIV)
    psi4   = exp(gama - rho)
    logPsi4= gama - rho
    acoup4 = exp(-sphi**2 * B_coup)
    gurr   = exp(-2.e0_wp * alpha)
    twist  = psi4 * r2_2d * m1_2d * ww
    gutt   = -exp(gama + rho)
    rhoH   = (energy + pressure) / (1.e0_wp - velocity_sq) - pressure
    
    call op%init()
    grad_logPsi4 = op%grad(logPsi4)
    grad_ww      = op%grad(ww)
    grad_sphi    = op%grad(sphi)
    grad_alp     = op%grad(alpha)

    ricci_lap = op%lap2(alpha) + 0.5e0_wp * op%laplacian(logPsi4) !op%laplacian(alpha + 0.5e0_wp * logPsi4)
    ricci_scal = 0.25e0_wp * op%scal(grad_logPsi4, grad_logPsi4) + 0.5e0_wp * op%divr(grad_logPsi4%r)
    ricci = -2.e0_wp * gurr * (ricci_lap + ricci_scal )

    check = op%scal(grad_logPsi4, grad_logPsi4) 
    if (has_scalar .and. abs(B_coup) > 1.e-30_wp) then
      dphidphi = op%scal(grad_sphi, grad_sphi) * gurr
      Vphi     = 0.5e0_wp * mphi_r**2 * sphi**2 / B_coup
    else
      dphidphi = 0.e0_wp
      Vphi     = 0.e0_wp
    end if

    !KK = 0.5e0_wp * gurr / psi4 / max(lapsesq, 1.e-30_wp) * op%scal(grad_ww, grad_ww) &
    !   / max(r2_2d, 1.e-30_wp) / max(m1_2d, 1.e-30_wp)
    KK = 0.5e0_wp * gurr * exp(-2.e0_wp * rho) * r2_2d * m1_2d * op%scal(grad_ww, grad_ww)
    ham = ricci + KK - 16.e0_wp * pi * rhoH * acoup4 + gutt * (dphidphi + 2.e0_wp * Vphi)

    hamL2 = 4.e0_wp * pi * r_e**3 * DS * DM * &
            sum(ham**2 * spread(s_gp**2 / max(1.e-30_wp, 1.e0_wp - s_gp)**4, dim=2, ncopies=MDIV))
    !write(*,'(16es8.1)') ricci(:,1)
    !write(*,*) " "
    !write(*,'(16es8.1)') 16.e0_wp * pi * rhoH(:,1) * acoup4(:,1)
    !write(*,*) " "
    !write(*,'(16es8.1)') Vphi(:,1)

    if (.true.) then
      open(newunit=unit, file="./Cont/hamiltonain.dat", status="replace", action="write", iostat=ios)
      write(unit,"(2(i0,2X))") SDIV, MDIV
      if (ios == 0) then
        do s = 1, SDIV
          do m = 1, MDIV
            write(unit, "(3es24.15)") s_gp(s), mu(m), ham(s,m)
          end do
        end do
        close(unit)
      else
        write(*,*) "hamiltonian: failed to open ./Cont/hamiltonain.dat"
      end if
    end if
  end subroutine hamiltonian

end module constrain_mod
