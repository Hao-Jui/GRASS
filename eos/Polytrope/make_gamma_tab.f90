!===============================================================================
! Program: gamma_tab
! Purpose: Generate polytropic equation of state table in geometric units
! Author: Hao-Jui Kuan
!
! Description:
!   Computes a polytropic EOS table with P = K * rho^Gamma
!   Outputs: energy density, pressure, specific enthalpy, number density
!   in geometric units suitable for relativistic hydrodynamics
!===============================================================================

program gamma_tab
  implicit none

  ! Physical constants (CGS units)
  real(8), parameter :: MB = 1.6749286d-24          ! Baryon mass [g]
  real(8), parameter :: clite = 2.99792458d10       ! Speed of light [cm/s]

  ! Unit conversion factors (CGS to geometric)
  real(8), parameter :: l_uni = 1.4769994423016508d0    ! Length [km]
  real(8), parameter :: rho_uni = 1.61930347d-18        ! Density [g/cm^3]
  real(8), parameter :: prs_uni = 1.80171810d-39        ! Pressure [dyne/cm^2]

  ! EOS parameters
  real(8), parameter :: Gamma = 1.0d0 + 1.d0 / 1.0d0    ! Polytropic index
  real(8), parameter :: K_km  = 100.d0
  real(8), parameter :: K_geo = 100.d0 !K_km / l_uni**((Gamma-1.d0)*2.d0)

  ! Grid parameters
  real(8), parameter :: log_rho_min = 5.0d0             ! log10(rho_min) [g/cm^3]
  real(8), parameter :: log_rho_max = 16.0d0            ! log10(rho_max) [g/cm^3]
  integer, parameter :: n_points = 10001                ! Number of table points

  ! Variables
  character(len=32) :: output_file = "Gam2_K100.dat"
  real(8) :: prs_geo, rho_geo, ee_geo, hh_geo
  real(8) :: step, log_rho, QB
  integer :: i, ios, unit

  !---------------------------------------------------------------------------
  ! Main computation
  !---------------------------------------------------------------------------

  ! Print header information
  call print_header()

  ! Initialize logarithmic density grid
  step = (log_rho_max - log_rho_min) / real(n_points - 1, 8)
  log_rho = log_rho_min

  ! Open output file
  open(newunit=unit, file="../"//trim(adjustl(output_file)), status='replace', &
     action='write', iostat=ios)

  if (ios /= 0) then
    write(*,'(A)') 'ERROR: Cannot open output file: ' // trim(output_file)
    stop 1
  end if

  ! Write number of points (table header)
  write(unit, '(I8)') n_points

  ! Generate EOS table
  do i = 1, n_points
    ! Calculate density at current grid point
    rho_geo = 10.0d0**log_rho * rho_uni

    ! Compute thermodynamic quantities in geometric units
    prs_geo = K_geo * rho_geo**Gamma
    ee_geo = prs_geo / (Gamma - 1.0d0) + rho_geo
    hh_geo = (ee_geo + prs_geo) / rho_geo

    ! Compute baryon number density
    QB = rho_geo / rho_uni / MB

    ! Write to file: energy density, pressure, specific enthalpy, number density
    write(unit, '(4ES27.19)') ee_geo / rho_uni, prs_geo / prs_uni, &
                  hh_geo, QB

    ! Increment log(density)
    log_rho = log_rho + step
  end do

  ! Close output file
  close(unit)

  ! Print completion message
  call print_footer()

contains

  !---------------------------------------------------------------------------
  ! Subroutine: print_header
  ! Purpose: Display program information and parameters
  !---------------------------------------------------------------------------
  subroutine print_header()
    implicit none

    write(*, '(A)') repeat('=', 70)
    write(*, '(A)') 'POLYTROPIC EOS TABLE GENERATOR'
    write(*, '(A)') repeat('=', 70)
    write(*, '(A)') ''
    write(*, '(A)') 'EOS Parameters:'
    write(*, '(A, F8.3)') '  Polytropic index (Gamma)    : ', Gamma
    write(*, '(A,F18.9)') '  K (geometric units)         : ', K_geo
    write(*, '(A,F18.9)') '  K (length    units)         : ', K_km
    write(*, '(A)') ''
    write(*, '(A)') 'Grid Parameters:'
    write(*, '(A, ES10.3, A, ES10.3, A)') &
      '  Density range               : ', 10.0d0**log_rho_min, &
      ' to ', 10.0d0**log_rho_max, ' g/cm^3'
    write(*, '(A, I8)') '  Number of points            : ', n_points
    write(*, '(A)') ''
    write(*, '(A)') 'Output:'
    write(*, '(A, A)') '  Filename                    : ', trim(output_file)
    write(*, '(A)') '  Columns: eps/rho_uni, P/prs_uni, h, n_B'
    write(*, '(A)') ''
    write(*, '(A)') 'Generating table...'

  end subroutine print_header

  !---------------------------------------------------------------------------
  ! Subroutine: print_footer
  ! Purpose: Display completion message
  !---------------------------------------------------------------------------
  subroutine print_footer()
    implicit none

    write(*, '(A)') ''
    write(*, '(A)') 'Table generation completed successfully!'
    write(*, '(A)') repeat('=', 70)

  end subroutine print_footer

end program gamma_tab
