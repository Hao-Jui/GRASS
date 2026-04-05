! ============================================================================
! Piecewise Polytropic Equation of State Generator
! ============================================================================

module physical_constants
  implicit none
  private
  public :: clite, MB, rho_uni, gravi
  
  real(8), parameter :: clite = 2.99792458d10      ! Speed of light [cm/s]
  real(8), parameter :: MB = 1.6749286d-24         ! Baryon mass [g]
  real(8), parameter :: rho_uni = 1.d0 / 1.61930347d-18  ! Unit density conversion
  real(8), parameter :: gravi = 6.67408d-20        ! Gravitational constant [cgs]
end module physical_constants

! ============================================================================

module crust_parameters
  implicit none
  private
  public :: crust_Gamma1, crust_Gamma2, crust_Gamma3, crust_Gamma4
  public :: crust_K1, crust_K2, crust_K3, crust_K4
  public :: crust_division1, crust_division2, crust_division3
  
  ! Adiabatic indices for crust regions
  real(8), parameter :: crust_Gamma1 = 1.58425d0
  real(8), parameter :: crust_Gamma2 = 1.28733d0
  real(8), parameter :: crust_Gamma3 = 0.62223d0
  real(8), parameter :: crust_Gamma4 = 1.35692d0
  
  ! Polytropic constants for crust regions
  real(8), parameter :: crust_K1 = 6.80110d-9
  real(8), parameter :: crust_K2 = 1.06186d-6
  real(8), parameter :: crust_K3 = 5.32697d1
  real(8), parameter :: crust_K4 = 3.99874d-8
  
  ! Density division boundaries (log10 scale)
  real(8), parameter :: crust_division1 = log10(2.44034d7)
  real(8), parameter :: crust_division2 = log10(3.78358d11)
  real(8), parameter :: crust_division3 = log10(2.62780d12)
end module crust_parameters

! ============================================================================

module core_parameters
  implicit none
  private
  public :: eos_name, num_points
  public :: dens_division1, dens_division2
  public :: Gamma1, Gamma2, Gamma3, log_p1
  public :: rho_match, K1, K2, K3
  public :: aa, aaa, dd
  
  character(len=128) :: eos_name = "Zdunik"
  integer, parameter :: num_points = 5931
  
  ! Core density divisions (log10 scale)
  real(8), parameter :: dens_division1 = 14.7d0
  real(8), parameter :: dens_division2 = 15.0d0
  
  ! Core EOS parameters (set by setup_eos)
  real(8) :: Gamma1, Gamma2, Gamma3, log_p1
  real(8) :: rho_match
  real(8) :: K1, K2, K3
  
  ! Working variables
  real(8) :: aa                    ! Current integration constant
  real(8) :: aaa(7)                ! Integration constants for all regions
  real(8) :: dd(6)                 ! Density boundaries
end module core_parameters

! ============================================================================

module eos_database
  implicit none
  private
  public :: setup_eos
  
contains

  subroutine setup_eos(name, Gamma1, Gamma2, Gamma3, log_p1, success)
    use core_parameters, only: eos_name
    character(len=*), intent(in) :: name
    real(8), intent(out) :: Gamma1, Gamma2, Gamma3, log_p1
    logical, intent(out) :: success
    
    success = .true.
    
    select case (trim(adjustl(name)))
      case ('PAL6')
        Gamma1 = 2.227d0; Gamma2 = 2.189d0; Gamma3 = 2.159d0; log_p1 = 34.380d0
      case ('SLy')
        Gamma1 = 3.005d0; Gamma2 = 2.988d0; Gamma3 = 2.851d0; log_p1 = 34.384d0
      case ('APR1')
        Gamma1 = 2.442d0; Gamma2 = 3.256d0; Gamma3 = 2.908d0; log_p1 = 33.943d0
      case ('APR2')
        Gamma1 = 2.643d0; Gamma2 = 3.014d0; Gamma3 = 2.945d0; log_p1 = 34.126d0
      case ('APR3')
        Gamma1 = 3.166d0; Gamma2 = 3.573d0; Gamma3 = 3.281d0; log_p1 = 34.392d0
      case ('APR4')
        Gamma1 = 2.830d0; Gamma2 = 3.445d0; Gamma3 = 3.348d0; log_p1 = 34.269d0
      case ('FPS')
        Gamma1 = 2.985d0; Gamma2 = 2.863d0; Gamma3 = 2.600d0; log_p1 = 34.283d0
      case ('WFF1')
        Gamma1 = 2.519d0; Gamma2 = 3.791d0; Gamma3 = 3.660d0; log_p1 = 34.031d0
      case ('WFF2')
        Gamma1 = 2.888d0; Gamma2 = 3.475d0; Gamma3 = 3.517d0; log_p1 = 34.233d0
      case ('WFF3')
        Gamma1 = 3.329d0; Gamma2 = 2.952d0; Gamma3 = 2.589d0; log_p1 = 34.283d0
      case ('BBB2')
        Gamma1 = 3.418d0; Gamma2 = 2.835d0; Gamma3 = 2.832d0; log_p1 = 34.331d0
      case ('ENG')
        Gamma1 = 3.514d0; Gamma2 = 3.130d0; Gamma3 = 3.168d0; log_p1 = 34.437d0
      case ('MPA1')
        Gamma1 = 3.446d0; Gamma2 = 3.572d0; Gamma3 = 2.887d0; log_p1 = 34.495d0
      case ('MS1')
        Gamma1 = 3.224d0; Gamma2 = 3.033d0; Gamma3 = 1.325d0; log_p1 = 34.858d0
      case ('MS2')
        Gamma1 = 2.447d0; Gamma2 = 2.184d0; Gamma3 = 1.855d0; log_p1 = 34.605d0
      case ('MS1b')
        Gamma1 = 3.456d0; Gamma2 = 3.011d0; Gamma3 = 1.425d0; log_p1 = 34.855d0
      case ('PS')
        Gamma1 = 2.216d0; Gamma2 = 1.640d0; Gamma3 = 2.365d0; log_p1 = 34.671d0
      case ('BGN1H1')
        Gamma1 = 3.258d0; Gamma2 = 1.472d0; Gamma3 = 2.464d0; log_p1 = 34.623d0
      case ('GNH3')
        Gamma1 = 2.664d0; Gamma2 = 2.194d0; Gamma3 = 2.304d0; log_p1 = 34.648d0
      case ('H1')
        Gamma1 = 2.595d0; Gamma2 = 1.845d0; Gamma3 = 1.897d0; log_p1 = 34.564d0
      case ('H2')
        Gamma1 = 2.775d0; Gamma2 = 1.855d0; Gamma3 = 1.858d0; log_p1 = 34.617d0
      case ('H3')
        Gamma1 = 2.787d0; Gamma2 = 1.951d0; Gamma3 = 1.901d0; log_p1 = 34.646d0
      case ('H4')
        Gamma1 = 2.909d0; Gamma2 = 2.246d0; Gamma3 = 2.144d0; log_p1 = 34.669d0
      case ('H5')
        Gamma1 = 2.793d0; Gamma2 = 1.974d0; Gamma3 = 1.915d0; log_p1 = 34.609d0
      case ('H7')
        Gamma1 = 2.621d0; Gamma2 = 2.048d0; Gamma3 = 2.006d0; log_p1 = 34.559d0
      case ('PCL2')
        Gamma1 = 2.554d0; Gamma2 = 1.880d0; Gamma3 = 1.977d0; log_p1 = 34.507d0
      case ('ALF1')
        Gamma1 = 2.013d0; Gamma2 = 3.389d0; Gamma3 = 2.033d0; log_p1 = 34.055d0
      case ('ALF2')
        Gamma1 = 4.070d0; Gamma2 = 2.411d0; Gamma3 = 1.890d0; log_p1 = 34.616d0
      case ('ALF3')
        Gamma1 = 2.883d0; Gamma2 = 2.653d0; Gamma3 = 1.952d0; log_p1 = 34.283d0
      case ('ALF4')
        Gamma1 = 3.009d0; Gamma2 = 3.438d0; Gamma3 = 1.803d0; log_p1 = 34.314d0
      case ('HB')
        Gamma1 = 3.000d0; Gamma2 = 3.000d0; Gamma3 = 3.000d0; log_p1 = 34.40364d0
      case ('H')
        Gamma1 = 3.000d0; Gamma2 = 3.000d0; Gamma3 = 3.000d0; log_p1 = 34.50364d0
      case ('B')
        Gamma1 = 3.000d0; Gamma2 = 3.000d0; Gamma3 = 3.000d0; log_p1 = 34.30364d0
      case ('15H')
        Gamma1 = 3.000d0; Gamma2 = 3.000d0; Gamma3 = 3.000d0; log_p1 = 34.70364d0
      case ('125H')
        Gamma1 = 3.000d0; Gamma2 = 3.000d0; Gamma3 = 3.000d0; log_p1 = 34.60364d0
      case ('Zdunik')
        Gamma1 = 3.000d0; Gamma2 = 1.050d0; Gamma3 = 5.100d0; log_p1 = 34.384d0
      case ('PSt')
        Gamma1 = 2.216d0; Gamma2 = 1.300d0; Gamma3 = 4.000d0; log_p1 = 34.671d0
      case default
        success = .false.
    end select
  end subroutine setup_eos

end module eos_database

! ============================================================================

module eos_calculator
  use physical_constants
  use crust_parameters
  use core_parameters
  implicit none
  private
  public :: calculate_matching_density, calculate_K_continuity
  public :: calculate_integration_constants, get_eos_parameters
  
contains

  subroutine calculate_matching_density()
    ! Calculate the matching density between crust and core
    rho_match = (log_p1 - 2.d0 * log10(clite) &
                 - (log10(crust_K4) + crust_Gamma4 * crust_division3) &
                 - Gamma1 * dens_division1 + crust_Gamma4 * crust_division3) &
                / (crust_Gamma4 - Gamma1)
  end subroutine calculate_matching_density

  subroutine calculate_K_continuity()
    ! Ensure K constants are continuous across boundaries
    K1 = crust_K4 * (1.d1**rho_match)**(crust_Gamma4 - Gamma1)
    K2 = K1 * (1.d1**dens_division1)**(Gamma1 - Gamma2)
    K3 = K2 * (1.d1**dens_division2)**(Gamma2 - Gamma3)
  end subroutine calculate_K_continuity

  subroutine calculate_integration_constants()
    ! Calculate integration constants for all regions
    dd(1) = 2.44034d7
    dd(2) = 3.78358d11
    dd(3) = 2.62780d12
    dd(4) = 1.d1**rho_match
    dd(5) = 1.d1**dens_division1
    dd(6) = 1.d1**dens_division2

    aaa(1) = 0.d0
    aaa(2) = compute_next_constant(aaa(1), dd(1), crust_K1, crust_Gamma1, crust_K2, crust_Gamma2)
    aaa(3) = compute_next_constant(aaa(2), dd(2), crust_K2, crust_Gamma2, crust_K3, crust_Gamma3)
    aaa(4) = compute_next_constant(aaa(3), dd(3), crust_K3, crust_Gamma3, crust_K4, crust_Gamma4)
    aaa(5) = compute_next_constant(aaa(4), dd(4), crust_K4, crust_Gamma4, K1, Gamma1)
    aaa(6) = compute_next_constant(aaa(5), dd(5), K1, Gamma1, K2, Gamma2)
    aaa(7) = compute_next_constant(aaa(6), dd(6), K2, Gamma2, K3, Gamma3)
  end subroutine calculate_integration_constants

  function compute_next_constant(a_prev, d, K_prev, G_prev, K_next, G_next) result(a_next)
    real(8), intent(in) :: a_prev, d, K_prev, G_prev, K_next, G_next
    real(8) :: a_next
    
    a_next = ((1.d0 + a_prev) * d + K_prev * d**G_prev / (G_prev - 1.d0)) / d - 1.d0 &
             - K_next * d**(G_next - 1.d0) / (G_next - 1.d0)
  end function compute_next_constant

  subroutine get_eos_parameters(log_rho, G_practical, K_practical)
    ! Determine which EOS parameters to use based on density
    real(8), intent(in) :: log_rho
    real(8), intent(out) :: G_practical, K_practical
    
    if (log_rho <= crust_division1) then
      aa = aaa(1); G_practical = crust_Gamma1; K_practical = crust_K1
    elseif (log_rho <= crust_division2) then
      aa = aaa(2); G_practical = crust_Gamma2; K_practical = crust_K2
    elseif (log_rho <= crust_division3) then
      aa = aaa(3); G_practical = crust_Gamma3; K_practical = crust_K3
    elseif (log_rho <= rho_match) then
      aa = aaa(4); G_practical = crust_Gamma4; K_practical = crust_K4
    elseif (log_rho <= dens_division1) then
      aa = aaa(5); G_practical = Gamma1; K_practical = K1
    elseif (log_rho <= dens_division2) then
      aa = aaa(6); G_practical = Gamma2; K_practical = K2
    else
      aa = aaa(7); G_practical = Gamma3; K_practical = K3
    endif
  end subroutine get_eos_parameters

end module eos_calculator

! ============================================================================

module output_module
  use physical_constants
  use crust_parameters
  use core_parameters
  implicit none
  private
  public :: print_summary, write_eos_table
  
contains

  subroutine print_summary()
    write(*,'(A)') '========================================'
    write(*,'(A,A)') 'EOS: ', trim(eos_name)
    write(*,'(A)') '========================================'
    write(*,'(A,2ES18.9)') 'Matching density: ', 10.d0**rho_match, 10.d0**rho_match / rho_uni
    write(*,*)
    write(*,'(A)') 'Density boundaries (g/cm^3):'
    write(*,'(10ES18.9)') dd
    write(*,*)
    write(*,'(A)') 'Integration constants:'
    write(*,'(10ES18.9)') aaa
    write(*,*)
    write(*,'(A)') 'Polytropic constants:'
    write(*,'(10ES18.9)') crust_K1, crust_K2, crust_K3, crust_K4, K1, K2, K3
    write(*,*)
    write(*,'(A)') 'Parameters for SACRA:'
    write(*,'(20ES18.9)') crust_Gamma1, crust_Gamma2, crust_Gamma3, crust_Gamma4, &
                          Gamma1, Gamma2, Gamma3, dd * 1.61930347d-18, &
                          crust_K1 * 1.61930347d-18**(1.d0 - crust_Gamma1)
    write(*,'(A)') '========================================'
  end subroutine print_summary

  subroutine write_eos_table()
    use eos_calculator
    integer :: unit, ios, ll
    real(8) :: log_rho, rho, step
    real(8) :: G_practical, K_practical, prs, hh, eps, cspeed
    character(len=256) :: filename
    
    ! Create output filename
    filename = '../' // trim(adjustl(eos_name)) // '.dat'
    
    ! Open file
    open(newunit=unit, file=filename, status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,I0)') 'ERROR: Cannot open output file. IOSTAT = ', ios
      stop
    endif
    
    write(unit,*) num_points
    
    ! Set up density grid
    log_rho = 3.d0
    step = (log10(5.d16) - log_rho) / dble(num_points - 1)
    
    ! Calculate EOS for each density point
    do ll = 1, num_points
      rho = 10.d0**log_rho
      
      call get_eos_parameters(log_rho, G_practical, K_practical)
      
      ! Calculate thermodynamic quantities
      prs = K_practical * rho**G_practical
      eps = rho * (1.d0 + aa) + prs / (G_practical - 1.d0)
      cspeed = sqrt(G_practical * prs / (eps + prs))
      hh = (1.d0 + aa) + K_practical * rho**(G_practical - 1.d0) * G_practical / (G_practical - 1.d0)
      
      write(unit,'(5ES18.11)') eps, prs * clite**2, hh, rho / MB, cspeed**2
      log_rho = log_rho + step
    enddo
    
    close(unit)
    write(*,'(A,A)') 'EOS table written to: ', trim(filename)
  end subroutine write_eos_table

end module output_module

! ============================================================================

program piecewise_polytropic_eos
  use eos_database
  use eos_calculator
  use output_module
  use core_parameters
  implicit none
  
  logical :: success
  
  ! Initialize EOS parameters
  call setup_eos(eos_name, Gamma1, Gamma2, Gamma3, log_p1, success)
  if (.not. success) then
    write(*,'(A,A)') 'ERROR: Unknown EOS name: ', trim(eos_name)
    stop
  endif
  
  ! Calculate matching conditions
  call calculate_matching_density()
  call calculate_K_continuity()
  call calculate_integration_constants()
  
  ! Print summary
  call print_summary()
  
  ! Generate and write EOS table
  call write_eos_table()
  
  write(*,*)
  write(*,'(A)') 'Program completed successfully!'
  
end program piecewise_polytropic_eos
