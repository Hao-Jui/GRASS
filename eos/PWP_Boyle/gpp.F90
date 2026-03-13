!-*-f90-*-
module physical_constants
  implicit none
  real(8), parameter :: clite = 2.99792458d10  ! Speed of light [cm/s]
  real(8), parameter :: MB = 1.6749286d-24     ! Neutron mass [g]
end module physical_constants

module crust_parameters
  implicit none
  
  ! Crust density divisions (log10 scale)
  real(8), parameter :: crust_devision1 = log10(6.285d5)
  real(8), parameter :: crust_devision2 = log10(1.826d8)
  real(8), parameter :: crust_devision3 = log10(3.350d11)
  real(8), parameter :: crust_devision4 = log10(5.317d11)
  
  ! Crust polytropic parameters
  real(8), parameter :: crust_Gamma1 = 1.611d0
  real(8), parameter :: crust_Gamma2 = 1.440d0
  real(8), parameter :: crust_Gamma3 = 1.269d0
  real(8), parameter :: crust_Gamma4 = -1.841d0
  real(8), parameter :: crust_Gamma5 = 1.382d0
  
  real(8), parameter :: crust_K1 = 5.214d-9
  real(8), parameter :: crust_K2 = 5.726d-8
  real(8), parameter :: crust_K3 = 1.662d-6
  real(8), parameter :: crust_K4 = -7.957d29
  real(8), parameter :: crust_K5 = 1.746d-8
  
  real(8), parameter :: crust_Lambda1 = 0.0d0
  real(8), parameter :: crust_Lambda2 = -1.354d0
  real(8), parameter :: crust_Lambda3 = -6.025d3
  real(8), parameter :: crust_Lambda4 = 1.193d9
  real(8), parameter :: crust_Lambda5 = 7.077d8
  
  real(8), parameter :: crust_a1 = 0.0d0
  real(8), parameter :: crust_a2 = -1.861d-5
  real(8), parameter :: crust_a3 = -5.278d-4
  real(8), parameter :: crust_a4 = 1.035d-2
  real(8), parameter :: crust_a5 = 8.208d-3
end module crust_parameters

module core_parameters
  implicit none
  
  character(len=128) :: eos_name = "APR4"
  integer, parameter :: num_points = 5931
  
  ! Core density divisions (log10 scale)
  real(8), parameter :: dens_division1 = 14.87d0
  real(8), parameter :: dens_division2 = 14.99d0
  
  ! Core polytropic parameters (to be read from table)
  real(8) :: Gamma1, Gamma2, Gamma3
  real(8) :: K1, K2, K3
  real(8) :: Lambda1, Lambda2, Lambda3
  real(8) :: rho_match
  real(8) :: aaa(3)  ! Energy shift parameters for core regions
end module core_parameters

module eos_setup
  use physical_constants
  use crust_parameters
  use core_parameters
  implicit none
  
contains
  
  subroutine read_eos_table()
    implicit none

    ! Load tabulated EOS parameters for the chosen eos_name
    call setup_table()

    write(*,*) 'EOS table setup complete for ', trim(eos_name)
  end subroutine read_eos_table

  subroutine setup_table()
    implicit none
    character(len=20) :: valid_eos(13)
    integer :: i
    logical :: eos_found

    ! Define valid EOS names
    data valid_eos / 'BHF', 'FPS', 'H4', 'QHC19', 'RS', 'KDE0V', &
                     'MS1', 'MS1b', 'APR4', 'MPA1', 'SKOP', 'SKI4', 'SKI2' /
    
    eos_found = .false.
    
    select case (trim(eos_name))
      
      case ('BHF')
        K1 = 10.d0**(-35.016d0)
        Gamma1 = 3.284d0
        Gamma2 = 2.774d0
        Gamma3 = 2.618d0
        eos_found = .true.
        
      case ('FPS')
        K1 = 10.d0**(-32.985d0)
        Gamma1 = 3.147d0
        Gamma2 = 2.652d0
        Gamma3 = 2.199d0
        eos_found = .true.
        
      case ('H4')
        K1 = 10.d0**(-23.310d0)
        Gamma1 = 2.514d0
        Gamma2 = 2.333d0
        Gamma3 = 1.562d0
        eos_found = .true.
        
      case ('QHC19')
        K1 = 10.d0**(-36.879d0)
        Gamma1 = 3.419d0
        Gamma2 = 2.760d0
        Gamma3 = 2.017d0
        eos_found = .true.
        
      case ('RS')
        K1 = 10.d0**(-25.150d0)
        Gamma1 = 2.636d0
        Gamma2 = 2.677d0
        Gamma3 = 2.0647d0
        eos_found = .true.
        
      case ('KDE0V')
        K1 = 10.d0**(-30.250d0)
        Gamma1 = 2.967d0
        Gamma2 = 2.835d0
        Gamma3 = 2.803d0
        eos_found = .true.
        
      case ('MS1')
        K1 = 10.d0**(-30.170d0)
        Gamma1 = 2.998d0
        Gamma2 = 2.123d0
        Gamma3 = 1.955d0
        eos_found = .true.
        
      case ('MS1b')
        K1 = 10.d0**(-33.774d0)
        Gamma1 = 3.241d0
        Gamma2 = 2.136d0
        Gamma3 = 1.963d0
        eos_found = .true.
        
      case ('APR4')
        K1 = 10.d0**(-33.210d0)
        Gamma1 = 3.169d0
        Gamma2 = 3.452d0
        Gamma3 = 3.310d0
        eos_found = .true.
        
      case ('MPA1')
        K1 = 10.d0**(-40.301d0)
        Gamma1 = 3.662d0
        Gamma2 = 3.057d0
        Gamma3 = 2.298d0
        eos_found = .true.
        
      case ('SKOP')
        K1 = 10.d0**(-26.089d0)
        Gamma1 = 2.693d0
        Gamma2 = 2.660d0
        Gamma3 = 2.579d0
        eos_found = .true.
        
      case ('SKI4')
        K1 = 10.d0**(-31.008d0)
        Gamma1 = 3.029d0
        Gamma2 = 2.759d0
        Gamma3 = 2.651d0
        eos_found = .true.
        
      case ('SKI2')
        K1 = 10.d0**(-24.202d0)
        Gamma1 = 2.575d0
        Gamma2 = 2.639d0
        Gamma3 = 2.656d0
        eos_found = .true.
        
      case default
        write(*,'(A)') 'ERROR: Unknown equation of state: ' // trim(eos_name)
        write(*,'(A)') 'Valid options are:'
        do i = 1, size(valid_eos)
          write(*,'(A,A)') '  - ', trim(valid_eos(i))
        end do
        stop 'Execution terminated: invalid EOS name'
        
    end select

  end subroutine setup_table
  
  subroutine calculate_matching_parameters()
    implicit none
    
    ! Calculate matching density between crust and core
    rho_match = log10(K1 * Gamma1 / crust_K5 / crust_Gamma5) / (crust_Gamma5 - Gamma1)
    
    ! Calculate K parameters for continuity
    K2 = K1 * (1.d1**dens_division1)**(Gamma1 - Gamma2) * Gamma1 / Gamma2
    K3 = K2 * (1.d1**dens_division2)**(Gamma2 - Gamma3) * Gamma2 / Gamma3
    
    ! Calculate Lambda parameters for continuity
    Lambda1 = crust_Lambda5 + (1.0d0 - crust_Gamma5/Gamma1) * crust_K5 * (1.d1**rho_match)**crust_Gamma5
    Lambda2 = Lambda1 + (1.0d0 - Gamma1/Gamma2) * K1 * (1.d1**dens_division1)**Gamma1
    Lambda3 = Lambda2 + (1.0d0 - Gamma2/Gamma3) * K2 * (1.d1**dens_division2)**Gamma2
    
    ! Calculate energy shift parameters
    aaa(1) = crust_a5 + crust_Gamma5 * (Gamma1 - crust_Gamma5) / (crust_Gamma5 - 1.0d0) / &
             (Gamma1 - 1.0d0) * crust_K5 * (1.d1**rho_match)**(crust_Gamma5 - 1.0d0)
    aaa(2) = aaa(1) + Gamma1 * (Gamma2 - Gamma1) / (Gamma1 - 1.0d0) / (Gamma2 - 1.0d0) * &
             K1 * (1.d1**dens_division1)**(Gamma1 - 1.0d0)
    aaa(3) = aaa(2) + Gamma2 * (Gamma3 - Gamma2) / (Gamma2 - 1.0d0) / (Gamma3 - 1.0d0) * &
             K2 * (1.d1**dens_division2)**(Gamma2 - 1.0d0)
  end subroutine calculate_matching_parameters
  
  subroutine print_eos_summary()
    implicit none
    
    write(*,*) '==== EOS Configuration Summary ===='
    write(*,*) 'Density partitions [g/cm^3]:'
    write(*,"(8es15.6)") 0.d0, 1.d1**crust_devision1, 1.d1**crust_devision2, &
                         1.d1**crust_devision3, 1.d1**crust_devision4, &
                         1.d1**rho_match, 1.d1**dens_division1, 1.d1**dens_division2
    write(*,*)
    write(*,*) 'Energy shift parameters (a):'
    write(*,"(8es15.6)") crust_a1, crust_a2, crust_a3, crust_a4, crust_a5, aaa
    write(*,*)
    write(*,*) 'Polytropic K parameters:'
    write(*,"(8es15.6)") crust_K1, crust_K2, crust_K3, crust_K4, crust_K5, K1, K2, K3
    write(*,*)
    write(*,*) 'Lambda parameters:'
    write(*,"(8es15.6)") crust_Lambda1, crust_Lambda2, crust_Lambda3, crust_Lambda4, &
                         crust_Lambda5, Lambda1, Lambda2, Lambda3
    write(*,*)
    write(*,*) 'Gamma parameters:'
    write(*,"(8es15.6)") crust_Gamma1, crust_Gamma2, crust_Gamma3, crust_Gamma4, &
                         crust_Gamma5, Gamma1, Gamma2, Gamma3
    write(*,*) '==================================='
    write(*,*)
  end subroutine print_eos_summary
  
end module eos_setup

module eos_calculations
  use physical_constants
  use crust_parameters
  use core_parameters
  implicit none
  
contains
  
  subroutine get_region_parameters(log_rho, aa, lmd, G_practical, K_practical)
    implicit none
    real(8), intent(in) :: log_rho
    real(8), intent(out) :: aa, lmd, G_practical, K_practical
    
    if (log_rho <= crust_devision1) then
      aa = crust_a1
      lmd = crust_Lambda1
      G_practical = crust_Gamma1
      K_practical = crust_K1
    elseif (log_rho <= crust_devision2) then
      aa = crust_a2
      lmd = crust_Lambda2
      G_practical = crust_Gamma2
      K_practical = crust_K2
    elseif (log_rho <= crust_devision3) then
      aa = crust_a3
      lmd = crust_Lambda3
      G_practical = crust_Gamma3
      K_practical = crust_K3
    elseif (log_rho <= crust_devision4) then
      aa = crust_a4
      lmd = crust_Lambda4
      G_practical = crust_Gamma4
      K_practical = crust_K4
    elseif (log_rho <= rho_match) then
      aa = crust_a5
      lmd = crust_Lambda5
      G_practical = crust_Gamma5
      K_practical = crust_K5
    elseif (log_rho <= dens_division1) then
      aa = aaa(1)
      lmd = Lambda1
      G_practical = Gamma1
      K_practical = K1
    elseif (log_rho <= dens_division2) then
      aa = aaa(2)
      lmd = Lambda2
      G_practical = Gamma2
      K_practical = K2
    else
      aa = aaa(3)
      lmd = Lambda3
      G_practical = Gamma3
      K_practical = K3
    endif
  end subroutine get_region_parameters
  
  subroutine calculate_thermodynamic_quantities(rho, aa, lmd, G_practical, K_practical, &
                                                eps, prs, cspead, hh)
    implicit none
    real(8), intent(in) :: rho, aa, lmd, G_practical, K_practical
    real(8), intent(out) :: eps, prs, cspead, hh
    real(8) :: rho_gamma, denominator
    
    rho_gamma = rho**G_practical
    
    ! Pressure
    prs = K_practical * rho_gamma + lmd
    
    ! Energy density
    eps = rho * (1.0d0 + aa) + K_practical * rho_gamma / (G_practical - 1.0d0) - lmd
    
    ! Sound speed squared
    denominator = (1.0d0 + aa) + K_practical * G_practical * rho**(G_practical - 1.0d0) / &
                  (G_practical - 1.0d0)
    cspead = K_practical * G_practical * rho**(G_practical - 1.0d0) / denominator
    
    ! Specific enthalpy
    hh = (eps + prs) / rho
  end subroutine calculate_thermodynamic_quantities
  
end module eos_calculations

program piecewise_polytropic_eos
  use physical_constants
  use crust_parameters
  use core_parameters
  use eos_setup
  use eos_calculations
  implicit none
  
  real(8) :: log_rho, rho, step
  real(8) :: aa, lmd, G_practical, K_practical
  real(8) :: eps, prs, cspead, hh
  integer :: i, unit, ios
  character(len=256) :: output_filename
  
  ! Initialize EOS
  call read_eos_table()
  call calculate_matching_parameters()
  call print_eos_summary()
  
  ! Open output file
  output_filename = '../' // trim(adjustl(eos_name)) // '_7pc.dat'
  open(newunit=unit, file=output_filename, status='replace', action='write', iostat=ios)
  if (ios /= 0) then
    write(*,*) 'Error: Cannot open output file ', trim(output_filename)
    write(*,*) 'IOSTAT = ', ios
    stop
  endif
  
  write(*,*) 'Writing EOS table to: ', trim(output_filename)
  write(unit,*) num_points
  
  ! Generate EOS table
  step = (log10(7.d15) - 3.0d0) / dble(num_points)
  
  do i = 0, num_points - 1
    log_rho = 3.0d0 + i * step
    rho = 10.0d0**log_rho
    
    ! Get parameters for current density region
    call get_region_parameters(log_rho, aa, lmd, G_practical, K_practical)
    
    ! Calculate thermodynamic quantities
    call calculate_thermodynamic_quantities(rho, aa, lmd, G_practical, K_practical, &
                                           eps, prs, cspead, hh)
    
    ! Write to file: energy density, pressure, enthalpy, baryon density, sound speed squared
    write(unit,"(5es18.11)") eps, prs * clite**2, hh, rho / MB, cspead
  enddo
  
  close(unit)
  write(*,*) 'EOS table generation complete!'
  
end program piecewise_polytropic_eos
