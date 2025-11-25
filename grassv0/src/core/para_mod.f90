module para_mod
  implicit none
  ! -- Theory selection ------------------------------------------------------
  integer, parameter :: THEORY_GR = 0
  integer, parameter :: THEORY_ST = 1
  integer :: active_theory = THEORY_ST

  character(len=20) :: relaxation_scheme = "anderson" ! anderson, newton (much slower per iteration)

  ! -- Running option --------------------------------------------------------
  integer, parameter :: MODE_REGRID  = 1
  integer, parameter :: MODE_DEFAULT = 2

  integer :: run_mode = MODE_REGRID 

  ! -- Rotation configuration ------------------------------------------------
  character(len=20) :: solver_type = "uniform"

  ! -- Solver state ----------------------------------------------------------
  logical :: output = .false.
  logical :: use_shoot_1d = .true.      ! adjust hc while keeping rep constant
  character(len=20) :: FIX1 = "Mb_goal"
  character(len=20) :: FIX2 = "J_goal"

  ! -- Resolutions -----------------------------------------------------------
  integer, parameter :: res  = 300
  integer, parameter :: s_pwr = 1
  integer :: SDIV = 2 * res + 1
  integer :: MDIV = 2 * res + 1

  ! -- Target quantities -----------------------------------------------------
  character(len=128) :: eos_file = "MPA1"
  real(8) :: M_goal   = 1.8d0
  real(8) :: Mb_goal  = 1.6d0
  real(8) :: J_goal   = 1.d0
  real(8) :: chi_goal = 0.63d0
  real(8) :: omc_goal = 30.d0

  real(8) :: B_goal   = 1.d1
  real(8) :: mphi_goal = 0.d-1

  ! -- Rotation-law parameters (advanced modes currently disabled) ----------
  real(8) :: A_diff  = 10.d0
  real(8) :: lambda1 = 2.d0
  real(8) :: lambda2 = 0.5d0
  integer :: uyru_p  = 1
  integer :: uyru_q  = 3

  real(8) :: parA = 1.d0
  real(8) :: parB = 1.d0

  real(8) :: cofA = 30.d0
  real(8) :: cofB = 0.d0
  real(8) :: cofp = 0.8d0
  real(8) :: cofq = 0.9d0

  ! -- Equation of state -----------------------------------------------------
  logical :: phase_transition = .false.
  integer :: num_tab = 0
  integer :: p_at_PT = 0
  real(8), allocatable :: log_p(:), log_e(:), log_h(:), log_n0(:)
  real(8) :: p_center = 0.d0
  real(8) :: h_center = 0.d0
  real(8) :: e_center = 0.d0
  real(8) :: enthalpy_min = 0.d0

  ! -- Grid configuration ----------------------------------------------------
  integer, parameter :: LMAX = 10
  real(8) :: SMAX  = 1.d0 - (1.d-1)**(9.d0 / dble(s_pwr))
  integer, parameter :: RDIV = 1800

  real(8) :: DS = 0.d0
  real(8) :: DM = 0.d0
  real(8), parameter :: s_e = 0.5d0

  real(8), allocatable :: s_gp(:), mu(:), sin_theta(:)

  ! -- Disk helper quantities (kept for compatibility) ----------------------
  logical :: disk_present = .false.
  real(8) :: edge_in = 800.d0
  real(8) :: s_inner = 0.5d0
  real(8) :: j_disk  = 4.5d0
  real(8) :: p_max_disk = 0.d0
  real(8) :: h_max_disk = 0.d0
  integer :: i_isco_p = 1
  integer :: i_isco_m = 1

  ! Fluid
  real(8), allocatable :: pressure(:,:), enthalpy(:,:), velocity_sq(:,:), &
                          energy(:,:), omg(:,:), F_j(:,:)
  real(8), allocatable :: v_plus(:), v_minus(:), V_rr_p(:), V_rr_m(:), sound_speed(:)

  ! Metric
  real(8), allocatable :: gama(:,:), rho(:,:), ww(:,:), alpha(:,:), sphi(:,:)

  ! Scalar field
  real(8) :: B_coup = 0.d0
  real(8) :: mphi_r = 0.d0
  real(8) :: sphi_c = 0.d0
  real(8) :: sphi_m = 0.d0

  real(8) :: B_burn_init = 15.d0
  real(8) :: mphi_burn_seed = 0.01d0
  integer, parameter :: scalar_burn_max_iter = 200
  real(8), parameter :: mphi_burn_threshold = 0.05d0

  ! Bulk properties
  real(8) :: Omega_c = 0.d0
  real(8) :: Omega_e = 0.d0
  real(8) :: Omega_K = 0.d0
  real(8) :: r_e     = 1.d0
  real(8) :: r_ratio = 1.d0
  real(8) :: r_circ  = 0.d0
  real(8) :: ang_mom = 0.d0
  real(8) :: mass    = 0.d0
  real(8) :: mass_0  = 0.d0
  real(8) :: mass_p  = 0.d0
  real(8) :: chi     = 0.d0
  real(8) :: T_kin   = 0.d0
  real(8) :: Fmax_h  = 0.d0
  real(8) :: F_equator_h = 0.d0

  integer :: n_of_relaxation_steps = 0
  ! Multipole information
  real(8) :: M2 = 0.d0
  real(8) :: M4 = 0.d0
  real(8) :: S3 = 0.d0

  ! Spectral helpers
  real(8), allocatable :: f_rho(:,:,:), f_gama(:,:,:)
  real(8), allocatable :: P_2n(:,:), P1_2n_1(:,:), sin_2n_1_theta(:,:)

  ! Timing
  real(8) :: start = 0.d0
  real(8) :: finish = 0.d0

  ! -- Physical constants ----------------------------------------------------
  real(8), parameter :: C    = 2.99792458d10
  real(8), parameter :: G    = 6.67408d-8
  real(8), parameter :: MSUN = 1.98847d33
  real(8), parameter :: MB   = 1.6749286d-24
  real(8), parameter :: pi   = acos(-1.d0)

  real(8), parameter :: accuracy  = 1.d-6
  real(8), parameter :: tov_rmin  = 1.d-15
  real(8), parameter :: KAPPA     = 1.d-15 * C**2 / G
  real(8), parameter :: KSCALE    = KAPPA * G / C**4
  real(8), parameter :: e_surface = 7.8d0 * C**2 * KSCALE
  real(8), parameter :: p_surface = 1.01d8 * KSCALE
  real(8), parameter :: rho_uni   = 1.61930347d-18
  real(8), parameter :: f_uni     = 2.029739818539300d5
  real(8), parameter :: hbar      = 6.582119569d-16
  real(8), parameter :: l_uni     = 1.4769994423016508d0
  real(8), parameter :: n_sat     = 2.7d14
  real(8), parameter :: scalarton = hbar * C / l_uni / 1.d5
contains

  pure function to_lower_str(str) result(out)
    character(*), intent(in) :: str
    character(len(str)) :: out
    integer :: i, code

    out = str
    do i = 1, len(str)
      code = iachar(out(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        out(i:i) = achar(code + 32)
      end if
    end do
  end function to_lower_str

  subroutine initialize_theory()
    select case (active_theory)
    case (THEORY_GR)
      call apply_gr_defaults()
    case (THEORY_ST)
    case default
      stop "initialize_theory: unknown theory mode"
    end select

    solver_type = trim(to_lower_str(adjustl(solver_type)))

    DS   = SMAX / (dble(SDIV) - 1.d0)
    DM   = 1.d0  / (dble(MDIV) - 1.d0)
    s_inner = edge_in**(1.d0 / dble(s_pwr)) / (edge_in**(1.d0 / dble(s_pwr)) + 1.d0)

    call allocate_fields()
  end subroutine initialize_theory

  subroutine apply_gr_defaults()
    B_goal  = 0.d0
    mphi_goal = 0.d0
    B_coup = 0.d0
    mphi_r = 0.d0
    sphi_c = 0.d0
    sphi_m = 0.d0
  end subroutine apply_gr_defaults

  subroutine allocate_fields()
    call deallocate_fields()

    allocate(s_gp(SDIV), source=0.d0)
    allocate(mu(MDIV), source=0.d0)
    allocate(sin_theta(MDIV), source=0.d0)

    allocate(pressure(SDIV,MDIV), source=0.d0)
    allocate(enthalpy(SDIV,MDIV), source=0.d0)
    allocate(velocity_sq(SDIV,MDIV), source=0.d0)
    allocate(energy(SDIV,MDIV), source=0.d0)
    allocate(omg(SDIV,MDIV), source=0.d0)
    allocate(F_j(SDIV,MDIV), source=0.d0)
    allocate(v_plus(SDIV), source=0.d0)
    allocate(v_minus(SDIV), source=0.d0)
    allocate(V_rr_p(SDIV), source=0.d0)
    allocate(V_rr_m(SDIV), source=0.d0)
    allocate(sound_speed(SDIV), source=0.d0)

    allocate(gama(SDIV,MDIV), source=0.d0)
    allocate(rho(SDIV,MDIV), source=0.d0)
    allocate(ww(SDIV,MDIV), source=0.d0)
    allocate(alpha(SDIV,MDIV), source=0.d0)
    allocate(sphi(SDIV,MDIV), source=0.d0)

    allocate(f_rho(SDIV, LMAX+1, SDIV), source=0.d0)
    allocate(f_gama(SDIV, LMAX+1, SDIV), source=0.d0)
    allocate(P_2n(MDIV, LMAX+1), source=0.d0)
    allocate(P1_2n_1(MDIV, LMAX+1), source=0.d0)
    allocate(sin_2n_1_theta(MDIV, LMAX), source=0.d0)
  end subroutine allocate_fields

  subroutine deallocate_fields()
    if (allocated(s_gp))             deallocate(s_gp)
    if (allocated(mu))               deallocate(mu)
    if (allocated(sin_theta))        deallocate(sin_theta)
    if (allocated(pressure))         deallocate(pressure)
    if (allocated(enthalpy))         deallocate(enthalpy)
    if (allocated(velocity_sq))      deallocate(velocity_sq)
    if (allocated(energy))           deallocate(energy)
    if (allocated(omg))              deallocate(omg)
    if (allocated(F_j))              deallocate(F_j)
    if (allocated(v_plus))           deallocate(v_plus)
    if (allocated(v_minus))          deallocate(v_minus)
    if (allocated(V_rr_p))           deallocate(V_rr_p)
    if (allocated(V_rr_m))           deallocate(V_rr_m)
    if (allocated(gama))             deallocate(gama)
    if (allocated(rho))              deallocate(rho)
    if (allocated(ww))               deallocate(ww)
    if (allocated(alpha))            deallocate(alpha)
    if (allocated(sphi))             deallocate(sphi)
    if (allocated(f_rho))            deallocate(f_rho)
    if (allocated(f_gama))           deallocate(f_gama)
    if (allocated(P_2n))             deallocate(P_2n)
    if (allocated(P1_2n_1))          deallocate(P1_2n_1)
    if (allocated(sin_2n_1_theta))   deallocate(sin_2n_1_theta)
  end subroutine deallocate_fields

end module para_mod
