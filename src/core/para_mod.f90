module para_mod
  use precision_mod, only: wp
  implicit none
  ! -- Theory selection ------------------------------------------------------
  integer, parameter :: THEORY_GR = 0
  integer, parameter :: THEORY_ST = 1
  integer :: active_theory = THEORY_ST

  ! hybrid / anderson
  ! -- Running option --------------------------------------------------------
  integer, parameter :: MODE_REGRID  = 1
  integer, parameter :: MODE_DEFAULT = 2

  integer :: run_mode = MODE_REGRID

  ! -- Rotation configuration ------------------------------------------------
  ! uniform / const_j / uryu
  character(len=20) :: solver_type = "uniform"

  ! -- Solver state ----------------------------------------------------------
  logical :: output = .false.
  logical :: timing = .false.
  logical :: use_shoot_1d = .true.      ! adjust hc while keeping rep constant
  character(len=20) :: FIX1 = "Mb_goal"
  character(len=20) :: FIX2 = "chi_goal"

  ! -- Resolutions -----------------------------------------------------------
  integer, parameter :: res  = 10
  integer, parameter :: s_pwr = 1
  integer :: SDIV = 60 * res + 1
  integer :: MDIV = 4 * res + 1

  ! -- Target quantities -----------------------------------------------------
  character(len=128) :: eos_file = "MPA1"
  real(wp) :: M_goal   = 1.2e0_wp
  real(wp) :: Mb_goal  = 1.8e0_wp
  real(wp) :: J_goal   = 1.6e0_wp
  real(wp) :: chi_goal = 0.1e0_wp
  real(wp) :: omc_goal = 30.e0_wp

  real(wp) :: B_goal   = 16.e2_wp
  real(wp) :: mphi_goal = 1.e0_wp

  ! -- Rotation-law parameters (KEH, Uryu enabled) --------------------------
  real(wp) :: A_diff  = 0.5e0_wp
  real(wp) :: lambda1 = 1.5e0_wp
  real(wp) :: lambda2 = 0.3e0_wp
  integer :: uyru_p  = 1
  integer :: uyru_q  = 3

  real(wp) :: parA = 1.e0_wp
  real(wp) :: parB = 1.e0_wp

  real(wp) :: cofA = 30.e0_wp
  real(wp) :: cofB = 0.e0_wp
  real(wp) :: cofp = 0.8e0_wp
  real(wp) :: cofq = 0.9e0_wp

  ! -- Equation of state -----------------------------------------------------
  logical :: phase_transition = .false.
  integer :: num_tab = 0
  integer :: p_at_PT = 0
  real(wp), allocatable :: log_p(:), log_e(:), log_h(:), log_n0(:)
  real(wp) :: p_center = 0.e0_wp
  real(wp) :: h_center = 0.e0_wp
  real(wp) :: e_center = 0.e0_wp
  real(wp) :: enthalpy_min = 0.e0_wp

  ! -- Grid configuration ----------------------------------------------------
  integer, parameter :: LMAX = 10
  real(wp) :: SMAX  = 1.e0_wp - (1.e-1_wp)**(9.e0_wp / dble(s_pwr))
  integer, parameter :: RDIV = 1800

  real(wp) :: DS = 0.e0_wp
  real(wp) :: DM = 0.e0_wp
  real(wp), parameter :: s_e = 0.5e0_wp

  real(wp), allocatable :: s_gp(:), mu(:), sin_theta(:)

  ! -- Disk helper quantities (kept for compatibility) ----------------------
  logical :: disk_present = .false.
  real(wp) :: edge_in = 800.e0_wp
  real(wp) :: s_inner = 0.5e0_wp
  real(wp) :: j_disk  = 4.5e0_wp
  real(wp) :: p_max_disk = 0.e0_wp
  real(wp) :: h_max_disk = 0.e0_wp
  integer :: i_isco_p = 1
  integer :: i_isco_m = 1

  ! Fluid
  real(wp), allocatable :: pressure(:,:), enthalpy(:,:), velocity_sq(:,:), &
                          energy(:,:), omg(:,:), F_j(:,:)
  real(wp), allocatable :: v_plus(:), v_minus(:), V_rr_p(:), V_rr_m(:), sound_speed(:)

  ! Metric
  real(wp), allocatable :: gama(:,:), rho(:,:), ww(:,:), alpha(:,:), sphi(:,:)

  ! Scalar field
  logical :: has_scalar = .false.
  real(wp) :: B_coup = 0.e0_wp
  real(wp) :: mphi_r = 0.e0_wp
  real(wp) :: sphi_c     = 0.e0_wp
  real(wp) :: sphi_m     = 0.e0_wp
  real(wp) :: r_sphi_max = 0.e0_wp  ! physical equatorial radius at max(sphi)

  real(wp) :: B_burn_init = 13.e0_wp
  real(wp) :: mphi_burn_seed = 0.05e0_wp
  integer, parameter :: scalar_burn_max_iter = 200
  real(wp), parameter :: mphi_burn_threshold = 0.05e0_wp

  ! Bulk properties
  real(wp) :: Omega_c = 0.e0_wp
  real(wp) :: Omega_e = 0.e0_wp
  real(wp) :: Omega_K = 0.e0_wp
  real(wp) :: r_e     = 1.e0_wp
  real(wp) :: r_ratio = 1.e0_wp
  real(wp) :: r_circ  = 0.e0_wp
  real(wp) :: ang_mom = 0.e0_wp
  real(wp) :: mass    = 0.e0_wp
  real(wp) :: mass_0  = 0.e0_wp
  real(wp) :: mass_p  = 0.e0_wp
  real(wp) :: chi     = 0.e0_wp
  real(wp) :: T_kin   = 0.e0_wp
  real(wp) :: I_inertia = 0.e0_wp
  real(wp) :: Love2   = 0.e0_wp
  real(wp) :: Fmax_h  = 0.e0_wp
  real(wp) :: F_equator_h = 0.e0_wp
  ! Multipole information
  real(wp) :: M2 = 0.e0_wp
  real(wp) :: M4 = 0.e0_wp
  real(wp) :: S3 = 0.e0_wp

  integer :: n_of_relaxation_steps = 0

  ! Green's functions (Legendre weights — precomputed on grid)
  real(wp), allocatable :: P_2n(:,:), P1_2n_1(:,:), sin_2n_1_theta(:,:)

  ! Timing
  real(wp) :: start = 0.e0_wp
  real(wp) :: finish = 0.e0_wp

  ! -- Physical constants ----------------------------------------------------
  real(wp), parameter :: C    = 2.99792458e10_wp
  real(wp), parameter :: G    = 6.67408e-8_wp
  real(wp), parameter :: MSUN = 1.98847e33_wp
  real(wp), parameter :: MB   = 1.6749286e-24_wp
  real(wp), parameter :: pi   = acos(-1.e0_wp)

  real(wp), parameter :: accuracy  = 1.e-7_wp
  real(wp), parameter :: tov_rmin  = 1.e-15_wp
  real(wp), parameter :: KAPPA     = 1.e-15_wp * C**2 / G
  real(wp), parameter :: KSCALE    = KAPPA * G / C**4
  real(wp), parameter :: e_surface = 7.8e0_wp * C**2 * KSCALE
  real(wp), parameter :: p_surface = 1.01e8_wp * KSCALE
  real(wp), parameter :: rho_uni   = 7.4259154861063358e-19_wp
  real(wp), parameter :: prs_uni   = rho_uni/(C*1.e5_wp)**2
  real(wp), parameter :: f_uni     = 2.029739818539300e5_wp
  real(wp), parameter :: hbar      = 6.582119569e-16_wp
  real(wp), parameter :: l_uni     = 1.4769994423016508e0_wp
  real(wp), parameter :: n_sat     = 2.7e14_wp
  real(wp), parameter :: scalarton = hbar * C / l_uni / 1.e5_wp
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
      call apply_st_defaults()
    case default
      stop "initialize_theory: unknown theory mode"
    end select

    solver_type = trim(to_lower_str(adjustl(solver_type)))

    DS   = SMAX / (dble(SDIV) - 1.e0_wp)
    DM   = 1.e0_wp  / (dble(MDIV) - 1.e0_wp)
    s_inner = edge_in**(1.e0_wp / dble(s_pwr)) / (edge_in**(1.e0_wp / dble(s_pwr)) + 1.e0_wp)

    call allocate_fields()
  end subroutine initialize_theory

  subroutine apply_gr_defaults()
    has_scalar = .false.
    B_goal  = 0.e0_wp
    mphi_goal = 0.e0_wp
    B_coup = 0.e0_wp
    mphi_r = 0.e0_wp
    sphi_c = 0.e0_wp
    sphi_m = 0.e0_wp
  end subroutine apply_gr_defaults

  subroutine apply_st_defaults()
    has_scalar = .true.
  end subroutine apply_st_defaults

  subroutine allocate_fields()
    call deallocate_fields()

    allocate(s_gp(SDIV), source=0.e0_wp)
    allocate(mu(MDIV), source=0.e0_wp)
    allocate(sin_theta(MDIV), source=0.e0_wp)

    allocate(pressure(SDIV,MDIV), source=0.e0_wp)
    allocate(enthalpy(SDIV,MDIV), source=0.e0_wp)
    allocate(velocity_sq(SDIV,MDIV), source=0.e0_wp)
    allocate(energy(SDIV,MDIV), source=0.e0_wp)
    allocate(omg(SDIV,MDIV), source=0.e0_wp)
    allocate(F_j(SDIV,MDIV), source=0.e0_wp)
    allocate(v_plus(SDIV), source=0.e0_wp)
    allocate(v_minus(SDIV), source=0.e0_wp)
    allocate(V_rr_p(SDIV), source=0.e0_wp)
    allocate(V_rr_m(SDIV), source=0.e0_wp)
    allocate(sound_speed(SDIV), source=0.e0_wp)

    allocate(gama(SDIV,MDIV), source=0.e0_wp)
    allocate(rho(SDIV,MDIV), source=0.e0_wp)
    allocate(ww(SDIV,MDIV), source=0.e0_wp)
    allocate(alpha(SDIV,MDIV), source=0.e0_wp)
    allocate(sphi(SDIV,MDIV), source=0.e0_wp)

    allocate(P_2n(MDIV, LMAX+1), source=0.e0_wp)
    allocate(P1_2n_1(MDIV, LMAX+1), source=0.e0_wp)
    allocate(sin_2n_1_theta(MDIV, LMAX), source=0.e0_wp)
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
    if (allocated(P_2n))             deallocate(P_2n)
    if (allocated(P1_2n_1))          deallocate(P1_2n_1)
    if (allocated(sin_2n_1_theta))   deallocate(sin_2n_1_theta)
  end subroutine deallocate_fields

end module para_mod
