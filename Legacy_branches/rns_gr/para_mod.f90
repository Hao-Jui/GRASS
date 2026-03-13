module para_mod

  ! goals
  character(128), parameter :: eos_file = "H4"
  real(8), parameter :: M_goal = 1.35d0
  real(8), parameter :: Mb_goal= 3.2d0
  real(8), parameter :: J_goal = 1.45d0
  real(8), parameter :: chi_goal = 0.45d0
  real(8), parameter :: omc_goal = 30.d0 ! Hz

  ! for const-j law
  real(8), parameter :: A_diff = 1.d1
  ! for Uryu law
  real(8), parameter :: lambda1 = 1.1d0
  real(8), parameter :: lambda2 = 0.5d0
  integer, parameter :: uyru_p = 1
  integer, parameter :: uyru_q = 3

  ! for poly law: parA * (Omega_c - Omg)^parB
  real(8), parameter :: parA = .01d0 ! the smaller, the more differential
  real(8), parameter :: parB = .2d0 ! the greater, the more differential

  ! for counter law: cofA * omg^cofp - cofB*tanh(Omg/Omega_c^cofq)
  real(8), parameter :: cofA = 30.d0
  real(8) :: cofB
  real(8), parameter :: cofp = 0.8d0 
  real(8), parameter :: cofq = 0.9d0

  ! EoS
  integer num_tab
  real(8), allocatable, dimension(:) :: log_p, log_e, log_h, log_n0
  real(8) :: p_center, h_center
  integer :: p_at_PT
  real(8) :: M2, M4, S3

  ! grid
  integer, parameter :: LMAX = 10
  integer, parameter :: res  = 100 ! 200 should be minimal for scientific run
  real(8), parameter :: SMAX = 1.d0 - 1.d-9

  integer, parameter :: SDIV = 2*res+1
  integer, parameter :: MDIV = 2*res+1
  integer, parameter :: RDIV = 1800 ! how many points to save for static NS
  real(8) :: s_gp(SDIV), mu(MDIV)
  real(8) :: sin_theta(MDIV)
  
  logical :: output=.false.

  ! fluid
  real(8), dimension(SDIV,MDIV) :: pressure, enthalpy, velocity_sq, energy, omg, F_j
  real(8), dimension(SDIV) :: v_plus, v_minus, V_rr_p, V_rr_m

  ! metric
  real(8), dimension(SDIV,MDIV) :: gama, rho, ww, alpha

  ! bulk properties
  real(8) :: Omega_c, Omega_e, r_e, ang_mom, Omega_K, r_circ
  real(8) :: mass, mass_0, chi, T_kin, mass_p
  real(8) :: Fmax_h, F_equator_h

  ! special functions
  real(8), dimension(SDIV,LMAX+1,SDIV) :: f_rho, f_gama
  real(8), dimension(MDIV, LMAX+1) :: P_2n, P1_2n_1
  real(8) :: sin_2n_1_theta(MDIV, LMAX)

  ! parameters to be specified in the main file
  real(8) :: e_center, enthalpy_min

  real(8) :: r_ratio

  real(8), parameter :: DM = (1.d0/(dble(MDIV)-1.d0))
  real(8), parameter :: DS = (SMAX/(dble(SDIV)-1.d0))
  real(8), parameter :: s_e = 5.d-1
  real(8), parameter :: cf  = .9d0
  real(8), parameter :: C = 2.99792458d10 ! cgs
  real(8), parameter :: G = 6.67408d-8    ! cgs
  real(8), parameter :: MSUN = 1.98847d33
  real(8), parameter :: MB = 1.6749286d-24  ! baryon mass
  real(8), parameter :: pi = acos(-1.0)
  real(8), parameter :: accuracy = 1.d-5
  real(8), parameter :: tov_rmin = 1.d-15
  real(8), parameter :: KAPPA = 1.d-15 * C**2 / G
  real(8), parameter :: KSCALE = KAPPA * G / C**4
  real(8), parameter :: e_surface = 7.8 * C**2 * KSCALE
  real(8), parameter :: p_surface = 1.01 * 1.d8 * KSCALE
  real(8), parameter :: rho_uni = 1.61930347d-18

end module para_mod
