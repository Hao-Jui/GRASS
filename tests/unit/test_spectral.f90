program test_spectral
  use precision_mod, only: wp
  use test_utils
  use spectral_hub_mod
  implicit none

  integer, parameter :: NGL = 4, NLOB = 5
  real(wp) :: x_gl(NGL), w_gl(NGL)
  real(wp) :: x_lob(NLOB), w_lob(NLOB)
  real(wp) :: wsum, integral
  real(wp), parameter :: tol_sum  = 1.0e-14_wp
  real(wp), parameter :: tol_int  = 1.0e-13_wp
  real(wp), parameter :: tol_endp = 1.0e-14_wp
  integer :: i

  ! -----------------------------------------------------------------
  ! Test 1: Gauss-Legendre weights sum to 1 on [0,1]
  ! The module maps GL points from [-1,1] to [0,1] and halves weights.
  ! -----------------------------------------------------------------
  call gauss_legendre(NGL, x_gl, w_gl)
  wsum = sum(w_gl)
  call assert_near("GL weights sum = 1", 1.0_wp, wsum, tol_sum)

  ! -----------------------------------------------------------------
  ! Test 2: 4-point GL integrates x^7 exactly on [0,1]
  ! A 4-point rule is exact for polynomials up to degree 2*4-1 = 7.
  ! integral_0^1 x^7 dx = 1/8
  ! Since GL points are on [0,1] and weights sum to 1, we compute:
  !   sum_i w_i * f(x_i)
  ! -----------------------------------------------------------------
  integral = 0.0_wp
  do i = 1, NGL
    integral = integral + w_gl(i) * x_gl(i)**7
  end do
  call assert_near("GL 4-pt exactness x^7 = 1/8", 0.125_wp, integral, tol_int)

  ! -----------------------------------------------------------------
  ! Test 3: Gauss-Lobatto endpoints are 0 and 1 on [0,1]
  ! -----------------------------------------------------------------
  call gauss_lobatto(NLOB, x_lob, w_lob)
  call assert_near("Lobatto left endpoint = 0", 0.0_wp, x_lob(1),    tol_endp)
  call assert_near("Lobatto right endpoint = 1", 1.0_wp, x_lob(NLOB), tol_endp)

  ! -----------------------------------------------------------------
  ! Test 4: Gauss-Lobatto weights sum to 1 on [0,1]
  ! -----------------------------------------------------------------
  wsum = sum(w_lob)
  call assert_near("Lobatto weights sum = 1", 1.0_wp, wsum, tol_sum)

  call test_summary("spectral_hub_mod")
end program test_spectral
