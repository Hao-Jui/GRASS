! spline_mod — not-a-knot cubic spline construction
!
! Given knots xp(1:n) and values yp(1:n), computes coefficients
! b, c, d such that on interval [xp(i), xp(i+1)]:
!
!   S(x) = yp(i) + b(i)*t + c(i)*t^2 + d(i)*t^3
!
! where t = x - xp(i).
!
! Not-a-knot BCs: third derivative continuous at xp(2) and
! xp(n-1), giving 4th-order accuracy on smooth data.
!
! Cost: O(n) build via Thomas algorithm on a reduced
! (n-2)x(n-2) tridiagonal system.
module spline_mod
  use precision_mod, only: wp
  implicit none
  private

  type, public :: spline_coeff
    integer :: n = 0
    real(wp), allocatable :: b(:), c(:), d(:)
  end type

  public :: build_spline, build_spline_segmented

contains

  ! ================================================================
  !  Build not-a-knot cubic spline on a single segment
  ! ================================================================
  subroutine build_spline(xp, yp, n, coeff)
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: xp(n), yp(n)
    type(spline_coeff), intent(out) :: coeff
    real(wp), allocatable :: h(:), sigma(:)
    integer :: i, nm

    coeff%n = n
    nm = n - 1

    if (n <= 1) then
      allocate(coeff%b(0), coeff%c(0), coeff%d(0))
      return
    end if

    allocate(coeff%b(nm), coeff%c(nm), coeff%d(nm))

    if (n == 2) then
      coeff%b(1) = (yp(2) - yp(1)) / (xp(2) - xp(1))
      coeff%c(1) = 0.0_wp
      coeff%d(1) = 0.0_wp
      return
    end if

    if (n == 3) then
      call build_3pt(xp, yp, coeff%b, coeff%c, coeff%d)
      return
    end if

    allocate(h(nm), sigma(n))
    do i = 1, nm
      h(i) = xp(i+1) - xp(i)
    end do

    call solve_nak(yp, n, h, sigma)

    do i = 1, nm
      coeff%c(i) = sigma(i)
      coeff%d(i) = (sigma(i+1) - sigma(i)) / (3.0_wp * h(i))
      coeff%b(i) = (yp(i+1) - yp(i)) / h(i) &
        - h(i) * (2.0_wp * sigma(i) + sigma(i+1)) / 3.0_wp
    end do
  end subroutine build_spline

  ! ================================================================
  !  Solve for sigma(1:n) = S''(xp)/2 with not-a-knot BCs.
  !
  !  Interior equations (i=2..n-1):
  !    h(i-1)*sigma(i-1) + 2(h(i-1)+h(i))*sigma(i)
  !      + h(i)*sigma(i+1) = rhs(i)
  !  where rhs(i) = 3*((yp(i+1)-yp(i))/h(i)
  !                   - (yp(i)-yp(i-1))/h(i-1))
  !
  !  Not-a-knot: d(1)=d(2) and d(n-2)=d(n-1), giving
  !    sigma(1) = ((h1+h2)*sigma(2) - h1*sigma(3)) / h2
  !    sigma(n) = ((h_{n-2}+h_{n-1})*sigma(n-1)
  !              - h_{n-1}*sigma(n-2)) / h_{n-2}
  !
  !  Substitute into rows 2 and n-1 to get a (n-2)x(n-2)
  !  tridiagonal system on sigma(2:n-1), solved by Thomas.
  ! ================================================================
  subroutine solve_nak(yp, n, h, sigma)
    integer,  intent(in)  :: n
    real(wp), intent(in)  :: yp(n), h(n-1)
    real(wp), intent(out) :: sigma(n)
    integer  :: i, m
    real(wp), allocatable :: dd(:), su(:), sb(:), r(:)
    real(wp) :: w

    m = n - 2
    allocate(dd(m), su(m), sb(m), r(m))
    sb(1) = 0.0_wp; su(m) = 0.0_wp

    ! Fill interior equations, mapped j=1..m ↔ i=2..n-1
    do i = 2, n - 1
      dd(i-1) = 2.0_wp * (h(i-1) + h(i))
      r(i-1)  = 3.0_wp * ((yp(i+1) - yp(i)) / h(i) &
        - (yp(i) - yp(i-1)) / h(i-1))
      if (i < n - 1) su(i-1) = h(i)
      if (i > 2)     sb(i-1) = h(i-1)
    end do

    ! Fold not-a-knot left BC into row j=1 (i=2):
    ! sigma(1) = ((h1+h2)*sigma(2) - h1*sigma(3)) / h2
    ! Original: h(1)*sigma(1) + ... → substitution adds:
    dd(1) = dd(1) + h(1) * (h(1) + h(2)) / h(2)
    if (m >= 2) su(1) = su(1) - h(1)**2 / h(2)

    ! Fold not-a-knot right BC into row j=m (i=n-1):
    ! sigma(n) = ((h_{n-2}+h_{n-1})*sigma(n-1)
    !            - h_{n-1}*sigma(n-2)) / h_{n-2}
    dd(m) = dd(m) &
      + h(n-1) * (h(n-2) + h(n-1)) / h(n-2)
    if (m >= 2) sb(m) = sb(m) - h(n-1)**2 / h(n-2)

    ! Thomas forward sweep
    do i = 2, m
      w = sb(i) / dd(i-1)
      dd(i) = dd(i) - w * su(i-1)
      r(i)  = r(i)  - w * r(i-1)
    end do

    ! Back substitution → sigma(2:n-1)
    sigma(n-1) = r(m) / dd(m)
    do i = m - 1, 1, -1
      sigma(i+1) = (r(i) - su(i) * sigma(i+2)) / dd(i)
    end do

    ! Recover endpoints from not-a-knot
    sigma(1) = ((h(1) + h(2)) * sigma(2) &
      - h(1) * sigma(3)) / h(2)
    sigma(n) = ((h(n-2) + h(n-1)) * sigma(n-1) &
      - h(n-1) * sigma(n-2)) / h(n-2)
  end subroutine solve_nak

  ! ================================================================
  !  3-point: quadratic (not-a-knot = linear 3rd deriv = 0)
  ! ================================================================
  subroutine build_3pt(xp, yp, b, c, d)
    real(wp), intent(in)  :: xp(3), yp(3)
    real(wp), intent(out) :: b(2), c(2), d(2)
    real(wp) :: h1, h2, dd1, dd2, dd12

    h1  = xp(2) - xp(1)
    h2  = xp(3) - xp(2)
    dd1 = (yp(2) - yp(1)) / h1
    dd2 = (yp(3) - yp(2)) / h2
    dd12 = (dd2 - dd1) / (h1 + h2)

    b(1) = dd1 - dd12 * h1; c(1) = dd12; d(1) = 0.0_wp
    b(2) = dd1 + dd12 * h1; c(2) = dd12; d(2) = 0.0_wp
  end subroutine build_3pt

  ! ================================================================
  !  Segmented spline: independent splines per smooth segment,
  !  separated at phase-transition indices.
  !
  !  pt_idx(1:n_pt) are the PT indices in the table.
  !  The interval (pt_idx(k)-1, pt_idx(k)) is the PT gap
  !  and gets linear interpolation (c=d=0).
  ! ================================================================
  subroutine build_spline_segmented(xp, yp, n, pt_idx, n_pt, &
      coeff)
    integer,  intent(in)  :: n, n_pt
    real(wp), intent(in)  :: xp(n), yp(n)
    integer,  intent(in)  :: pt_idx(n_pt)
    type(spline_coeff), intent(out) :: coeff
    integer :: seg_start, seg_end, ip, nm, i
    type(spline_coeff) :: sc

    coeff%n = n
    nm = n - 1
    allocate(coeff%b(nm), coeff%c(nm), coeff%d(nm))
    coeff%b = 0.0_wp; coeff%c = 0.0_wp; coeff%d = 0.0_wp

    if (n_pt == 0) then
      call build_spline(xp, yp, n, sc)
      coeff%b(1:nm) = sc%b(1:nm)
      coeff%c(1:nm) = sc%c(1:nm)
      coeff%d(1:nm) = sc%d(1:nm)
      return
    end if

    seg_start = 1
    do ip = 1, n_pt
      seg_end = pt_idx(ip) - 1
      call fill_segment(xp, yp, n, seg_start, seg_end, coeff)

      ! PT gap: linear between seg_end and pt_idx(ip)
      if (seg_end >= 1 .and. seg_end < n) then
        i = seg_end
        coeff%b(i) = (yp(i+1) - yp(i)) &
          / (xp(i+1) - xp(i))
        coeff%c(i) = 0.0_wp
        coeff%d(i) = 0.0_wp
      end if

      seg_start = pt_idx(ip)
    end do

    ! Last segment
    call fill_segment(xp, yp, n, seg_start, n, coeff)

  contains
    subroutine fill_segment(xp_a, yp_a, na, is, ie, co)
      integer,  intent(in)    :: na, is, ie
      real(wp), intent(in)    :: xp_a(na), yp_a(na)
      type(spline_coeff), intent(inout) :: co
      type(spline_coeff) :: tmp
      integer :: ns, j

      ns = ie - is + 1
      if (ns < 2) return
      call build_spline(xp_a(is:ie), yp_a(is:ie), ns, tmp)
      do j = 1, ns - 1
        co%b(is + j - 1) = tmp%b(j)
        co%c(is + j - 1) = tmp%c(j)
        co%d(is + j - 1) = tmp%d(j)
      end do
    end subroutine
  end subroutine build_spline_segmented

end module spline_mod
