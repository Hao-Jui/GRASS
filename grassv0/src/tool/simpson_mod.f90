module simpson_mod
  implicit none

contains

  !--------------------------------------
  ! Vectorized Simpson's rule for multiple function arrays
  !--------------------------------------
  function simpson_1d(fvals, a, b) result(integrals)
    implicit none
    real(8), intent(in) :: fvals(:,:)  ! shape (Npoints, Nfuncs)
    real(8), intent(in) :: a, b
    real(8) :: integrals(size(fvals,2))
    real(8) :: h
    integer :: npts
    real(8) :: w(size(fvals,1)) ! Simpson weights

    npts = size(fvals,1)
    if (mod(npts,2) == 0) then
       print *, "Simpson requires odd number of points. Increment N by 1."
       stop
    endif

    ! Step size
    h = (b - a)/dble(npts-1)

    ! Build Simpson weights
    w = 2.d0
    w(1) = 1.d0
    w(npts) = 1.d0
    w(2:npts-1:2) = 4.d0

    ! Compute integrals
    integrals = matmul(w, fvals) * (h/3.d0)
  end function simpson_1d

end module simpson_mod
