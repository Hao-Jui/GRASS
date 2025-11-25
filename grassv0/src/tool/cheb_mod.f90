module cheb_mod
  implicit none
  public
  
  real(8), parameter, private :: pi = acos(-1.d0)
contains

  subroutine cheb_diff_matrix(N, D)
  ! Trefethen’s
    integer, intent(in) :: N
    real(8), intent(out) :: D(0:N,0:N)
    real(8) :: x(0:N), c(0:N), factor
    integer :: i, j

    do i = 0, N
      x(i) = cos( pi * dble(i) / dble(N) )
      if (i == 0 .or. i == N) then
        c(i) = 2.d0
      else
        c(i) = 1.d0
      end if
    end do

    do i = 0, N
      do j = 0, N
        if (i == j) then
          if (i == 0) then
            D(i,j) = ( 2.d0 * N**2 + 1.d0 ) / 6.d0
          else if (i == N) then
            D(i,j) = -( 2.d0 * N**2 + 1.d0 ) / 6.d0
          else
            D(i,j) = -x(i) / ( 2.d0 * ( 1.d0 - x(i)**2 ) )
          end if
        else
          factor = (c(i)/c(j)) * (-1.d0)**(i+j)
          D(i,j) = factor / (x(i) - x(j))
        end if
      end do
    end do
  end subroutine cheb_diff_matrix
  
end module cheb_mod
