module lapack_interfaces_mod
  implicit none

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(lda, *), b(ldb, *)
      double precision, intent(inout) :: c(ldc, *)
    end subroutine dgemm

    subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, lda, incx, incy
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(lda, *), x(*)
      double precision, intent(inout) :: y(*)
    end subroutine dgemv

    subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: info
      double precision, intent(inout) :: a(lda, *), b(ldb, *)
    end subroutine dposv

    subroutine dpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
      character(len=1), intent(in) :: uplo
      integer, intent(in) :: n, kd, nrhs, ldab, ldb
      integer, intent(out) :: info
      double precision, intent(inout) :: ab(ldab, *), b(ldb, *)
    end subroutine dpbsv

    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      integer, intent(in) :: n, nrhs, lda, ldb
      integer, intent(out) :: info, ipiv(*)
      double precision, intent(inout) :: a(lda, *), b(ldb, *)
    end subroutine dgesv

    subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
      character(len=1), intent(in) :: trans
      integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
      integer, intent(out) :: info
      double precision, intent(inout) :: a(lda, *), b(ldb, *)
      double precision, intent(inout) :: work(*)
    end subroutine dgels
  end interface

end module lapack_interfaces_mod
