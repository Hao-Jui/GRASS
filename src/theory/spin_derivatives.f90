module spin_derivatives
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use para_mod, only: r_ratio
  implicit none
  private
  public :: deriv_s_vec, deriv_m_vec, deriv_sm_vec, deriv_s_sub, deriv_m_sub

  interface
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      character(len=1), intent(in) :: transa, transb
      integer, intent(in) :: m, n, k, lda, ldb, ldc
      double precision, intent(in) :: alpha, beta
      double precision, intent(in) :: a(lda,*), b(ldb,*)
      double precision, intent(inout) :: c(ldc,*)
    end subroutine dgemm
  end interface

contains

  pure function deriv_s_vec(f) result(df_ds)
    use para_mod, only : SDIV, MDIV
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_ds
    call deriv_s_sub(f, df_ds)
  end function deriv_s_vec

  function deriv_m_vec(f) result(df_dm)
    use para_mod, only : SDIV, MDIV
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_dm
    call deriv_m_sub(f, df_dm)
  end function deriv_m_vec

  pure subroutine deriv_s_sub(f, df_ds)
    use para_mod, only : SDIV, MDIV, DS
    real(wp), dimension(SDIV,MDIV), intent(in)  :: f
    real(wp), dimension(SDIV,MDIV), intent(out) :: df_ds
    integer :: m
    real(wp) :: inv2ds, inv12ds

    inv2ds = 1.0_wp / (2.0_wp * DS)
    inv12ds = 1.0_wp / (12.0_wp * DS)

    if (SDIV < 5) then
      do m = 1, MDIV
        df_ds(1,m) = (f(2,m) - f(1,m)) / DS
        if (SDIV > 2) df_ds(2:SDIV-1,m) = (f(3:SDIV,m) - f(1:SDIV-2,m)) * inv2ds
        df_ds(SDIV,m) = (f(SDIV,m) - f(SDIV-1,m)) / DS
      end do
      return
    end if

    do m = 1, MDIV
      df_ds(1,m) = (-25.0_wp*f(1,m) + 48.0_wp*f(2,m) - 36.0_wp*f(3,m) + 16.0_wp*f(4,m) - 3.0_wp*f(5,m)) * inv12ds
      df_ds(2,m) = ( -3.0_wp*f(1,m) - 10.0_wp*f(2,m) + 18.0_wp*f(3,m) - 6.0_wp*f(4,m) + f(5,m)) * inv12ds
      df_ds(SDIV-1,m) = (3.0_wp*f(SDIV,m) + 10.0_wp*f(SDIV-1,m) - 18.0_wp*f(SDIV-2,m) + 6.0_wp*f(SDIV-3,m) - f(SDIV-4,m)) * inv12ds
      df_ds(SDIV,m) = (25.0_wp*f(SDIV,m) - 48.0_wp*f(SDIV-1,m) + 36.0_wp*f(SDIV-2,m) - 16.0_wp*f(SDIV-3,m) + 3.0_wp*f(SDIV-4,m)) * inv12ds
      df_ds(3:SDIV-2,m) = (-f(5:SDIV,m) + 8.0_wp*f(4:SDIV-1,m) - 8.0_wp*f(2:SDIV-3,m) + f(1:SDIV-4,m)) * inv12ds
    end do
  end subroutine deriv_s_sub

  subroutine deriv_m_sub(f, df_dm)
    use para_mod, only : SDIV, MDIV, DM, angular_collocation, COLLOCATION_UNI, D_mu_t
    real(wp), dimension(SDIV,MDIV), intent(in)  :: f
    real(wp), dimension(SDIV,MDIV), intent(out) :: df_dm
    if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
      df_dm = 0.e0_wp
      return
    end if
    if (angular_collocation /= COLLOCATION_UNI) then
      call dgemm('N', 'N', SDIV, MDIV, MDIV, 1.0_wp, f, SDIV, D_mu_t, MDIV, 0.0_wp, df_dm, SDIV)
      return
    end if
    if (MDIV < 5) then
      df_dm(:,1) = (f(:,2) - f(:,1)) / DM
      if (MDIV > 2) df_dm(:,2:MDIV-1) = (f(:,3:MDIV) - f(:,1:MDIV-2)) / (2.e0_wp * DM)
      df_dm(:,MDIV) = (f(:,MDIV) - f(:,MDIV-1)) / DM
      return
    end if

    df_dm(:,1) = (-25.e0_wp*f(:,1)+48.e0_wp*f(:,2)-36.e0_wp*f(:,3)+16.e0_wp*f(:,4)-3.e0_wp*f(:,5)) / (12.e0_wp*DM)
    df_dm(:,2) = ( -3.e0_wp*f(:,1)-10.e0_wp*f(:,2)+18.e0_wp*f(:,3)-6.e0_wp*f(:,4)+f(:,5)) / (12.e0_wp*DM)
    df_dm(:,MDIV-1) = (3.e0_wp*f(:,MDIV)+10.e0_wp*f(:,MDIV-1)-18.e0_wp*f(:,MDIV-2)+6.e0_wp*f(:,MDIV-3)-f(:,MDIV-4)) / (12.e0_wp*DM)
    df_dm(:,MDIV) = (25.e0_wp*f(:,MDIV)-48.e0_wp*f(:,MDIV-1)+36.e0_wp*f(:,MDIV-2)-16.e0_wp*f(:,MDIV-3)+3.e0_wp*f(:,MDIV-4)) / (12.e0_wp*DM)
    df_dm(:,3:MDIV-2) = (-f(:,5:MDIV)+8.e0_wp*f(:,4:MDIV-1)-8.e0_wp*f(:,2:MDIV-3)+f(:,1:MDIV-4)) / (12.e0_wp*DM)
  end subroutine deriv_m_sub

  function deriv_sm_vec(f) result(df_dsm)
    use para_mod, only : SDIV, MDIV
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_dsm, temp
    call deriv_s_sub(f, temp)
    call deriv_m_sub(temp, df_dsm)
  end function deriv_sm_vec

end module spin_derivatives
