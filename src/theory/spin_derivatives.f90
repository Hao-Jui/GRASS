module spin_derivatives
  use, intrinsic :: iso_fortran_env, only: wp => real64
  use para_mod, only: r_ratio
  implicit none
  private
  public :: deriv_s_vec, deriv_m_vec, deriv_sm_vec, deriv_s_sub, deriv_m_sub

contains

  pure function deriv_s_vec(f) result(df_ds)
    use para_mod, only : SDIV, MDIV
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_ds
    call deriv_s_sub(f, df_ds)
  end function deriv_s_vec

  pure function deriv_m_vec(f) result(df_dm)
    use para_mod, only : SDIV, MDIV
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_dm
    call deriv_m_sub(f, df_dm)
  end function deriv_m_vec

  pure subroutine deriv_s_sub(f, df_ds)
    use para_mod, only : SDIV, MDIV, DS
    real(wp), dimension(SDIV,MDIV), intent(in)  :: f
    real(wp), dimension(SDIV,MDIV), intent(out) :: df_ds
    real(wp) :: inv60DS
    if (SDIV < 5) then
      df_ds(1,:) = (f(2,:) - f(1,:)) / DS
      if (SDIV > 2) df_ds(2:SDIV-1,:) = (f(3:SDIV,:) - f(1:SDIV-2,:)) / (2.e0_wp * DS)
      df_ds(SDIV,:) = (f(SDIV,:) - f(SDIV-1,:)) / DS
      return
    end if

    if (SDIV < 7) then
      df_ds(1,:) = (-25.e0_wp*f(1,:)+48.e0_wp*f(2,:)-36.e0_wp*f(3,:)+16.e0_wp*f(4,:)-3.e0_wp*f(5,:)) / (12.e0_wp*DS)
      df_ds(2,:) = ( -3.e0_wp*f(1,:)-10.e0_wp*f(2,:)+18.e0_wp*f(3,:)-6.e0_wp*f(4,:)+f(5,:)) / (12.e0_wp*DS)
      df_ds(SDIV-1,:) = (3.e0_wp*f(SDIV,:)+10.e0_wp*f(SDIV-1,:)-18.e0_wp*f(SDIV-2,:)+6.e0_wp*f(SDIV-3,:)-f(SDIV-4,:)) / (12.e0_wp*DS)
      df_ds(SDIV,:) = (25.e0_wp*f(SDIV,:)-48.e0_wp*f(SDIV-1,:)+36.e0_wp*f(SDIV-2,:)-16.e0_wp*f(SDIV-3,:)+3.e0_wp*f(SDIV-4,:)) / (12.e0_wp*DS)
      df_ds(3:SDIV-2,:) = (-f(5:SDIV,:)+8.e0_wp*f(4:SDIV-1,:)-8.e0_wp*f(2:SDIV-3,:)+f(1:SDIV-4,:)) / (12.e0_wp*DS)
      return
    end if

    inv60DS = 1.e0_wp / (60.e0_wp * DS)
    df_ds(1,:) = (-147.e0_wp*f(1,:)+360.e0_wp*f(2,:)-450.e0_wp*f(3,:)+400.e0_wp*f(4,:) &
                  -225.e0_wp*f(5,:)+ 72.e0_wp*f(6,:)- 10.e0_wp*f(7,:)) * inv60DS
    df_ds(2,:) = ( -10.e0_wp*f(1,:)- 77.e0_wp*f(2,:)+150.e0_wp*f(3,:)-100.e0_wp*f(4,:) &
                  +  50.e0_wp*f(5,:)- 15.e0_wp*f(6,:)+  2.e0_wp*f(7,:)) * inv60DS
    df_ds(3,:) = (   2.e0_wp*f(1,:)- 24.e0_wp*f(2,:)- 35.e0_wp*f(3,:)+ 80.e0_wp*f(4,:) &
                  -  30.e0_wp*f(5,:)+  8.e0_wp*f(6,:)-       f(7,:)) * inv60DS
    df_ds(SDIV-2,:) = (        f(SDIV-6,:)-  8.e0_wp*f(SDIV-5,:)+ 30.e0_wp*f(SDIV-4,:)- 80.e0_wp*f(SDIV-3,:) &
                      + 35.e0_wp*f(SDIV-2,:)+ 24.e0_wp*f(SDIV-1,:)-  2.e0_wp*f(SDIV,:)) * inv60DS
    df_ds(SDIV-1,:) = (  -2.e0_wp*f(SDIV-6,:)+ 15.e0_wp*f(SDIV-5,:)- 50.e0_wp*f(SDIV-4,:)+100.e0_wp*f(SDIV-3,:) &
                      -150.e0_wp*f(SDIV-2,:)+ 77.e0_wp*f(SDIV-1,:)+ 10.e0_wp*f(SDIV,:)) * inv60DS
    df_ds(SDIV,:)   = (  10.e0_wp*f(SDIV-6,:)- 72.e0_wp*f(SDIV-5,:)+225.e0_wp*f(SDIV-4,:)-400.e0_wp*f(SDIV-3,:) &
                      +450.e0_wp*f(SDIV-2,:)-360.e0_wp*f(SDIV-1,:)+147.e0_wp*f(SDIV,:)) * inv60DS
    df_ds(4:SDIV-3,:) = (-f(1:SDIV-6,:) + 9.e0_wp*f(2:SDIV-5,:) - 45.e0_wp*f(3:SDIV-4,:) &
                       + 45.e0_wp*f(5:SDIV-2,:) - 9.e0_wp*f(6:SDIV-1,:) + f(7:SDIV,:)) * inv60DS
  end subroutine deriv_s_sub

  pure subroutine deriv_m_sub(f, df_dm)
    use para_mod, only : SDIV, MDIV, DM, angular_collocation, COLLOCATION_UNI, D_mu
    real(wp), dimension(SDIV,MDIV), intent(in)  :: f
    real(wp), dimension(SDIV,MDIV), intent(out) :: df_dm
    real(wp) :: inv60DM
    if (abs(r_ratio - 1.e0_wp) < epsilon(r_ratio)) then
      df_dm = 0.e0_wp
      return
    end if
    if (angular_collocation /= COLLOCATION_UNI) then
      df_dm = matmul(f, transpose(D_mu))
      return
    end if
    if (MDIV < 5) then
      df_dm(:,1) = (f(:,2) - f(:,1)) / DM
      if (MDIV > 2) df_dm(:,2:MDIV-1) = (f(:,3:MDIV) - f(:,1:MDIV-2)) / (2.e0_wp * DM)
      df_dm(:,MDIV) = (f(:,MDIV) - f(:,MDIV-1)) / DM
      return
    end if

    if (MDIV < 7) then
      df_dm(:,1) = (-25.e0_wp*f(:,1)+48.e0_wp*f(:,2)-36.e0_wp*f(:,3)+16.e0_wp*f(:,4)-3.e0_wp*f(:,5)) / (12.e0_wp*DM)
      df_dm(:,2) = ( -3.e0_wp*f(:,1)-10.e0_wp*f(:,2)+18.e0_wp*f(:,3)-6.e0_wp*f(:,4)+f(:,5)) / (12.e0_wp*DM)
      df_dm(:,MDIV-1) = (3.e0_wp*f(:,MDIV)+10.e0_wp*f(:,MDIV-1)-18.e0_wp*f(:,MDIV-2)+6.e0_wp*f(:,MDIV-3)-f(:,MDIV-4)) / (12.e0_wp*DM)
      df_dm(:,MDIV) = (25.e0_wp*f(:,MDIV)-48.e0_wp*f(:,MDIV-1)+36.e0_wp*f(:,MDIV-2)-16.e0_wp*f(:,MDIV-3)+3.e0_wp*f(:,MDIV-4)) / (12.e0_wp*DM)
      df_dm(:,3:MDIV-2) = (-f(:,5:MDIV)+8.e0_wp*f(:,4:MDIV-1)-8.e0_wp*f(:,2:MDIV-3)+f(:,1:MDIV-4)) / (12.e0_wp*DM)
      return
    end if

    inv60DM = 1.e0_wp / (60.e0_wp * DM)
    df_dm(:,1) = (-147.e0_wp*f(:,1)+360.e0_wp*f(:,2)-450.e0_wp*f(:,3)+400.e0_wp*f(:,4) &
                  -225.e0_wp*f(:,5)+ 72.e0_wp*f(:,6)- 10.e0_wp*f(:,7)) * inv60DM
    df_dm(:,2) = ( -10.e0_wp*f(:,1)- 77.e0_wp*f(:,2)+150.e0_wp*f(:,3)-100.e0_wp*f(:,4) &
                  +  50.e0_wp*f(:,5)- 15.e0_wp*f(:,6)+  2.e0_wp*f(:,7)) * inv60DM
    df_dm(:,3) = (   2.e0_wp*f(:,1)- 24.e0_wp*f(:,2)- 35.e0_wp*f(:,3)+ 80.e0_wp*f(:,4) &
                  -  30.e0_wp*f(:,5)+  8.e0_wp*f(:,6)-       f(:,7)) * inv60DM
    df_dm(:,MDIV-2) = (        f(:,MDIV-6)-  8.e0_wp*f(:,MDIV-5)+ 30.e0_wp*f(:,MDIV-4)- 80.e0_wp*f(:,MDIV-3) &
                      + 35.e0_wp*f(:,MDIV-2)+ 24.e0_wp*f(:,MDIV-1)-  2.e0_wp*f(:,MDIV)) * inv60DM
    df_dm(:,MDIV-1) = (  -2.e0_wp*f(:,MDIV-6)+ 15.e0_wp*f(:,MDIV-5)- 50.e0_wp*f(:,MDIV-4)+100.e0_wp*f(:,MDIV-3) &
                      -150.e0_wp*f(:,MDIV-2)+ 77.e0_wp*f(:,MDIV-1)+ 10.e0_wp*f(:,MDIV)) * inv60DM
    df_dm(:,MDIV)   = (  10.e0_wp*f(:,MDIV-6)- 72.e0_wp*f(:,MDIV-5)+225.e0_wp*f(:,MDIV-4)-400.e0_wp*f(:,MDIV-3) &
                      +450.e0_wp*f(:,MDIV-2)-360.e0_wp*f(:,MDIV-1)+147.e0_wp*f(:,MDIV)) * inv60DM
    df_dm(:,4:MDIV-3) = (-f(:,1:MDIV-6) + 9.e0_wp*f(:,2:MDIV-5) - 45.e0_wp*f(:,3:MDIV-4) &
                       + 45.e0_wp*f(:,5:MDIV-2) - 9.e0_wp*f(:,6:MDIV-1) + f(:,7:MDIV)) * inv60DM
  end subroutine deriv_m_sub

  function deriv_sm_vec(f) result(df_dsm)
    use para_mod, only : SDIV, MDIV
    real(wp), dimension(SDIV,MDIV), intent(in) :: f
    real(wp), dimension(SDIV,MDIV) :: df_dsm, temp
    call deriv_s_sub(f, temp)
    call deriv_m_sub(temp, df_dsm)
  end function deriv_sm_vec

end module spin_derivatives
