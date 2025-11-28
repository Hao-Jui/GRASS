module rotation_uniform
  use toolkit_mod
  use para_mod
  use spin_helper
  implicit none
contains

subroutine spin
  implicit none
  integer :: m, s, n_of_it
  real(8) :: r_e_old, r_e_new, r_e_new_sq
  real(8) :: gama_pole_h, gama_center_h, gama_equator_h
  real(8) :: rho_pole_h, rho_center_h, rho_equator_h, ww_equator_h
  real(8) :: sphi_pole_h, sphi_center_h, sphi_equator_h
  real(8) :: root_mphi_re
  logical :: zero_scalar_mode
  character(32) :: fil1, fil2, fil3, fil4, fil5, fil6

  ! externs used in precompute
  real(8) :: n0_at_e

  zero_scalar_mode = merge(.true., .false., active_theory == THEORY_GR)

  dif = 1.d0
  n_of_it = 0
  r_e_new = r_e
  r_e_new_sq = r_e_new**2

  if ( maxval(sphi*sqrt(B_coup)) < 1.d-3 ) sphi = sphi*1.d1
  if (zero_scalar_mode) sphi = 0.d0
  if ( any(isnan(sphi)) ) stop "NaN found in sphi"
  
  ! ---------------------------------------------------------------
  ! Iteration
  ! ---------------------------------------------------------------
  !write(*,*) r_ratio, h_center
  call allocate_workspace

  do while( dif > 1.d-7 .or. n_of_it < 2 )
    if (zero_scalar_mode) sphi = 0.d0
    sphi_m   = maxval( sphi(:,1) * sqrt(B_coup) )
    call rescale_metric(r_e_new_sq)

    ! --- Compute r_e ---
    r_e_old = r_e_new
    call update_equatorial_radius( r_e_old, &                            ! input 
         r_e_new, dif, &                                                 ! output           
         sphi_pole_h, gama_pole_h, rho_pole_h, &                         ! output
         gama_equator_h, rho_equator_h, ww_equator_h, sphi_equator_h, &  ! output
         sphi_center_h, gama_center_h, rho_center_h )                    ! output   

    r_e_new_sq = r_e_new**2

    if ( n_of_it > 50 .and. mod(n_of_it,50)==0 ) then 
      write(*,'(A,i4,A,es12.4,A,2es10.2,A,3es10.2,4es18.9)') 'iter = ', n_of_it, ', diff :', dif, &
        ' sphi :', sphi_center_h*r_e_old*sqrt(B_coup), sphi_m, &
        '  |', gama_center_h, rho_center_h, alpha(1,1), r_e_old, r_e_new
    endif

    call update_angular_velocity(r_e_new, gama_pole_h, rho_pole_h, gama_equator_h, rho_equator_h, &
         sphi_pole_h, sphi_equator_h, ww_equator_h)

    call update_eos_and_velocity(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h)

    root_mphi_re = sqrt(mphi_r * r_e_new_sq)

    call get_all_targets(r_e_new, gama_pole_h, rho_pole_h, sphi_pole_h, root_mphi_re, &
                         target_rho, target_gama, target_ww, target_sphi)
    
    call relaxation(r_e_new, target_rho, target_gama, target_ww, target_sphi, root_mphi_re, n_of_it)

    ! ---------------------------------------------------------------
    ! Divergence check & rigid rotation enforcement
    ! ---------------------------------------------------------------
    if (abs(rho(2,1))>100.d0 .or. abs(gama(2,1))>300.d0 .or. abs(ww(2,1))>100.d0 &
        .or. abs(sphi(2,1))>10.d0) then
      write(*,"(i5,4es18.9)") n_of_it, rho(2,1), gama(2,1), ww(2,1), sphi(2,1)
      stop "something diverged"
    end if

    if (r_ratio == 1.d0) call impose_rigid_rotation()

    ! ---------------------------------------------------------------
    ! Fourth equation (alpha), reuse caches where possible
    ! ---------------------------------------------------------------
    call update_alpha_potential(r_e_new, dg_s_cache, dg_m_cache, dr_s_cache, dr_m_cache, dww_s_cache, dww_m_cache, &
         ds_s_cache, ds_m_cache, d2g_ss_cache, d2g_mm_cache, e_rsm_cache)
    
    n_of_it = n_of_it + 1
    if ( n_of_it > 2000 ) stop "Probably won't converge"
  enddo
  !write(*,*) r_ratio, h_center, n_of_it
  call deallocate_workspace
  n_of_relaxation_steps = n_of_relaxation_steps + n_of_it
  ! --- End of iteration

  ! compute omega & outputs (unchanged)
  omg(:,:) = Omega_c / r_e_new
  Omega_c  = Omega_c / r_e_new
  Omega_e  = Omega_c
  r_e      = r_e_new
  if (zero_scalar_mode) then
    sphi = 0.d0
    sphi_c = 0.d0
    sphi_m = 0.d0
  else
    sphi_c   = sphi(1,1) * sqrt(B_coup)
    sphi_m   = maxval( sphi(:,1) * sqrt(B_coup) )
  end if

  if (output) call output_helper(D2_metric_rho, D2_metric_omega)

contains
  subroutine rescale_metric(factor)
    implicit none
    real(8), intent(in) :: factor
    real(8) :: inv_factor, sqrt_factor
    integer :: s, m

    inv_factor = 1.d0 / factor
    sqrt_factor = sqrt(factor)

    do s = 1, SDIV
      do m = 1, MDIV
        rho  (s,m) = rho  (s,m) * inv_factor
        gama (s,m) = gama (s,m) * inv_factor
        alpha(s,m) = alpha(s,m) * inv_factor
        ww   (s,m) = ww   (s,m) * sqrt_factor
        sphi (s,m) = sphi (s,m) / sqrt_factor
      end do
    end do
  end subroutine rescale_metric

  subroutine impose_rigid_rotation()
    real(8) :: rho_s1, sphi_s1, gama_s1
    integer :: s
    do s = 1, SDIV
      rho_s1 = rho(s,1)
      sphi_s1 = sphi(s,1)
      gama_s1 = gama(s,1)
      rho(s,:)  = rho_s1
      sphi(s,:) = sphi_s1
      gama(s,:) = gama_s1
      ww(s,:)   = 0.d0
    end do
  end subroutine impose_rigid_rotation
end subroutine spin

end module rotation_uniform
