program gamma_tab

    implicit none
    character(16) :: fname = "Gam2_K100.dat"
    real(8), parameter :: Gamma = 2.d0
    real(8) :: prs, rho, ee, hh
    real(8) :: step, log_rho, QB
    integer :: i
    real(8), parameter :: MB = 1.6749286d-24
    real(8), parameter :: rho_uni = 1.61930347d-18
    real(8), parameter :: prs_uni = 1.80171810d-39
    real(8), parameter :: clite = 2.99792458d10
    real(8), parameter :: K = 1.d2 * rho_uni**Gamma / prs_uni / clite**2

    log_rho = 9.d0
    step = (log10(1.d16)-log_rho)/10001

    ! write(*,"(A15,es27.18)") "K in geo. units:", K / rho_uni**Gamma * prs_uni * clite**2

    open( 221, file = trim(adjustl(fname)) )
    write(221,"(i8)") 10001
    do while(log_rho < log10(1.d16))
        rho = 1.d1**(log_rho)
        QB  = rho * rho_uni
        prs = K * QB**Gamma / prs_uni
        ee  = QB**Gamma * K / (Gamma-1.d0) + QB 
        hh  = (ee/rho_uni + prs/clite**2) / rho
        write(221,"(4es27.19)") ee/rho_uni, prs, log(hh), rho/MB
        log_rho = log_rho + step
    enddo
    close(221)
end program