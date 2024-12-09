subroutine d01gaf(X,Y,N,ANS,ER,IFAIL)
!    THIS SUBROUTINE INTEGRATES A FUNCTION (Y) SPECIFIED NUMERICALLY AT N POINTS (X), WHERE N IS AT LEAST 4.
!     THE POINTS NEED NOT BE EQUALLY SPACED, BUT SHOULD BE DISTINCT AND IN ASCENDING OR DESCENDING ORDER.
!     AN ERROR ESTIMATE IS RETURNED. THE METHOD IS DUE TO GILL AND MILLER.
!
!     NAG COPYRIGHT 1975
!     MARK 5 RELEASE
!     MARK 7 REVISED IER-154 (DEC 1978)
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
	implicit none
	integer,intent(in) :: N
	real(8),intent(in) ::  X(N), Y(N)
	real(8),intent(out) :: ANS
	real(8) :: ER
	integer :: IFAIL, ii, NN
	real(8) ::  C, D1, D2, D3, H1, H2, H3, H4, R1, R2, R3, R4, S


	if( size(X) .ne. size(Y) ) stop "dimension not match, d01gaf"
	ANS = 0.d0
  	ER  = 0.d0

  	if (N.ge.4) then
		! CHECK POINTS are Strictly INCREASING or DECREASING
		H2 = x(2) - x(1)
		do ii = 3, N
			H3 = X(ii) - X(ii-1)
			if ( H2*H3 < 0.d0) then
				write(*,*) x
				write(*,*) ii,x(ii),x(ii-1)
				stop "Integration grid should be monotonic in D01GAF"
			elseif ( H2*H3 == 0.d0 ) then
				write(*,*) ii,x(ii),x(ii-1)
				stop "Integration grid should not be zero in D01GAF"
			endif
		enddo
		!  INTEGRATE OVER INITIAL INTERVAL
		D3 = (Y(2)-Y(1)) / H2
		H3 = X(3) - X(2)
		D1 = (Y(3)-Y(2)) / H3
		H1 = H2 + H3
		D2 = (D1-D3) / H1
		H4 = X(4) - X(3)
		R1 = (Y(4)-Y(3)) / H4
		R2 = (R1-D1) / (H4+H3)
		H1 = H1 + H4
		R3 = (R2-D2) / H1
		ANS = H2 * ( Y(1) + H2 * (  D3/2.d0 - H2*(D2/6.d0-(H2+2.d0*H3)*R3/12.d0)  ) )
		S = -(H2**3) * (H2*(3.d0*H2+5.d0*H4)+10.d0*H3*H1)/60.d0
		R4 = 0.d0
		!  INTEGRATE OVER CENTRAL PORTION OF RANGE
		NN = N - 1
		do ii = 3, NN
			ANS = ANS + H3*((Y(ii)+Y(ii-1))/2.d0-H3*H3*(D2+R2+(H2-H4)*R3)/12.d0)
			C = H3**3*(2.d0*H3*H3+5.d0*(H3*(H4+H2)+2.d0*H4*H2))/120.d0
			ER = ER + (C+S)*R4
			if (ii.ne.3) S = C
			if (ii.eq.3) S = S + 2.d0*C
			if ((ii-N+1).ne.0) then
				H1 = H2
				H2 = H3
				H3 = H4
				D1 = R1
				D2 = R2
				D3 = R3
				H4 = X(ii+2) - X(ii+1)
				R1 = (Y(ii+2)-Y(ii+1))/H4
				R4 = H4 + H3
				R2 = (R1-D1)/R4
				R4 = R4 + H2
				R3 = (R2-D2)/R4
				R4 = R4 + H1
				R4 = (R3-D3)/R4
			else
				! INTEGRATE OVER FINAL INTERVAL
				ANS = ANS + H4*(Y(N)-H4*(R1/2.d0+H4*(R2/6.d0 &
					+(2.d0*H3+H4)*R3/12.d0)))
				ER = ER - H4**3*R4*(H4*(3.d0*H4+5.d0*H2) &
					+10.d0*H3*(H2+H3+H4))/60.d0 + S*R4
				ANS = ANS + ER
				IFAIL = 0
			endif
		enddo
	else
		stop "N < 4 for D01GAF"
	endif

end subroutine d01gaf
