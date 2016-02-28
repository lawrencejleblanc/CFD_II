!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Name: Lawrence LeBlanc							    !
!	Date: Feb 4, 2016							    !
!	Module name: Solvers  							    !
!	Description: This module contains subroutines for all of the mathematical   !
!	             solvers needed for various assignments. It is essentially a    !
!		     library that will be called upon frequently.                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module Solvers
	contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Name: LUdecomp								    !
!	Inputs: Vector size, n, and 2D coefficient matrix array			    !
!	Outputs: Coefficient matrix A						    !
!	Description: Subroutine takes in the coefficient matrix that make up the    !
!		     linear equations and performs LU decomposition.		    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE LUdecomp(n, A)
	
	IMPLICIT NONE

	! Variable declaration
	INTEGER :: k, i, j
	REAL*8 :: lam = 0
	INTEGER, INTENT(IN) :: n
	REAL*8, dimension(0:n-1,0:n-1), INTENT(INOUT) :: A

	DO k = 0, n-2, 1
	
		DO j = k + 1, n-1, 1

			IF (A(j,k) .ne. 0) THEN
				lam = A(j,k) / A(k,k)
				DO i = k + 1, n-1, 1
					A(j,i) = A(j,i) - lam * A(k,i)
					A(j,k) = lam
				ENDDO 
			ENDIF
		ENDDO
	ENDDO

	RETURN
	END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Name: LUsolve                                                               !
!       Inputs: Vector size, 2D coefficient matrix array, and source   vector       !
!       Outputs: Solution vector                                                    !
!       Description: Subroutine takes in the coefficient matrix that make up the    !
!                    linear equations and performs LU decomposition.                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE LUsolve(n, A, b)
	
	IMPLICIT NONE

	! Variable declaration
	INTEGER :: k, i
	INTEGER, INTENT(IN) :: n
	REAL*8 :: dotproduct
	REAL*8, dimension(0:n-1,0:n-1), INTENT(IN) :: A
	REAL*8, INTENT(INOUT) :: b(0:n-1)

	DO k = 1, n-1, 1

	! Calculate dot product
		dotproduct = 0.0
		DO i =  0, k-1, 1
			dotproduct = dotproduct + A(k,i)*b(i)
		ENDDO
		b(k) = b(k) - dotproduct
	ENDDO

	b(n-1) = b(n-1) / A(n-1,n-1)

	DO k = n-2, 0, -1
		dotproduct = 0.0
		DO i = k+1, n-1, 1
			dotproduct = dotproduct + A(k,i) * b(i)
		ENDDO
		b(k) = (b(k) - dotproduct) / A(k,k)
	ENDDO

	RETURN
	END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Name: gauss_seidel                                                          !
!       Inputs: x and y direction mesh size, temperature vector, matrix coefficients! 
!	        for present, west, south, north, and east points, heat generation   ! 
!               divided by thermal conductivity (source), values for each boundary  !
!		, temperature matrix to store temperatures for Gauss-Seidel, total  !
!		 number of mesh points, and boundary condition types                !
!	Outputs: Temperature values for the mesh points				    !
!	Description: This subroutine iteratively solves a system of linear          !
!		     equations using the Gauss-Seidel method.			    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gauss_seidel(xsize, ysize, An, As, Aw, Ae, Ap, source, Wbc, Ebc, Nbc, Sbc, T, area, wtype, etype, stype, ntype)

	IMPLICIT NONE

	! Variable declaration
	INTEGER, INTENT(IN) :: xsize, ysize               ! x and y size of coefficient matrix
	INTEGER :: i, ii, j, n = 0                        ! loop  counters
	INTEGER :: converge = 0
	INTEGER, INTENT(IN) :: wtype, stype, ntype, etype ! B.C. type determination
	INTEGER, INTENT(IN) :: area                       ! total number of mesh points
	REAL*8 :: T_prev(1:area)                          ! previous coefficient matrix for                                                                     ! convergence  check
	REAL*8, INTENT(INOUT) :: T(1:xsize*ysize)         ! coefficient matrix
	REAL*8, INTENT(IN) :: An, Ae, Aw, As, Ap          ! internal matrix points coefficients
	REAL*8, INTENT(IN) :: source                      ! source term  consisting of heat generation
                                                          ! over thermal conductivity
	REAL*8 :: criteria = 0.00001
	REAL*8, INTENT(IN) :: Wbc, Ebc, Nbc, Sbc          ! Boundary conditions
        REAL*8 :: Tn, Ts, Te, Tw                          ! Neighboring temperatures 

	! Initialize temperature previous vectors to 0
	DO i = 1, area
		T_prev(i) = 0
	ENDDO

	!DO i = 1, area, 1
	! Loop until convergence  criteria met or specified number of iterations
	DO WHILE (converge .eq. 0 .AND. n .lt. 10000)
		converge = 1  !Assume converged
		
		! Loop over west boundary and assign B.C. values
		
		DO i = 1, (area-xsize), xsize
			IF (wtype .eq. 1) THEN
				T(i) = Wbc
			ELSE IF (wtype .eq. 0) THEN
				T(i) = T(i+1)
			ENDIF
		ENDDO
		! Loop over east boundary and set B.C. values
		!IF
		DO i = xsize, area, xsize
			IF (etype .eq. 1) THEN
				T(i) = Ebc
			ELSE IF (etype .eq. 0) THEN
				T(i) = T(i-1)
			ENDIF
		ENDDO
		! Loop over south boundary and assign B.C. values
		DO i = 1, xsize
			IF (stype .eq. 1) THEN
				T(i) = Sbc
			ELSE IF (stype .eq. 0) THEN
				T(i) = T(i+xsize)
			END IF
		ENDDO
		! Loop over north boundary and assign B.C. values
		DO i = area, area - (xsize-1), -1
			IF (ntype .eq. 1) THEN
				T(i) = Nbc
			ELSE IF (ntype .eq. 0) THEN
				T(i) = T(i-xsize)
			END IF
		ENDDO
		! Loop over internal points
		DO i = xsize + 2, area - xsize, xsize
			j = 0
			DO WHILE (j .le. (xsize - 3))
				ii = i + j
				Tn = T(ii+xsize)
				Ts = T(ii-xsize)
				Te = T(ii+1)
				Tw = T(ii-1)
				T(ii) = 1/Ap*(-Ae*Te - As*Ts - Aw*Tw - An*Tn - source)
				If (ABS(T(ii) - T_prev(ii)) .gt. criteria) THEN
					converge = 0  
				ENDIF
				! Set T previous to newly found T value
				T_prev(ii) = T(ii)
				j = j+1
				
			ENDDO
		ENDDO
		n = n + 1
 
	ENDDO
	!ENDDO
	RETURN

	END SUBROUTINE

END Module
