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

	SUBROUTINE LUsolve(n, A, source)
	
	IMPLICIT NONE

	! Variable declaration
	INTEGER :: k, i
	INTEGER, INTENT(IN) :: n
	REAL*8 :: dotproduct
	REAL*8, dimension(0:n-1,0:n-1), INTENT(IN) :: A
	REAL*8, INTENT(INOUT) :: source(0:n-1)

	DO k = 1, n-1, 1

	! Calculate dot product
		dotproduct = 0.0
		DO i =  0, k-1, 1
			dotproduct = dotproduct + A(k,i)*source(i)
		ENDDO
		source(k) = source(k) - dotproduct
	ENDDO

	source(n-1) = source(n-1) / A(n-1,n-1)

	DO k = n-2, 0, -1
		dotproduct = 0.0
		DO i = k+1, n-1, 1
			dotproduct = dotproduct + A(k,i) * b(i)
		ENDDO
		source(k) = (source(k) - dotproduct) / A(k,k)
	ENDDO

	RETURN
	END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Name: gauss_seidel                                                          !
!       Inputs: x and y direction mesh size, temperature vector, 
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE gauss_seidel(xsize, ysize, An, As, Aw, Ae, Ap, source, Wbc, Ebc, Nbc, Sbc, T)

	IMPLICIT NONE

	! Variable declaration
	INTEGER, INTENT(IN) :: xsize, ysize          ! x and y size of coefficient matrix
	INTEGER :: i, ii, j, n = 0                   ! loop  counters
	INTEGER :: converge = 0
	INTEGER :: area = xsize * ysize              ! total number of mesh points
	REAL*8 :: T_prev(1:area)                     ! temp solution for convergence check
	REAL*8 :: INTENT(INOUT) :: T(1:xsize*ysize)  ! coefficient matrix
	INTEGER, INTENT(IN) :: An, Ae, Aw, As        ! internal matrix points coefficients
	REAL*8, INTENT(IN) :: source(1:ysize)        ! source vector consisting of heat generation
                                                     ! over thermal conductivity
	REAL*8 :: criteria = 0.0001
	REAL*8, INTENT(IN) :: Wbc, Ebc, Nbc, Sbc     ! Boundary conditions
        REAL*8, INTENT(IN) :: Tn, Ts, Te, Tw         ! Neighboring temperatures 

	! Initialize temperature previous vectors to 0
	DO i = 1, area
		T_prev(i) = 0
	ENDDO

	! Loop until convergence  criteria met or specified number of iterations
	DO WHILE (converge .ne. 1 .or. n .lt. 100)
		! Loop over west boundary and assign B.C. values
		DO i = 0, (area-xsize), xsize
			T(i) = Wbc
		ENDDO
		! Loop over south boundary and assign B.C. values
		DO i = 2, xsize
			T(i) = Sbc
		ENDDO 
		! Loop over east boundary and set B.C. values
		DO i = xsize, area, xsize
			T(i) = Ebc
		ENDDO
		! Loop over north boundary and assign B.C. values
		DO i = area, area - xsize, -1
			T(i) = Nbc
		ENDDO
		! Loop over internal points
		DO i = xsize + 2, area - xsize, xsize
			j = 0
			DO WHILE (j .lt. xsize - 3)
				ii = i + j
				Tn = T(ii+xsize)
				Ts = T(ii-xsize)
				Te = T(ii+1)
				Tw = T(ii-1)
				T(ii) = 1/Ap*(-Ae*Te - As*Ts - Aw*Tw - An*Tn - source(ii))
				if (ABS(T(ii) - T_prev(ii)) .gt. criteria)
					  
				ENDIF
				j=j+1
				
			ENDDO
		ENDDO
		n = n + 1
	ENDDO

	END SUBROUTINE

END Module
