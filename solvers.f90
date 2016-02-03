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
	
END Module
