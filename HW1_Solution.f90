!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Name: Lawrence LeBlanc
!Date: Feb 4, 2016
!Program name: HW1_Solution.f90
!Program description: Calls the solver submodules to solve a matrix using the LU decomposition method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program HW1
	Use Solvers
	IMPLICIT NONE
	REAL*8, dimension(:,:), allocatable :: A
	REAL*8, allocatable :: b(:), x(:)
	INTEGER :: ssize = 3, i, j

	!!Allocate array sizes -- note that I start the indicies at 0 and not 1 -- default is 1
	allocate(A(0:ssize-1,0:ssize-1))
	allocate(b(0:ssize-1))
	allocate(x(0:ssize-1))

	!Set the coefficient matrix and the source vector
	A = reshape((/-2, 1, 0, 1, -2, 1, 0, 1, -2/), shape(A))
	
	!print*, shape(A)

	DO j = 0, ssize-1
		print *, A(j,0), A(j,1), A(j,2)
	ENDDO
	
	b = reshape((/-100, 0, -10/), shape(b))

	!Call the LU Decomposition functions	
	CALL LUdecomp(ssize, A)
	CALL LUsolve(ssize, A, b)

	DO j =0, ssize-1
		x(j) = b(j)
		print *, x(j)
	ENDDO

END PROGRAM HW1
