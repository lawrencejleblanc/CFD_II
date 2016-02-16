!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Name: Lawrence LeBlanc
!Date: Feb 4, 2016
!Program name: HW1_Solution.f90
!Program description: Calls the solver submodules to solve a matrix using the LU decomposition method.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program HW1
	!Use Solvers
	IMPLICIT NONE
	REAL*8, dimension(:,:), allocatable :: A
	REAL*8, allocatable :: source(:), x(:)
	INTEGER :: ssize = 3, i, j, OpenStatus
	REAL*8 :: xsize, ysize, dx, dy, q_gen

	!! Allocate array sizes -- note that I start the indicies at 0 and not 1 -- default is 1
	!allocate(A(0:xsize-1,0:ysize-1))
	!allocate(source(0:ysize-1))
	!allocate(x(0:ysize-1))
	
	! Opens input file
	OPEN (UNIT = 1, FILE = "input.dat", STATUS = "OLD", ACTION = "READWRITE", IOSTAT = OpenStatus)

	IF (OpenStatus .gt. 0) STOP "** Cannot open input file **"
	
	! Reads values from input file
	READ (1, *) xsize, ysize, dx, dy, q_gen
	!7 FORMAT(1X, F5.2, /, 1X, F5.2, /, 1X, F5.2, /, 1X, F5.2, / ,1X, F5.2)
	WRITE (*,8) "xsize =", xsize, "ysize =", ysize, "dx =", dx, "dy =", dy, "q_gen =", q_gen 
	8 FORMAT(1X, A7, F6.2, /, 1X, A7, F6.2, /, 1X, A4, F6.2, /, 1X, A4, F6.2, /, 1X, A5, F6.2)  
	! Set the coefficient matrix and the source vector
	!A = reshape((/-2, 1, 0, 1, -2, 1, 0, 1, -2/), shape(A))

	! Call coefficient matrix constructer subroutine
	!CALL coefficient_construct

	! Call solver suboutines
	!IF (solver .eq. 1)
	!	CALL gauss_seidel(xsize, ysize, An, As, Aw, Ae, Ap, source, Wbc, Ebc, Nbc, Sbc, T)
	!ELSEIF (solver .ne. 1)
	!	CALL coefficient_construct(A, source)
	!	CALL LUdecomp(ssize, A)
	!	CALL LUsolve(ssize, A, b)
	!	        DO j =0, ssize-1
	!			x(j) = source(j)
	!			print *, x(j)
	!	 ENDDO
	!ENDIF

	CLOSE (UNIT =1)

END PROGRAM HW1
