!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Name: Lawrence LeBlanc
!Date: Feb 4, 2016
!Program name: HW2_solution.f90
!Program description: 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program HW2
	Use Solvers
	IMPLICIT NONE
	REAL*8, dimension(:,:), allocatable :: A
	REAL*8, allocatable :: b(:), x(:), T(:)
	INTEGER :: i, j, OpenStatus, solver = 1
	REAL*8 :: xsize, ysize, dx, dy, q_gen, Wbc, Ebc, Nbc, Sbc

	! Opens input file
	!OPEN (UNIT = 1, FILE = "input.dat", STATUS = "OLD", ACTION = "READWRITE", IOSTAT = OpenStatus)

	!IF (OpenStatus .gt. 0) STOP "** Cannot open input file **"
	!ENDIF

	! Reads values from input file
	!READ (1, *) xsize, ysize, dx, dy, source, solver
	!7 FORMAT(1X, F5.2, /, 1X, F5.2, /, 1X, F5.2, /, 1X, F5.2, / ,1X, F5.2, /, 1X, I1)
	
	! Print the input file values to the screen
	!WRITE (*,8) "xsize =", xsize, "ysize =", ysize, "dx =", dx, "dy =", dy, "source =", source 
	!8 FORMAT(1X, A7, F6.2, /, 1X, A7, F6.2, /, 1X, A4, F6.2, /, 1X, A4, F6.2, /, 1X, A5, F6.2)  
	
	! Allocate array sizes
	allocate(T(1:xsize*ysize))
	!allocate(b(0:ysize-1))
	!allocate(A(0:xsize-1,0:ysize-1))
	!allocate(x(0:ysize-1))
	
	xsize = 4
	ysize = 4
	dx = 1
	dy = 1
	q_gen = 1
	Wbc = 1
	Sbc = 0
	Ebc = 1
	Nbc = 1

	DO i = 1, xsize*ysize
		T(i) = 0
	ENDDO

	As = 1/(dy**2)
	An = As
	Aw = 1/(dx**2)
	Ae = Aw
	Ap = 2/(dx**2)-2/(dy**2)

	! Call solver suboutines
	IF (solver .eq. 1)
		CALL gauss_seidel(xsize, ysize, An, As, Aw, Ae, Ap, q_gen, Wbc, Ebc, Nbc, Sbc, T)
	!ELSEIF (solver .ne. 1)
	!	CALL coefficient_construct(A, source, T)
	!	CALL LUdecomp(ssize, A)
	!	CALL LUsolve(ssize, A, b)
	!	        DO j =0, ssize-1
	!			x(j) = source(j)
	!			print *, x(j)
	!	 ENDDO
	ENDIF

	!CLOSE (UNIT =1)

END PROGRAM HW1
