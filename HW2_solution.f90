!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Name: Lawrence LeBlanc
!Date: Feb 4, 2016
!Program name: HW2_solution.f90
!Program description: This program solves the 2-D steady state heat equation with no generation.
!		      The program reads in mesh parameters and B.C. from an input file and solves
!	              the heat equation either iteratively with the gauss seidel method or directly
!		      with LU Decomposition. Both method's algorithms are found in the solvers file.
!		      A seperate file contains a subroutine to create the matrices needed for LU
! 		      Decomposition. The results are written to a tecplot file for analysis.  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program HW2
	Use Solvers !, utility
	IMPLICIT NONE

	! Variable declaration
	REAL*8, dimension(:,:), allocatable :: A
	REAL*8, allocatable :: b(:), x(:), T(:)
	INTEGER :: i, j, wtype, etype, stype, ntype, OpenStatus, solver = 1
	REAL*8 :: dx, dy, q_gen, Wbc, Ebc, Nbc, Sbc, As, An, Ae, Ap, Aw, y, xx
	INTEGER :: xsize, ysize, area, source_size

	! Opens input file and output tecplot file
	OPEN (UNIT = 1, FILE = "input.dat", IOSTAT = OpenStatus)
	OPEN (UNIT = 2, FILE = "output.dat", ACTION = "WRITE")
	IF (OpenStatus .gt. 0) THEN
		STOP "** Cannot open input file **"
	ENDIF

	! Reads values from input file
	! Input file must be in this order and format
	! can include comments on the same lines only. must be after
	! the required numbers/integers/floats
	READ (1, *) xsize
	READ (1, *) ysize
	READ (1, *) dx
	READ (1, *) dy
	READ (1, *) q_gen
	READ (1, *) Wbc, wtype
	READ (1, *) Nbc, ntype
	READ (1, *) Sbc, stype
	READ (1, *) Ebc, etype
	READ (1, *) solver 
	
	! Print the input file values to the screen
	WRITE (*,8) "xsize =", xsize, "ysize =", ysize, "dx =", dx, "dy =", dy 
	8 FORMAT(1X, A7, I2, /, 1X, A7, I2, /, 1X, A4, F6.2, /, 1X, A4, F6.2)
	
	area = xsize * ysize
	
	! Calculate matrix coefficients
	As = 1/(dy**2)
	An = As
	Aw = 1/(dx**2)
	Ae = Aw
	Ap = (-2/(dx**2))-(2/(dy**2))
	
	print*, As, An, Aw, Ae, Ap
	
	print*, ""

	! Call solver suboutines
	IF (solver .eq. 1)THEN
		allocate(T(1:xsize*ysize))
		! Initialize temperature array to 0. These values are sued as initial T_prev
		DO i = 1, xsize*ysize, 1
			T(i) = 0
		ENDDO
		
		CALL gauss_seidel(xsize, ysize, An, As, Aw, Ae, Ap, q_gen, Wbc, Ebc, Nbc, Sbc, T, area, wtype, etype, stype, ntype)
		
		! Write the header to the output tecplot file
		WRITE(2,*) 'VARIABLES = "X" "Y" "TEMPERATURE"'
		WRITE(2,9) "zone I=", xsize, "J=", ysize, " SOLUTIONTIME=", 0.00, ", F=POINT"
		9 FORMAT(2X, A7, I1, 2X, A3, I1, 2X, A14, F6.3, 2X, A7)
		
		! Write the x position, y position, and temperature values to the tecplot file
		DO j = 1, area, 1
			xx = (MOD(j, (xsize)) - 1)*dx
			IF (xx .lt. 0) THEN
				xx = xx + (xsize*dx)
			ENDIF
			y = ((j - 1)/(xsize))*dy
			WRITE(2,11) xx, y, T(j)
			11 FORMAT(F10.5, 3X, F10.5, 3X, F10.5)
		ENDDO
		
	!ELSE
		!ALLOCATE(A(0:xsize-1,0:ysize-1))
		!ALLOCATE(b(0:ysize-1))
		!ALLOCATE(x(0:ysize-1))
		!CALL coefficient_construct(A, b, T)
		!CALL LUdecomp(ysize, A)
		!CALL LUsolve(ysize, A, b)
		        !DO j =0, ysize-1
				!x(j) = b(j)
				!print *, x(j)
			!ENDDO
	ENDIF

	CLOSE (1) ! close input file
	CLOSE (2) ! close output file

END PROGRAM HW2
