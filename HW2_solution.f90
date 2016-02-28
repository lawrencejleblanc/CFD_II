!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Name: Lawrence LeBlanc
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
	Use Solvers
        Use matrix_coefficients 
	IMPLICIT NONE

	! Variable declaration
	REAL*8, dimension(:,:), allocatable :: A
	REAL*8, allocatable :: b(:), x(:), T(:)
	INTEGER :: i, j, wtype, etype, stype, ntype, OpenStatus, solver = 0
	REAL*8 :: dx, dy, q_gen, Wbc, Ebc, Nbc, Sbc, As, An, Ae, Ap, Aw, y, xx
	INTEGER :: xsize, ysize, area, source_size

	! Opens input file and output tecplot file
	OPEN (UNIT = 1, FILE = "input.dat", IOSTAT = OpenStatus)
	OPEN (UNIT = 2, FILE = "output.dat", ACTION = "WRITE")
	IF (OpenStatus .gt. 0) THEN
		STOP "** Cannot open input file **"
	ENDIF

	! Reads values from input file
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
	WRITE (*,8) "xsize =", xsize, "ysize =", ysize, "dx =", dx, "dy =", dy, "q_gen =", q_gen 
	8 FORMAT(1X, A7, I4, /, 1X, A7, I4, /, 1X, A4, F5.2, /, 1X, A4, F5.2, /, 1X, A7, F6.2)
	WRITE (*,7) "Wbc =", Wbc, "wtype =", wtype, "Nbc =", Nbc, "ntype =", ntype, "Sbc =", Sbc
	7 FORMAT(1X, A5, F5.2, /, 1X, A7, I1, /, 1X, A5, F5.2, /, 1X, A7, I1, /, 1X, A5, F5.2)
	WRITE (*,10) "stype =", stype, "Ebc =", Ebc, "etype =", etype, "solver =", solver
	10 FORMAT(1x, A7, I1, /, 1X, A5, F5.2, /, 1x, A7, I1, /, 1X, A8, I1)
	
	! Determines total number of mesh points in the domain
	area = xsize * ysize
	
	! Calculate matrix coefficients
	As = 1/(dy**2)
	An = As
	Aw = 1/(dx**2)
	Ae = Aw
	Ap = (-2/(dx**2))-(2/(dy**2))
	
	print*, ""

	! Call solver suboutines
	! Condition check for Gauss-Seidel solver
	IF (solver .eq. 1)THEN
		allocate(T(1:xsize*ysize))
		! Initialize temperature array to 0. These values are used as initial T_prev
		DO i = 1, xsize*ysize, 1
			T(i) = 0
		ENDDO
		
		CALL gauss_seidel(xsize, ysize, An, As, Aw, Ae, Ap, q_gen, Wbc, Ebc, Nbc, Sbc, T, area, wtype, etype, stype, ntype)
		
		! Write the header to the output tecplot file
		WRITE(2,*) 'Variables = "X" "Y" "TEMPERATURE"'
		WRITE(2,9) "zone I=", xsize, "J=", ysize, " SOLUTIONTIME=", 0.00, "F=POINT"
		9 FORMAT(2X, A7, I4, 2X, A3, I4, 2X, A14, F6.3, 2X, A7)
		
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
		
	! Condition check for LU Decomposition solver
	ELSEIF (solver .ne. 1)THEN
		! Matrice construction subroutine
		CALL coefficient_construct(A, b, Ap, Ae, An, Aw, As, q_gen, &
                &source_size, ntype, stype, wtype, etype, Wbc, Sbc, Nbc, Ebc, xsize, ysize)
              
		!DO i = 0, source_size - 1
                 !       DO j = 0, source_size - 1
                  !              print*, A(i,j), "", i, "", j
                   !     ENDDO
		!ENDDO
		!DO i = 0, source_size - 1
		!	print*,b(i)
		!ENDDO 
		! LU Decomposition solver
		CALL LUdecomp(source_size, A)
		CALL LUsolve(source_size, A, b)

		! Write the header to the output tecplot file
		WRITE(2,*) 'Variables = "X" "Y" "TEMPERATURE"'
		WRITE(2,13) "zone I=", xsize-2, "J=", ysize-2, " SOLUTIONTIME=", 0.00, "F=POINT"
		13 FORMAT(2X, A7, I4, 2X, A3, I4, 2X, A14, F6.3, 2X, A7)
		
		! Write the x-position, y-position, and temperature values to the tecplot file
		DO j = 0, source_size - 1, 1
			xx = (MOD(j, (xsize-2)) + 1)*dx
			IF (xx .lt. 0) THEN
				xx = xx + (xsize*dx)
			ENDIF
			y = ((j)/(xsize-2)+1)*dy
			WRITE(2,14) xx, y, b(j)
			14 FORMAT(F10.5, 3X, F10.5, 3X, F10.5)
		ENDDO
	ENDIF

	! Close input and output files
	CLOSE (1)
	Close (2)

END PROGRAM HW2
