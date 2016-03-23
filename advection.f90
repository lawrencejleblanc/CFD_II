!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Name: Lawrence LeBlanc
!Date: March 8, 2016
!Program name: heat_eq.f90
!Program description: This program solves the 1D advection equation. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM advection

	IMPLICIT NONE
	
	! Variable declaration
	INTEGER :: c = 1                 ! Assumed 1 for this assignment
	REAL*8  :: dt, dx                ! Time step size and grid spacing size
	INTEGER :: num_dt                ! Number of time steps
	INTEGER :: scheme                ! Determine numerical scheme (Upwinding, MacCormack ,Etc.)  
	INTEGER :: L                     ! Number of grid points in the domain
	REAL*8  :: Wbc, Ebc              ! 1-D boundary conditions
	
	REAL*8, allocatable :: T(:), T_old(:), T_exact(:)  ! Function vectors
	
	INTEGER :: x                     ! distance x along domain
	
	! Open output Tecplot file
	OPEN (UNIT = 2, FILE = "output.dat", ACTION = "WRITE")

	! Write tecplot header
	WRITE(2,*) 'Variables = "X" "T" "T_exact"'
	WRITE(2,9) "zone I=", L, "SOLUTIONTIME=", num_dt*dt, "F=POINT"
	9 FORMAT(2X, A7, I5, 2X, A13, F5.3, 2X, A7)
	
	! Specify general paramters for solving
	PRINT*,"Enter a real value for time step."
	READ*, dt
	PRINT*, ""
	PRINT*, "ENTER a real value for grid element size."
	READ*, dx
	PRINT*,""
	PRINT*, "Enter value for number of time steps"
	READ*, num_dt
	PRINT*, ""
	PRINT*, "Enter a value for number of grid points."
	READ*, L
	PRINT*, ""
	PRINT*, "Select numerical scheme to solve advection equation with"
	PRINT*, "1) Upwinding"
	PRINT*, "2) MacCormack"
	PRINT*, "3) Lax-Friedrich"
	READ*, scheme
	
	! Initialize function vectors to 0
	allocate(T(0:L-1))
	allocate(T_old(0:L-1))
	DO i = 0, L-1
		T(i) = 0
		T_old(i) = 0
	ENDDO
	
	! Initialize step function
	DO i = L = 0, L / 4
		T(i) = 2
	ENDDO 
	
	! Determine B.C. values
	PRINT*, "Left and right edge B.C. are periodic."
	PRINT*, "Enter a value for the left and right B.C. (Enter 1 value)"
	READ*, Wbc
	Ebc = Wbc
	
	! Upwinding scheme
	IF (scheme .eq. 1) THEN
		DO n = 0, num_dt
			!T_prev(0) = Wbc
			!T_prev(xsize-1) = Wbc
			DO i = 0, L-1
				T(i) = T_prev(i) - c*dt*n*(T_prev(i) - T_prev(i-1))/dx
				!! INSERT ANALYTICAL SOLUTION
				x = i * dx
				WRITE(2,*) x, T(i), T_exact(i)
			ENDDO
			! Reset new value as old value
			DO i = 0, L
				T_prev(i) = T(i)
			ENDDO
		ENDDO
	
	! MacCormack Method
	ELSEIF (scheme .eq. 2) THEN
	!!! CODE CODE CODE!!!
	!!! CODE CODE CODE!!!
	
	! Lax- Friedrichs Method
	ELSEIF (scheme .eq. 3) THEN
		DO n = 0, num_dt
			!!! B.C.
			DO i = 0, L-1
				T(i) = 0.5*(T_prev(i+1)-T_prev(i-1))-c*dt*n*&
				&((T_prev(i+1)-T_prev(i-1))/(2*dx))
				!!! INSERT ANALYTICAL SOLUTION
				x = i * dx
				WRITE(2,*) x, T(i), T_exact
			ENDDO
			! Reset new value as old value
			DO i = 0, L
				T_prev(i) = T(i)
			ENDDO
		ENDDO
	
	ELSE
		PRINT*,"Invalid entry for numerical scheme"
	ENDIF
	

END PROGRAM advection
