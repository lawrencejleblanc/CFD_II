!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Name: Lawrence LeBlanc
!Date: March 8, 2016
!Program name: heat_eq.f90
!Program description: This program solves the 1D heat equation explicitly and
!                     implicitly. For implicit, uses a Gauss-Seidel solver. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM heat_eq
	IMPLICIT NONE
	
	! Variable declaration
	INTEGER :: method            ! Determines time step numerical method (Implicit, Explicit, Etc.)
	REAL*8  :: dt                ! Time step size
	REAL*8 :: L                 ! Length of problem
	REAL*8, allocatable :: T(:)       ! 1-D vector for the variable to be solved (next time step)
	REAL*8, allocatable :: T_old(:)   ! 1-D vector for previous value (current time step)
	REAL*8, allocatable :: T_prev(:)  ! Prevous function vector for Gauss-Seidel
        REAL*8, allocatable  :: T_exact(:)                ! Exact solution
        REAL*8  :: alpha                  ! Thermal diffusivity for heat equation
	INTEGER :: i, n, m                ! Loop counters
	REAL*8  :: Wbc, Ebc               ! Boundary conditions
	INTEGER :: num_dt = 100           ! Number of time steps
	REAL*8  :: dx                     ! element size
        INTEGER :: num_points= 10         ! Number of mesh points 
	INTEGER :: converge = 0           ! Gauss-Seidel convergence check
	REAL*8  :: criteria = 0.00001     ! Gauss-Seidel convergence criteria               
	REAL*8  :: Fo                     ! Fourier number
        INTEGER :: iter                   ! iteration counter for Gauss-Seidel
	REAL*8  :: summation              ! exact solution individual summation term
        REAL*8  :: running_sum            ! Total summation term
        REAL*8  :: T_side = 10            ! B.C.
	REAL*8  :: T_initial = 1            ! Initial solution
        REAL*8  :: pi = 3.14159           ! pi
        REAL*8  :: x                      ! x-position along the domain for
                                          ! writing out to tecplot
	REAL*8  :: del_t                  ! Variable for maximum time step size
	REAL*8  :: A, B, C                ! Exact solution coefficients
	! Open output Tecplot file
	OPEN (UNIT = 2, FILE = "output.dat", ACTION = "WRITE")

	! Input general solution parameters
	PRINT*, "Enter a value for alpha."
	READ*, alpha
	PRINT*, "Enter a value for domain length."
	READ*, L
        dx = L / (num_points-1)
	del_t = (dx**2)/(2*alpha)
	PRINT*, "Enter a value for time step size that is less than or equal to", del_t
	READ*, dt

	! Initialize function vectors to initial solution
	allocate(T(0:num_points-1))
	allocate(T_old(0:num_points-1))
	allocate(T_prev(0:num_points-1))
        allocate(T_exact(0:num_points-1))
	DO i = 0, num_points-1
		T(i) = T_initial
		T_old(i) = T_initial
		T_prev(i) = T_initial
                T_exact(i) = T_initial
	ENDDO
	
	! Determine numerical method to use.
	PRINT*, "Determine time step numerical method."
	PRINT*, "1) Implicit"
	PRINT*, "2) Explicit"
	PRINT*, "3) Crank-Nicolson"
	READ* , method
	
	! Explicit Euler Method
	IF (method .eq. 2) THEN
		
		DO n = 0, num_dt
                        ! Write tecplot header
                        WRITE(2,*) 'Variables = "X" "T" "T_exact"'
                        WRITE(2,*) "zone I=", num_points, "SOLUTIONTIME=", n*dt, "F=POINT"
                        !9 FORMAT(2X, A7, I5, 2X, A13, F8.6, 2X, A7)
			T(0) = T_side
			T_old(0) = T_side
                        T_exact(0) = T_side
			T(num_points-1) = T_side
			T_old(num_points-1) = T_side
                        T_exact(num_points-1) = T_side
			
			! Calculate Fourier Number
			Fo = (alpha*dt*n)/(L**2)
			
			! Loop over internal points
			DO i = 1, num_points - 2
				T(i) = T_old(i) + alpha*dt*((T_old(i+1)-2*T_old(i)+T_old(i-1))/(dx**2))
				
				! Calculate Fourier number
				!Fo = (alpha*dt*n)/(L**2)
                                running_sum = 0
                                ! Calculating the sum term of exact solution
                                DO m = 1, 5
					A = (1-cos(m*pi))/(m*pi)
					B = sin((m*pi*dx*i)/L)
					C = -((m*pi)**2)*Fo
					summation = (EXP(C))*A*B
				       ! summation = EXP(-((m*pi)**2)&
                                       ! &*Fo)*((1-cos(m*pi))/m*pi)*sin((m*pi*dx*i)/L)

                                        running_sum = running_sum + summation	
					
                                ENDDO
				PRINT*, running_sum, "", Fo
				PRINT*,""
                                ! Calculate exact solution
                                T_exact(i) = T_side + 2*(T_initial - T_side)*running_sum
 
			ENDDO
    
			DO i = 0, num_points - 1
				T_old = T(i)
                                x = i*dx
                                WRITE(2,*) x, T(i), T_exact(i)
			ENDDO
		ENDDO
	
	! Implicit Euler Method
!	ELSEIF (method .eq. 1) THEN
		! Define matrix coefficients
!		Ap = 1+((2*alpha*dt)/(dx**2))
!		Ae = (-dt*alpha)/(dx**2)
!		Aw = Ae
!		DO WHILE (converge .eq. 0 .AND. iter .lt. 10000)
!			converge = 1    ! Assume convergance
			! Loop over time
!			DO n = 0, total_time
				! Loop over space and solve with Gauss-Seidel
!				DO i = 0, L-1
!					IF (i = 0) THEN
!						T(i) = T_side
!					ELSEIF (i = L-1) THEN
!						T(i) = T_side
!					ELSE 
!						T(i) = (T_old(i) - Ae*T(i+1) -
 !                                               & Aw*T(i-1)) * (1 / Ap)
!					ENDIF
!					IF (ABS(T(i) - T_prev(i)) .gt. criteria) THEN
!						converge = 0
!					ENDIF
!				ENDDO
!				! Set previous n+1 values to n
!				DO i = 0, L
!					T(i) = T_prev(i)
!				ENDDO
!			ENDDO
!		ENDDO
	
	! Crank-Nicolson Method
	!ELSEIF (method .eq. 3) THEN
	
!	ELSE
!		PRINT*, "Invalid entry for numerical method."
	ENDIF
	
END PROGRAM heat_eq
