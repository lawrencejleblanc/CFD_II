Module utility
	contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: coefficient_construct
! Inputs: A(coefficient matrix), source(source vector, T(Temperature vector)
! Outputs: A, source
! Description: Subroutine reads in a 1-D vector of temperatures and stores them
!	       in a 2-D matrix to be used in LU decomposition. Subroutine also
!	       reads in and creates the source vector for LU decomposition.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE coefficient_construct(A, b, Ap, Ae, An, Aw, As, q_gen, source_size, ntype, stype, wtype, etype, Wbc, Sbc, Nbc, Ebc)

	IMPLICIT NONE

	! Variable declaration
	REAL*8, INTENT(IN) :: q_gen, Ap, Ae, An, Aw, As, Wbc, Sbc, Nbc, Ebc
	INTEGER, INTENT(IN) :: xsize, ysize
	INTEGER :: i, j				!Loop counters
	! Source size represents the size of the source matrix for LU decomposition as well as the
	! the number of internal unknown temperatures in the domain.
	REAL*8, INTENT(OUT) :: source_size
	
!!!!!!!!WHAT SIZE DO YOU BRING IN????????
	REAL*8, dimension(:,:), allocatable, INTENT(INOUT) :: A
	REAL*8, allocatable, INTENT(INOUT) :: b(:) 
	REAL*8 :: squarert
	
	! Compute number of unknowns to be solved for
	source_size = (xsize*ysize) - (2*xsize) - ((ysize-2)*2)
	
	squarert = SQRT(source_size)
	
	! Create coefficient and source matrices
	allocate(A(0:source_size-1,0:source_size-1))
	allocate(b(0:source_size-1))

	DO i = 0, source_size - 1, 1
		DO j = 0, source_size - 1, 1
			! If all dirichlet B.C. 
			IF (wtype.eq.1 .AND. stype.eq.1 .AND. etype.eq.1 .AND. ntype.eq.1) THEN 
				IF(i .eq. j) THEN
					A(i,j) = Ap
				ELSEIF(j .eq. i+1 .AND. j .ne. squarert .AND. j .ne. (source_size - squarert)
					A(i,j) = Ae
				ELSEIF(j .eq. i-1 .AND. i .ne. squarert .AND. i .ne. (source_size - squarert)
					A(i,j) = Aw
				ELSEIF(j .eq. i+3) THEN
					A(i,j) = An
				ELSEIF(j .eq. i-3) THEN
					A(i,j) = As
				ELSE 
					A(i,j) = 0
				ENDIF
			ENDIF  
		ENDDO
	ENDDO


	RETURN
	END SUBROUTINE
END MODULE
