MODULE matrix_coefficients
        CONTAINS      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Name: coefficient_construct							 !
! Inputs: Coefficient matrix, source vector, B.C. types and values, mesh         !
!	  size in x and y direction, heat generation over thermal conductivity,  !
!	  and coefficients for present, north, south, east, and west temperatures!    
! Outputs: A, source, and size of the source vector				 !
! Description: Subroutine reads in a 1-D source vector, 2-D coefficient array    !
!	       and the associated boundary conditions and creates the coefficient!
!	       matrix and source vector to be inputted into LU decomposition for !
!	       solving for the interior temperature points of the mesh.		 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE coefficient_construct(A, b, Ap, Ae, An, Aw, As, &
&q_gen, source_size, ntype, stype, wtype, etype, Wbc, Sbc, Nbc, Ebc, xsize, ysize)
	IMPLICIT NONE

	! Variable declaration
	REAL*8, INTENT(IN) :: q_gen, Ap, Ae, An, Aw, As, Wbc, Sbc, Nbc, Ebc
	INTEGER, INTENT(IN) :: xsize, ysize
	INTEGER :: i, j                              !Loop counters
	INTEGER, INTENT(IN) :: ntype, stype, wtype, etype
	! Source size represents the size of the source matrix for LU decomposition as well as the
	! the number of internal unknown temperatures in the domain.
	INTEGER, INTENT(OUT) :: source_size
	
	! Allocate coefficient and source matrices
	REAL*8, dimension(:,:), allocatable, INTENT(INOUT) :: A
	REAL*8, allocatable, INTENT(INOUT) :: b(:) 
	
	! Compute number of unknowns to be solved for
	source_size = (xsize-2)*(ysize-2)
	! Create coefficient and source matrices
	allocate(A(0:source_size-1,0:source_size-1))
	allocate(b(0:source_size-1))

	! Set coefficient matrix values
	DO i = 0, source_size - 1, 1 
		DO j = 0, source_size - 1, 1 
			! If all dirichlet B.C. 
			IF (wtype.eq.1 .AND. stype.eq.1 .AND. etype.eq.1 .AND. ntype.eq.1) THEN
				! Sets Ap values 
				IF(i .eq. j) THEN
					A(i,j) = Ap
				! Sets Ae values
				ELSEIF(j .eq. i+1 .AND. MOD(j,(xsize-2)) .ne. 0) THEN
					A(i,j) = Ae
				! Sets Aw values
				ELSEIF(j .eq. i-1 .AND. MOD(i,(xsize-2)) .ne. 0) THEN
					A(i,j) = Aw
				! Sets An values
				ELSEIF(j .eq. i+(xsize-2)) THEN
					A(i,j) = An
				! Sets As values
				ELSEIF(j .eq. i-(xsize-2)) THEN
					A(i,j) = As
				! Zeroes all remaining values
				ELSE 
					A(i,j) = 0
				ENDIF
			
			! If at least 1 Neumann B.C.
			ELSE
				! Sets Ap values
				IF(i .eq. j) THEN
					A(i,j) = Ap
					IF (MOD(j,(xsize-2)) .eq. 0 .AND. wtype .eq. 0) THEN
						A(i,j) = A(i,j) + Aw
					ENDIF
					IF (i .lt. (xsize-2) .AND. stype .eq. 0) THEN
						A(i,j) = A(i,j) + As
					ENDIF
					IF (i .ge. (source_size - (xsize-2)) .AND. ntype .eq. 0) THEN
						A(i,j) = A(i,j) + An
					ENDIF
					IF (MOD((j+1),(xsize-2)) .eq. 0 .AND. etype .eq. 0) THEN
						A(i,j) = A(i,j) + Ae
					ENDIF
				! Sets Ae values
				ELSEIF(j .eq. i+1 .AND. MOD(j,(xsize-2)) .ne. 0) THEN
					A(i,j) = Ae
				! Sets Aw values
				ELSEIF(j .eq. i-1 .AND. MOD(i,(xsize-2)) .ne. 0) THEN
					A(i,j) = Aw
				! Sets An values
				ELSEIF(j .eq. i+(xsize-2)) THEN
					A(i,j) = An
				! Sets As values
				ELSEIF(j .eq. i-(xsize-2)) THEN
					A(i,j) = As
				! Zeroes all remaining values
				ELSE
					A(i,j) = 0
				ENDIF
			ENDIF   
		ENDDO
	ENDDO
	
	! Initialize source vector values to q_gen values
	DO i = 0, source_size -1, 1
		IF (q_gen .eq. 0) THEN
			b(i) = 0
		ELSE
			b(i) = -q_gen
		ENDIF
	ENDDO
	
	! Set source vector values
	DO i = 0, source_size -1, 1
		
		! If all dirichlet B.C.
		IF (wtype.eq.1 .AND. stype.eq.1 .AND. etype.eq.1 .AND. ntype.eq.1) THEN
			IF (MOD(i,(xsize-2)) .eq. 0) THEN
				b(i) = b(i) - (Wbc*Aw)
			ENDIF
			IF (i .lt. (xsize-2)) THEN
				b(i) = b(i) - (Sbc*As)
			ENDIF
			IF (MOD((i+1),(xsize-2)) .eq. 0) THEN
				b(i) = b(i) - (Ebc*Ae)
			ENDIF
			IF (i .ge. (source_size - (xsize-2))) THEN
				b(i) = b(i) - (Nbc*An)
			ENDIF
			
		! If at least 1 Neumann B.C.
		ELSE
			IF (MOD(i,(xsize-2)) .eq. 0 .AND. wtype .eq. 1) THEN
				b(i) = b(i) -(Wbc*Aw)
			ENDIF
			IF (i .lt. (xsize-2) .AND. stype .eq. 1) THEN
				b(i) = b(i) -(Sbc*As)
			ENDIF
			IF (MOD((i+1),(xsize-2)) .eq. 0 .AND. etype .eq. 1) THEN
				b(i) = b(i) - (Ebc*Ae)
			ENDIF
			IF (i .ge. (source_size - (xsize-2)) .AND. ntype .eq. 1) THEN
				b(i) = b(i) -(Nbc*An)
			ENDIF
		ENDIF
	ENDDO
	
	RETURN
	END SUBROUTINE
END MODULE
