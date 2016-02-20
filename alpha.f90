PROGRAM ALPHA
	IMPLICIT NONE
	INTEGER :: IOS
	INTEGER :: A, B, C, D

	PRINT *, "PROGRAM ALPHA"
	PRINT *, "PRESS ANY KEY TO CONTINUE"
	READ (*, *)

	OPEN ( unit=10, file="alpha_in.txt" )

	READ ( 10, *, iostat=IOS ) A, B, C
	PRINT *, "READ: 1"
	PRINT *, "ATTEMPTING TO READ FIRST 3 NUMBERS IN THE FIRST LINE"
	PRINT *, "NUMBERS ARE:", A, B, C

	IF ( IOS .eq. 0 ) THEN
		READ ( 10, *, iostat=IOS ) A, B
		PRINT *, "READ: 2"
		PRINT *, "NEXT PARTIAL LINE IS:", A, B
		!results - did not work as expected
		!        - after this read finishes, the file pointer
		!        - advances to the next line whether or not it
		!        - all of the stuff on the prev line
		PRINT *, "SKIPPED LAST ITEM. NEXT READ STARTS ON NEXT LINE"
		PRINT *, ""
	ENDIF

	IF ( IOS .eq. 0 ) THEN
		READ ( 10, *, iostat=IOS ) A, B, C, D
		PRINT *, "READ: 3"
		PRINT *, "ALL OF LN 3 AND PARTIAL NEXT LN 4:", A, B, C, D
		!results - when contents of current line were all read
		!        - the program advances automatically to the next line
	ENDIF

	CLOSE ( 10 )

	PRINT *, ""
	PRINT *, "END"
END PROGRAM
