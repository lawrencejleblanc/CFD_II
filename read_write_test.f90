Program read_write_test
	IMPLICIT NONE
	INTEGER :: file_integer, i
	CHARACTER :: stringstring(20)

	OPEN (unit = 1, file = "input.dat")

	!READ (1, *) file_integer
	!DO i = 1, 1
	READ (1, "(A)") stringstring
	PRINT *, "string\n", stringstring
		!IF (stringstring(1) .ne. "#") THEN
		!	READ( stringstring, *) file_integer
		!	print *, "found integer", file_integer
		!END IF
	!END DO

	CLOSE (1)

	PRINT *, "i found this in the file", file_integer

END PROGRAM read_write_test
