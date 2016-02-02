#see http://genius2k.is-programmer.com/posts/40301.html
#see http://www-k12.atmos.washington.edu/~ovens/junk/ifc7docs/f_ug/fced_mod.htm
#see http://faculty.washington.edu/rjl/classes/am583s2013/notes/fortran_modules.html

FC = gfortran
objects = HW1_Solution.o
MOD = solvers.o

%.o: %.f90
	$(FC) -c $<

%.mod: %.f90
	$(FC) -c $(MOD)

example.exe: $(MOD) $(objects)
	$(FC) $(objects) $(MOD) -o example.exe

clean:
	rm -f *.o *.mod
