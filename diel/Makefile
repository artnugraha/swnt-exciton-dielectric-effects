
# Fortran Compiler
FC  = ifort
F90 = ifort

all: findeii

# .f90: free-form Fortran 95/90 source files
%.o: %.f90
	$(F90) -c $(F90FLAGS) $< -o $@

# .f: fixed-form  Fortran source files
%.o: %.f
	$(FC) -c $(FFLAGS) $< -o $@

# .for: fixed-form  Fortran source files
%.o: %.for
	$(FC) -c $(FFLAGS) $< -o $@

# .ftn: fixed-form  Fortran source files
%.o: %.ftn
	$(FC) -c $(FFLAGS) $< -o $@


MODS_PROG = mdatEii.o
OBJS_PROG = eiiFind.o

findeii: $(MODS_PROG) $(OBJS_PROG)
	ifort $(MODS_PROG) $(OBJS_PROG) testeii.f90 -o testeii.out
