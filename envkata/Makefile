
# Fortran Compiler
FC  = ifort
F90 = ifort

# Extra flags to give to the Fortran compiler.
FFLAGS   = -w
F90FLAGS = -w

all: removedir createdir envkata

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

MODS_EXKATAMPI = common.o exkataparameter.o exkataif.o
OBJS_EXKATAMPI = tbpot.o tbtube.o tbtubecb.o tubpar.o  \
	check.o cutline.o cutoff.o dia.o diacb.o dsort.o \
	eg3.o exh.o dgaussls.o findcon.o findkc.o findkr.o \
	ftdia.o ftrans.o ftranscb.o gcd.o iniga-dia.o inik.o \
	mukrg.o num2str.o numvhs.o qline.o qpencal.o \
	tbcoord.o wfcoe.o wfcoecb.o  wfcoecbi.o wfline.o \
	wflinecb.o wflinecbi.o select.o itoc.o

removedir:
	rm -Rf data-pol data-qex data-wf data-kapp

createdir:
	mkdir data-pol; mkdir data-qex; mkdir data-wf; mkdir data-kapp

envkata: $(MODS_EXKATAMPI) $(OBJS_EXKATAMPI)
	mpif90 $(F90FLAGS) $(MODS_EXKATAMPI) $(OBJS_EXKATAMPI) \
	arguments-kappa.f90 envkata.f90 $(LAPACK) -o envkata.out

clean:
	rm -f *.a *.exe *.mod *.o *.out *~
