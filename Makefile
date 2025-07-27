
#MODEL
MODEL = -DSCAL 

#COMPILER
COMP = mpif90 $(OPT)

COMP_SERIAL = gfortran


# OPT = -g -fbacktrace  -fbounds-check -Waliasing -Wunderflow -Wsurprising -fbacktrace -fcheck=all -Wall -fcheck=all
#-O3 -funroll-loops -finline-functions -ftree-vectorize -fopt-info-vec-optimized -ffast-math -fwhole-file # -march=native 
OPT = -O2

#FLAG TO COMPILER (-cpp for pre processing of directives, INTEL -fpp)
FLAG =  -cpp  $(MODEL)

#LIBRARY PATHS
PATH2FFTW = $(FFTW_LIB)
PATH2OPENBLAS = $(OPENBLAS_LIB)

#LINK LIBRARIES
LINK_LIB = -L$(PATH2FFTW) -lfftw3_mpi -lfftw3 -lm -L$(PATH2OPENBLAS) -lopenblas

#INCLUDE PATH
INC = -I$(FFTW_INCLUDE) -I$(OPENBLAS_INCLUDE)


#OBJECT FILES USED BY POST
MODULES =  variables.o in_out.o transforms.o physics.o solvers.o  setup.o


#MAIN OBJECT
MAINOBJ = convflow.o

#EXECUTABLE
CONVFLOW:		$(MODULES) $(MAINOBJ);	 $(COMP) $(FLAG) -o CONVFLOW.x $(MODULES) $(MAINOBJ) $(LINK_LIB) $(INC)


#COMPILE MODULES for POST
variables.o:		variables.f90;		$(COMP) $(FLAG) -c -o variables.o variables.f90
in_out.o:			in_out.f90;			$(COMP) $(FLAG) -c -o in_out.o in_out.f90
transforms.o:		transforms.f90;		$(COMP) $(FLAG) -c -o transforms.o transforms.f90 $(LINK_LIB) $(INC)
solvers.o:			solvers.f90;		$(COMP) $(FLAG) -c -o solvers.o solvers.f90
physics.o:			physics.f90;		$(COMP) $(FLAG) -c -o physics.o physics.f90
setup.o:			setup.f90;			$(COMP) $(FLAG) -c -o setup.o setup.f90

#POST-PROCESSING: 
convflow.o:		main.f90; 	$(COMP) $(FLAG) -c -o  convflow.o main.f90 $(LINK_LIB) $(INC)

# CLEAN
clean:
	rm -f *.o *.mod CONVFLOW.x