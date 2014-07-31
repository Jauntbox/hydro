#################################################################

# STEP 1: set the LAPACK directory
 
FC = gfortran
FCbasic = -fno-range-check -fmax-errors=100 $(SPECIAL_FC_FLAGS) -fprotect-parens -fno-sign-zero
FCfree = -ffree-form $(FC_free_preprocess)
FCopenmp = -fopenmp

LOAD_LAPACK = -L/Users/Kevin/Downloads/lapack-3.5.0 -llapack -lblas
LOAD_FFTW = -lfftw3
INCLUDE_FFTW = -I/Users/Kevin/Downloads/fftw-3.3.4/api
LOAD_PGPLOT = -L$(PGPLOT_DIR) -lpgplot
LOAD_OTHER = 

#################################################################

# STEP 2: build the program

#For compiling the tridiagonal tester:
#PROG = test_tridiag
#PROG_OBJS = test_tridiag.o

#For compiling the FFTW tester:
#PROG = test_fftw
#PROG_OBJS = test_fftw.o

#For compiling the 2D hydro code:
PROG = test_hydro
PROG_OBJS = test_hydro.o

PROG_DIR = .

$(PROG) : $(PROG_OBJS)
	$(FC) $(FCbasic) $(FCopenmp) -o$(PROG_DIR)/$(PROG) $(PROG_OBJS) $(LOAD_LAPACK) $(LOAD_FFTW) $(LOAD_PGPLOT) $(LOAD_OTHER)

#################################################################

MY_FC_FLAGS = $(FCfree)
SRC_DIR = src

%.o: $(SRC_DIR)/%.f
	$(FC) $(FCbasic) $(FCopenmp) $(MY_FC_FLAGS) $(INCLUDE_FFTW) -c $<
