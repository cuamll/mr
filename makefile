# maggs-rossetto makefile
GF = gfortran
EXECNAME = mr_test

# libraries in different directories on Mac
UNAME = $(shell uname)

DEBUG = 1
ifeq ($(DEBUG), 1)
	DEBUGFLAG = -g
else
	DEBUG FLAG =
endif

ifeq ($(UNAME), Darwin)
	LIBS = -llapack -lfftw3 -L/opt/local/lib -L/usr/local/lib
	CFLAGS = 
	LFLAGS = $(DEBUGFLAG) $(LIBS)
endif
ifeq ($(UNAME), Linux)
	LIBS = -llapack -lfftw3
	CFLAGS =
	LFLAGS = $(DEBUGFLAG) $(LIBS)
endif

OBJECTS = fftw.o common.o io.o linear_solver.o setup.o
MODS = $(OBJECTS:.o=.mod)

$(EXECNAME) : $(OBJECTS)
	$(GF) $(LFLAGS) $(DEBUGFLAG) $(OBJECTS) main.f90 -o $(EXECNAME)

fftw.o : fftw.f03
	$(GF) $(DEBUGFLAG) -c -I/usr/local/include $<

%.o : %.f90
	$(GF) $(DEBUGFLAG) -c $<

# add this if need be, can't be arsed to fuck around rn
.PHONY: clean clean_obj clean_mod

clean:
	\rm -f $(OBJECTS) $(MODS) $(EXECNAME)

clean_obj:
	\rm -f $(OBJECTS)

clean_mods:
	\rm -f $(MODS)
