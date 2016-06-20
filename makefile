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
	LIBS = -llapack -L/opt/local/lib
	CFLAGS = 
	LFLAGS = $(DEBUGFLAG) $(LIBS)
endif
ifeq ($(UNAME), Linux)
	LIBS = -llapack
	CFLAGS =
	LFLAGS = $(DEBUGFLAG) $(LIBS)
endif

OBJECTS = common.o io.o linear_solver.o
MODS = $(OBJECTS:.o=.mod)

$(EXECNAME) : $(OBJECTS)
	$(GF) $(LFLAGS) $(OBJECTS) main.f90 -o $(EXECNAME)

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
