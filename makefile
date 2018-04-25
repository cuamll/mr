# maggs-rossetto makefile
GF = mpif90
EXECNAME = mr2d
MOD_DIR = mod
OBJ_DIR = obj
CFLAGS = -J$(MOD_DIR) -std=f2003 -fopenmp

$(shell mkdir -p $(MOD_DIR))
$(shell mkdir -p $(OBJ_DIR))
VPATH = $(OBJ_DIR)
VPATH = $(MOD_DIR)

# libraries in different directories on Mac
UNAME = $(shell uname)
RV = $(shell git rev-parse --short HEAD)

DEBUG = 0
ifeq ($(DEBUG), 1)
	# DEBUGFLAGS = -g -pg -fbacktrace -fopenmp -fbounds-check \
	# 	     -ffpe-trap=invalid,zero,denormal,underflow,overflow
	DEBUGFLAGS = -g -pg -fbacktrace -fopenmp -fbounds-check \
		     -ffpe-trap=invalid,zero
else
	DEBUGFLAGS = -O2
endif

ifeq ($(UNAME), Darwin)
	LIBS = -llapack -L/opt/local/lib
	LFLAGS = $(DEBUGFLAGS) $(LIBS)
endif
ifeq ($(UNAME), Linux)
	LIBS = -llapack
	LFLAGS = $(DEBUGFLAGS) $(LIBS)
endif

SOURCES = common.f03\
	  input.f03\
	  linear_solver.f03\
	  output.f03\
	  setup.f03\
	  update.f03

OBJ_T1 = $(patsubst %.f03, %.o,$(SOURCES))
OBJ_T2 = $(notdir $(OBJ_T1))
OBJECTS = $(patsubst %.o, $(OBJ_DIR)/%.o,$(OBJ_T2))
MOD_T1 = $(patsubst %.f03, %.mod,$(SOURCES))
MOD_T2 = $(notdir $(MOD_T1))
MODS = $(patsubst %.o, $(MOD_DIR)/%.o,$(MOD_T2))

$(EXECNAME) : $(OBJECTS)
	echo "character(len=7), parameter :: revision = '$(RV)'" > rev.inc
	$(GF) $(CFLAGS)  $(LFLAGS) $(OBJECTS) main.f03 -o $(EXECNAME)

$(OBJ_DIR)/%.o : %.f03
	$(GF) $(CFLAGS) $(DEBUGFLAGS) -o $@ -c $<

# add this if need be, can't be arsed to fuck around rn
.PHONY: clean cleaner all remake clean_obj clean_mod

clean:
	\rm -f $(EXECNAME) revision.inc

cleaner:
	\rm -rf $(OBJ_DIR) $(MOD_DIR) $(EXECNAME) revision.inc

remake:
	cleaner

clean_obj:
	\rm -f $(OBJECTS)

clean_mods:
	\rm -f $(MODS)
