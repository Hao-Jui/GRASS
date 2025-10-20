# program file name
PROG = a.out
# compiler
FXX = gfortran
#FXX = ifort
# comiler option
FFLAGS= -O3 -fopenmp -g -cpp #--warn-all
#FFLAGS = -warn nousage -m64 -xHOST -O3 -ip -prec-div -pg
# suffix rule
.SUFFIXES: .f90 .c .o .mod
# source file
SRC = make_gamma_tab.f90

OBJS = $(SRC:%.f90=%.o)

.PHONY: all
all: $(PROG)

# 1st rule
$(PROG): $(OBJS)
	$(FXX) $(FFLAGS) -o $(PROG) $(OBJS)

%.o: %.f90
	$(FXX) -c -o $@ $< $(FFLAGS)

%.mod: %.f90 %.o
	@test -f $@ || $(FXX) -c -o $(@:.mod=.o) $< $(FFLAGS)

# clean rule
.PHONY: clean

clean:
	$(RM) $(OBJS) fort* *.mod

# Module Dependence

