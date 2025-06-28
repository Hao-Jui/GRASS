# program file name
PROG = rns
# compiler
FXX = gfortran
FFLAGS = -O3 -cpp -Wconversion # -fopenmp -fbacktrace -g --warn-all
#FXX = ifort
#FFLAGS = -O3 -g -cpp -mcmodel=large

# suffix rule
.SUFFIXES: .f90 .c .o .mod
# source file
SRC = para_mod.f90 tools.f90 eos.f90 grid.f90 NAG.f90 \
	sphere.f90 spin.f90 spin_massless.f90 const_j.f90 uryu.f90 analysis_v2.f90 restart_read.f90 \
	brent.f90 J_seq.f90 shoot_unstable.f90 shoot_v2.f90 main.f90

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

