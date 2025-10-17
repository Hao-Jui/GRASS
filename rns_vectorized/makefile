# program file name
PROG = a.out
# compiler
FXX = gfortran
#FXX = ifort
FFLAGS = -O3 -ffast-math -funroll-loops -march=native -cpp -Wconversion #-fopenmp #-fbacktrace -g --warn-all -fmax-stack-var-size=10

# suffix rule
.SUFFIXES: .f90 .c .o .mod
# source file
SRC = para_mod.f90 simpson_mod.f90 newton_mod.f90 toolkit_mod.f90 eos.f90 grid.f90 NAG.f90 set_disk.f90 \
	sphere.f90 spin.f90 const_j.f90 uryu.f90 analysis_v2.f90 restart_read.f90 \
	brent.f90 shoot_v2.f90 main.f90

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

