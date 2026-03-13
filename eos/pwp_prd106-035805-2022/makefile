# Makefile

# program file name
PROG = a.out
# compiler
#FXX = ifort
# compiler option
#FFLAGS = -warn nousage -m64 -xHOST -O3 -ip -prec-div -pg
#FFLAGS = -warn nousage -m64 -xHOST -O3 -ipo -prec-div -pg
#FFLAGS = -warn nousage -m64 -xHOST -O3 -ipo -prec-div
#FFLAGS = -cm -w -xHOST -O3 -ip -prec-div
#FFLAGS = -cm -w -O0 -ip
OMPFLAGS = -qopenmp #-qopenmp-report2

# for gfortran
FXX = gfortran
FFLAGS = -O3 -Wconversion -cpp #-fbacktrace -g
OMPFLAGS = -fopenmp

# suffix rule
.SUFFIXES: .f90 .o

# source file
SRC = main.f90

OBJS = $(SRC:%.f90=%.o)
LIBS = -L/Applications/mesasdk/math-slots/default/lib -llapack -lblas

.PHONY: all
all: $(PROG)

# 1st rule
$(PROG): $(OBJS)
	$(FXX) $(FFLAGS) -o $(PROG) $(OMPFLAGS) $(OBJS) $(LIBS)

.f90.o:
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $<

# clean rule
.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
