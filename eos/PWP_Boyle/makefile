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
FFLAGS = -g -O3
OMPFLAGS = -fopenmp

# suffix rule
.SUFFIXES: .F90 .o

# source file
SRC = gpp.F90
OBJS = $(SRC:%.F90=%.o)

.PHONY: all
all: $(PROG)

# 1st rule
$(PROG): $(OBJS)
	$(FXX) $(FFLAGS) -o $(PROG) $(OMPFLAGS) $(OBJS)

.F90.o:
	$(FXX) $(FFLAGS) $(OMPFLAGS) -c $<

# clean rule
.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
