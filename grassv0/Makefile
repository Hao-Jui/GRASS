PROG    := a.out
FXX     := gfortran
FFLAGS  := -O3 -funroll-loops -march=native -cpp -Wconversion -I. #-fbacktrace

SRC_CORE   := src/core
SRC_THEORY := src/theory

SOURCES := \
  $(SRC_CORE)/para_mod.f90 \
  $(SRC_CORE)/simpson_mod.f90 \
  $(SRC_CORE)/toolkit_mod.f90 \
  $(SRC_THEORY)/spin_gr.f90 \
  $(SRC_THEORY)/spin_st.f90 \
  $(SRC_THEORY)/rotation_dispatch.f90 \
  $(SRC_CORE)/newton_mod.f90 \
  $(SRC_CORE)/newton1d_mod.f90 \
  $(SRC_CORE)/eos.f90 \
  $(SRC_CORE)/grid.f90 \
  $(SRC_CORE)/NAG.f90 \
  $(SRC_CORE)/sphere.f90 \
  $(SRC_CORE)/brent.f90 \
  $(SRC_CORE)/restart_read.f90 \
  $(SRC_CORE)/analysis_v2.f90 \
  $(SRC_CORE)/shoot_v2.f90 \
  src/main.f90

OBJECTS := $(SOURCES:.f90=.o)

.PHONY: all clean

all: $(PROG)

$(PROG): $(OBJECTS)
	$(FXX) $(FFLAGS) -o $@ $(OBJECTS)

%.o: %.f90
	$(FXX) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJECTS) $(PROG) *.mod
