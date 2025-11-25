PROG    := a.out
FC      = gfortran
FFLAGS  ?= -O2 -march=native -cpp -ftree-vectorize -I.
LDFLAGS ?=
LIBS    ?= -llapack -lblas
SRC_TOOL := src/tool
SRC_CORE   := src/core
SRC_THEORY := src/theory

SOURCES := \
  $(SRC_TOOL)/ad_mod.f90 \
  $(SRC_TOOL)/simpson_mod.f90 \
  $(SRC_TOOL)/nag_compat_mod.f90 \
  $(SRC_TOOL)/brent.f90 \
  $(SRC_TOOL)/gauss_quad.f90 \
  $(SRC_TOOL)/cheb_mod.f90 \
  $(SRC_CORE)/para_mod.f90 \
  $(SRC_TOOL)/toolkit_mod.f90 \
  $(SRC_CORE)/miscellaneous.f90 \
  $(SRC_CORE)/set_disk.f90 \
  $(SRC_THEORY)/spin_helper.f90 \
  $(SRC_THEORY)/spin.f90 \
  $(SRC_THEORY)/spin_gr.f90 \
  $(SRC_THEORY)/spin_massless.f90 \
  $(SRC_THEORY)/rotation_dispatch.f90 \
  $(SRC_CORE)/newton_mod.f90 \
  $(SRC_CORE)/newton1d_mod.f90 \
  $(SRC_CORE)/eos.f90 \
  $(SRC_CORE)/grid.f90 \
  $(SRC_CORE)/sphere.f90 \
  $(SRC_CORE)/restart_read.f90 \
  $(SRC_CORE)/analysis_v2.f90 \
  $(SRC_CORE)/shoot_v2.f90 \
  $(SRC_CORE)/MRcurve.f90 \
  src/main.f90

OBJECTS := $(SOURCES:.f90=.o)

.PHONY: all clean

all: $(PROG)

$(PROG): $(OBJECTS)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) $(OBJECTS) $(PROG) *.mod
