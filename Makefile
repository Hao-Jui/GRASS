PROG    := a.out
.DEFAULT_GOAL := all
FC      = gfortran
FFLAGS  ?= -O2 -march=native -ftree-vectorize -I.
LDFLAGS ?=
LIBS    ?= -llapack -lblas
BUILDDIR := build
OBJDIR   := $(BUILDDIR)/obj
MODDIR   := $(BUILDDIR)/mod
BINDIR   := $(BUILDDIR)/bin
TARGET   := $(BINDIR)/$(PROG)

FFLAGS  += -J$(MODDIR) -I$(MODDIR)
SRC_TOOL := src/tool
SRC_CORE   := src/core
SRC_THEORY := src/theory

SOURCES := \
  $(SRC_CORE)/precision_mod.f90 \
  $(SRC_TOOL)/ad_mod.f90 \
  $(SRC_TOOL)/simpson_mod.f90 \
  $(SRC_TOOL)/nag_compat_mod.f90 \
  $(SRC_TOOL)/brent.f90 \
  $(SRC_TOOL)/cheb_mod.f90 \
  $(SRC_CORE)/para_mod.f90 \
  $(SRC_TOOL)/toolkit_mod.f90 \
  $(SRC_CORE)/miscellaneous.f90 \
  $(SRC_CORE)/set_disk.f90 \
  $(SRC_THEORY)/relaxation_mod.f90 \
  $(SRC_THEORY)/rotational_law_mod.f90 \
  $(SRC_THEORY)/spin_helper.f90 \
  $(SRC_THEORY)/rotation_solver.f90 \
  $(SRC_CORE)/newton_mod.f90 \
  $(SRC_CORE)/newton1d_mod.f90 \
  $(SRC_CORE)/eos.f90 \
  $(SRC_CORE)/grid.f90 \
  $(SRC_CORE)/sphere.f90 \
  $(SRC_CORE)/restart_read.f90 \
  $(SRC_CORE)/analysis_v2.f90 \
  $(SRC_CORE)/scalar_burning.f90 \
  $(SRC_CORE)/starting_model.f90 \
  $(SRC_CORE)/shoot_v2.f90 \
  $(SRC_CORE)/MRcurve.f90 \
  src/main.f90

OBJECTS := $(patsubst %.f90,$(OBJDIR)/%.o,$(SOURCES))

# Ensure precision module is built before any source that imports it.
$(filter-out $(OBJDIR)/$(SRC_CORE)/precision_mod.o,$(OBJECTS)): $(OBJDIR)/$(SRC_CORE)/precision_mod.o
$(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o: $(OBJDIR)/$(SRC_TOOL)/ad_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/miscellaneous.o: $(OBJDIR)/$(SRC_TOOL)/nag_compat_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/set_disk.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/relaxation_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/rotational_law_mod.o: $(OBJDIR)/$(SRC_TOOL)/brent.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/spin_helper.o: $(OBJDIR)/$(SRC_TOOL)/nag_compat_mod.o $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/rotation_solver.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_helper.o
$(OBJDIR)/$(SRC_CORE)/newton_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver.o
$(OBJDIR)/$(SRC_CORE)/newton1d_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver.o
$(OBJDIR)/$(SRC_CORE)/eos.o: $(OBJDIR)/$(SRC_TOOL)/ad_mod.o $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/grid.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/sphere.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/restart_read.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/analysis_v2.o: $(OBJDIR)/$(SRC_TOOL)/ad_mod.o $(OBJDIR)/$(SRC_TOOL)/cheb_mod.o $(OBJDIR)/$(SRC_TOOL)/nag_compat_mod.o $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/miscellaneous.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/scalar_burning.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver.o
$(OBJDIR)/$(SRC_CORE)/starting_model.o: $(OBJDIR)/$(SRC_THEORY)/rotation_solver.o $(OBJDIR)/$(SRC_CORE)/miscellaneous.o $(OBJDIR)/$(SRC_CORE)/scalar_burning.o
$(OBJDIR)/$(SRC_CORE)/shoot_v2.o: $(OBJDIR)/$(SRC_CORE)/starting_model.o
$(OBJDIR)/$(SRC_CORE)/MRcurve.o: $(OBJDIR)/$(SRC_THEORY)/rotation_solver.o $(OBJDIR)/$(SRC_CORE)/starting_model.o
$(OBJDIR)/src/main.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(MODDIR)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) -r $(BUILDDIR)
