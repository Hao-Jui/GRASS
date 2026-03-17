PROG    := a.out
.DEFAULT_GOAL := all
FC      = gfortran

MODE    ?= Debug
BASE_FFLAGS := -I.
RELEASE_FFLAGS := -O3 -march=native
DEBUG_FFLAGS := -O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -fcheck=all -Wuninitialized -Wconversion -Wuse-without-only -finit-real=nan
ifeq ($(MODE),Release)
  MODE_FFLAGS := $(RELEASE_FFLAGS)
else ifeq ($(MODE),Debug)
  MODE_FFLAGS := $(DEBUG_FFLAGS)
else
  $(error Unsupported MODE='$(MODE)'. Use MODE=Release or MODE=Debug)
endif
FFLAGS  ?= $(BASE_FFLAGS) $(MODE_FFLAGS)

LDFLAGS ?= -Wl,-stack_size,0x4000000
LIBS    ?= -llapack -lblas
BUILDDIR := build
OBJDIR   := $(BUILDDIR)/obj
MODDIR   := $(BUILDDIR)/mod
BINDIR   := $(BUILDDIR)/bin
TARGET   := $(BINDIR)/$(PROG)

FFLAGS  += -J$(MODDIR) -I$(MODDIR)
SRCDIR := src
SRC_TOOL := $(SRCDIR)/tool
SRC_CORE := $(SRCDIR)/core
SRC_THEORY := $(SRCDIR)/theory

SOURCES := \
  $(SRC_CORE)/precision_mod.f90 \
  $(SRC_TOOL)/ad_mod.f90 \
  $(SRC_TOOL)/simpson_mod.f90 \
  $(SRC_TOOL)/nag_compat_mod.f90 \
  $(SRC_TOOL)/brent.f90 \
  $(SRC_TOOL)/cheb_mod.f90 \
  $(SRC_CORE)/para_mod.f90 \
  $(SRC_CORE)/Ope_eq.f90 \
  $(SRC_CORE)/constraint_eq.f90 \
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
$(OBJDIR)/$(SRC_CORE)/Ope_eq.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/constraint_eq.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/Ope_eq.o
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
$(OBJDIR)/src/main.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/constraint_eq.o

.PHONY: all clean release debug

all: $(TARGET)

release:
  $(MAKE) MODE=Release

debug:
  $(MAKE) MODE=Debug

$(TARGET): $(OBJECTS)
	@mkdir -p $(BINDIR)
	$(FC) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

$(OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(MODDIR)
	$(FC) $(FFLAGS) -c -o $@ $<

clean:
	$(RM) -r $(BUILDDIR)
