PROG    := a.out
.DEFAULT_GOAL := all
FC      = gfortran

MODE    ?= Release
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
FFTW_LIBS ?= $(shell pkg-config --libs fftw3 2>/dev/null || echo -lfftw3)
LIBS    ?= -llapack -lblas $(FFTW_LIBS)
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
  $(SRC_TOOL)/nag_compat_mod.f90 \
  $(SRC_TOOL)/brent_mod.f90 \
  $(SRC_TOOL)/cheb_mod.f90 \
  $(SRC_TOOL)/donutization_mod.f90 \
  $(SRC_TOOL)/spectral_hub_mod.f90 \
  $(SRC_CORE)/para_mod.f90 \
  $(SRC_CORE)/ope_eq_mod.f90 \
  $(SRC_CORE)/constraint_mod.f90 \
  $(SRC_TOOL)/toolkit_mod.f90 \
  $(SRC_CORE)/miscellaneous_mod.f90 \
  $(SRC_CORE)/set_disk.f90 \
  $(SRC_THEORY)/relaxation_mod.f90 \
  $(SRC_THEORY)/rotational_law_mod.f90 \
  $(SRC_THEORY)/spin_derivatives_mod.f90 \
  $(SRC_THEORY)/spin_updates_mod.f90 \
  $(SRC_THEORY)/spin_workspace_mod.f90 \
  $(SRC_THEORY)/spin_integration_mod.f90 \
  $(SRC_THEORY)/spin_relaxation_mod.f90 \
  $(SRC_THEORY)/rotation_solver_mod.f90 \
  $(SRC_CORE)/shoot_solver_2d_mod.f90 \
  $(SRC_CORE)/shoot_solver_1d_hc_mod.f90 \
  $(SRC_CORE)/shoot_solver_1d_r_ratio_mod.f90 \
  $(SRC_CORE)/eos_mod.f90 \
  $(SRC_TOOL)/exporter_mod.f90 \
  $(SRC_CORE)/grid_mod.f90 \
  $(SRC_CORE)/sphere_mod.f90 \
  $(SRC_CORE)/regrid_mod.f90 \
  $(SRC_CORE)/analysis_mod.f90 \
  $(SRC_CORE)/scalar_burning_mod.f90 \
  $(SRC_CORE)/starting_model_mod.f90 \
  $(SRC_CORE)/shoot_mod.f90 \
  $(SRC_CORE)/MRcurve_mod.f90 \
  src/main.f90

OBJECTS := $(patsubst %.f90,$(OBJDIR)/%.o,$(SOURCES))

# Ensure precision module is built before any source that imports it.
$(filter-out $(OBJDIR)/$(SRC_CORE)/precision_mod.o,$(OBJECTS)): $(OBJDIR)/$(SRC_CORE)/precision_mod.o
$(OBJDIR)/$(SRC_TOOL)/donutization_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o: $(OBJDIR)/$(SRC_TOOL)/ad_mod.o $(OBJDIR)/$(SRC_TOOL)/spectral_hub_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/ope_eq_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/constraint_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/ope_eq_mod.o
$(OBJDIR)/$(SRC_CORE)/miscellaneous_mod.o: $(OBJDIR)/$(SRC_TOOL)/nag_compat_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/set_disk.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/relaxation_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/rotational_law_mod.o: $(OBJDIR)/$(SRC_TOOL)/brent_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/spin_derivatives_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/spin_updates_mod.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotational_law_mod.o $(OBJDIR)/$(SRC_TOOL)/brent_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o
$(OBJDIR)/$(SRC_THEORY)/spin_workspace_mod.o: $(OBJDIR)/$(SRC_TOOL)/nag_compat_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_TOOL)/exporter_mod.o: $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_THEORY)/spin_integration_mod.o: $(OBJDIR)/$(SRC_THEORY)/spin_workspace_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_derivatives_mod.o $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_TOOL)/exporter_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o
$(OBJDIR)/$(SRC_THEORY)/spin_relaxation_mod.o: $(OBJDIR)/$(SRC_THEORY)/spin_workspace_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/relaxation_mod.o
$(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o: $(OBJDIR)/$(SRC_THEORY)/spin_workspace_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_integration_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_relaxation_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_updates_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o
$(OBJDIR)/$(SRC_CORE)/shoot_solver_2d_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o
$(OBJDIR)/$(SRC_CORE)/shoot_solver_1d_hc_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o
$(OBJDIR)/$(SRC_CORE)/shoot_solver_1d_r_ratio_mod.o: $(OBJDIR)/$(SRC_CORE)/shoot_solver_1d_hc_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o
$(OBJDIR)/$(SRC_CORE)/eos_mod.o: $(OBJDIR)/$(SRC_TOOL)/ad_mod.o $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/grid_mod.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_CORE)/sphere_mod.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/set_disk.o
$(OBJDIR)/$(SRC_CORE)/regrid_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/grid_mod.o
$(OBJDIR)/$(SRC_CORE)/analysis_mod.o: $(OBJDIR)/$(SRC_TOOL)/ad_mod.o $(OBJDIR)/$(SRC_TOOL)/cheb_mod.o $(OBJDIR)/$(SRC_TOOL)/nag_compat_mod.o $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_TOOL)/exporter_mod.o $(OBJDIR)/$(SRC_CORE)/miscellaneous_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o
$(OBJDIR)/$(SRC_CORE)/scalar_burning_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o
$(OBJDIR)/$(SRC_CORE)/starting_model_mod.o: $(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o $(OBJDIR)/$(SRC_CORE)/miscellaneous_mod.o $(OBJDIR)/$(SRC_CORE)/scalar_burning_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o $(OBJDIR)/$(SRC_CORE)/regrid_mod.o $(OBJDIR)/$(SRC_CORE)/sphere_mod.o
$(OBJDIR)/$(SRC_CORE)/shoot_mod.o: $(OBJDIR)/$(SRC_TOOL)/donutization_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_workspace_mod.o $(OBJDIR)/$(SRC_CORE)/starting_model_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o $(OBJDIR)/$(SRC_CORE)/shoot_solver_2d_mod.o $(OBJDIR)/$(SRC_CORE)/shoot_solver_1d_hc_mod.o $(OBJDIR)/$(SRC_CORE)/shoot_solver_1d_r_ratio_mod.o
$(OBJDIR)/$(SRC_CORE)/MRcurve_mod.o: $(OBJDIR)/$(SRC_TOOL)/donutization_mod.o $(OBJDIR)/$(SRC_THEORY)/spin_workspace_mod.o $(OBJDIR)/$(SRC_THEORY)/rotation_solver_mod.o $(OBJDIR)/$(SRC_CORE)/starting_model_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/analysis_mod.o
$(OBJDIR)/$(SRC_CORE)/set_disk.o: $(OBJDIR)/$(SRC_CORE)/eos_mod.o
$(OBJDIR)/src/main.o: $(OBJDIR)/$(SRC_TOOL)/toolkit_mod.o $(OBJDIR)/$(SRC_CORE)/para_mod.o $(OBJDIR)/$(SRC_CORE)/constraint_mod.o $(OBJDIR)/$(SRC_CORE)/eos_mod.o $(OBJDIR)/$(SRC_CORE)/grid_mod.o $(OBJDIR)/$(SRC_CORE)/shoot_mod.o

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

# --- Test targets ---
LIBGRASS := libgrass.a
LIB_OBJECTS := $(filter-out $(OBJDIR)/src/main.o,$(OBJECTS))

$(LIBGRASS): $(LIB_OBJECTS)
	ar rcs $@ $^

.PHONY: test test-unit test-integration test-debug

test: test-unit test-integration
	@echo "All tests passed."

ABSFFLAGS := -I$(CURDIR) -I$(CURDIR)/$(MODDIR) -J$(CURDIR)/$(MODDIR) \
  $(filter-out -I. -I$(MODDIR) -J$(MODDIR),$(FFLAGS))

test-unit: $(LIBGRASS)
	@$(MAKE) -C tests unit FC=$(FC) FFLAGS="$(ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(LIBGRASS) LIBS="$(LIBS)"

test-integration: $(LIBGRASS)
	@$(MAKE) -C tests integration FC=$(FC) FFLAGS="$(ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(LIBGRASS) LIBS="$(LIBS)"

test-debug:
	@$(MAKE) MODE=Debug $(LIBGRASS)
	@$(MAKE) -C tests sanitizer FC=$(FC) FFLAGS="$(ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(LIBGRASS) LIBS="$(LIBS)"
