PROG    := a.out
.DEFAULT_GOAL := all
FC      = gfortran

MODE    ?= Release
BASE_FFLAGS := -I. -ffree-line-length-none
RELEASE_FFLAGS := -O3 -march=native
DEBUG_FFLAGS := -O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -fcheck=all -Wuninitialized -Wconversion -Wuse-without-only -finit-real=nan
SANITIZER_FFLAGS := -O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface -fcheck=all -finit-real=nan -fsanitize=address,undefined
ifeq ($(MODE),Release)
  MODE_FFLAGS := $(RELEASE_FFLAGS)
else ifeq ($(MODE),Debug)
  MODE_FFLAGS := $(DEBUG_FFLAGS)
else ifeq ($(MODE),Sanitizer)
  MODE_FFLAGS := $(SANITIZER_FFLAGS)
else
  $(error Unsupported MODE='$(MODE)'. Use MODE=Release, MODE=Debug, or MODE=Sanitizer)
endif
FFLAGS  ?= $(BASE_FFLAGS) $(MODE_FFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  LDFLAGS ?= -Wl,-stack_size,0x4000000
  LAPACK_BLAS_LIBS ?= -framework Accelerate
else
  LDFLAGS ?=
  LAPACK_BLAS_LIBS ?= -llapack -lblas
endif
FFTW_LIBS ?= $(shell pkg-config --libs fftw3 2>/dev/null || echo -lfftw3)
LIBS    ?= $(LAPACK_BLAS_LIBS) $(FFTW_LIBS)
BUILDDIR := build
OBJDIR   := $(BUILDDIR)/obj
MODDIR   := $(BUILDDIR)/mod
BINDIR   := $(BUILDDIR)/bin
TARGET   := $(BINDIR)/$(PROG)
LIBGRASS := bin/libgrass.a
CI_LIBGRASS := bin/libgrass_ci.a

FFLAGS  += -J$(MODDIR) -I$(MODDIR)
SRCDIR := src
SRC_TOOL := $(SRCDIR)/tool
SRC_CORE := $(SRCDIR)/core
SRC_THEORY := $(SRCDIR)/theory
SRC_BH_TOROID := $(SRCDIR)/bh_toroid
SRC_UTIL := $(SRCDIR)/utilities

# SOURCES is topologically sorted so sequential make works without
# explicit dependency rules.  For parallel builds use CMake instead.
SOURCES := \
  $(SRC_CORE)/precision_mod.f90 \
  $(SRC_TOOL)/ad_mod.f90 \
  $(SRC_TOOL)/nag_compat_mod.f90 \
  $(SRC_TOOL)/lapack_interfaces_mod.f90 \
  $(SRC_TOOL)/brent_mod.f90 \
  $(SRC_TOOL)/cheb_mod.f90 \
  $(SRC_TOOL)/spline_mod.f90 \
  $(SRC_TOOL)/spectral_hub_mod.f90 \
  src/para_panel.f90 \
  $(SRC_BH_TOROID)/core_mod.f90 \
  $(SRC_BH_TOROID)/updates_mod.f90 \
  $(SRC_BH_TOROID)/solver_mod.f90 \
  $(SRC_TOOL)/donutization_mod.f90 \
  $(SRC_TOOL)/toolkit_mod.f90 \
  $(SRC_CORE)/ope_eq_mod.f90 \
  $(SRC_CORE)/constraint_mod.f90 \
  $(SRC_CORE)/eos_mod.f90 \
  $(SRC_CORE)/miscellaneous_mod.f90 \
  $(SRC_CORE)/set_disk.f90 \
  $(SRC_TOOL)/exporter_mod.f90 \
  $(SRC_CORE)/grid_mod.f90 \
  $(SRC_THEORY)/relaxation_mod.f90 \
  $(SRC_THEORY)/rotational_law_mod.f90 \
  $(SRC_THEORY)/spin_derivatives_mod.f90 \
  $(SRC_THEORY)/spin_workspace_mod.f90 \
  $(SRC_THEORY)/spin_updates_mod.f90 \
  $(SRC_CORE)/sphere_mod.f90 \
  $(SRC_CORE)/restart_format_mod.f90 \
  $(SRC_CORE)/regrid_mod.f90 \
  $(SRC_CORE)/analysis_mod.f90 \
  $(SRC_THEORY)/spin_integration_mod.f90 \
  $(SRC_THEORY)/spin_relaxation_mod.f90 \
  $(SRC_THEORY)/rotation_solver_mod.f90 \
  $(SRC_UTIL)/scalar_burning_mod.f90 \
  $(SRC_CORE)/shoot_solver_2d_mod.f90 \
  $(SRC_CORE)/shoot_solver_1d_hc_mod.f90 \
  $(SRC_CORE)/shoot_solver_1d_r_ratio_mod.f90 \
  $(SRC_CORE)/starting_model_mod.f90 \
  $(SRC_UTIL)/shoot_mod.f90 \
  $(SRC_UTIL)/MRcurve_mod.f90 \
  src/main.f90

OBJECTS := $(patsubst %.f90,$(OBJDIR)/%.o,$(SOURCES))
CI_PARA_FILE := tests/CI/para_panel.f90
CI_OBJDIR := $(BUILDDIR)/make_ci_obj
CI_MODDIR := $(BUILDDIR)/make_ci_mod
CI_SOURCES := $(patsubst src/para_panel.f90,$(CI_PARA_FILE),$(filter-out src/main.f90,$(SOURCES)))
CI_OBJECTS := $(patsubst %.f90,$(CI_OBJDIR)/%.o,$(CI_SOURCES))
CI_FFLAGS := $(filter-out -I$(MODDIR) -J$(MODDIR),$(FFLAGS)) -I$(CI_MODDIR) -J$(CI_MODDIR)


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

$(CI_OBJDIR)/%.o: %.f90
	@mkdir -p $(dir $@) $(CI_MODDIR)
	$(FC) $(CI_FFLAGS) -c -o $@ $<

clean:
	$(RM) -r $(BUILDDIR) $(LIBGRASS) $(CI_LIBGRASS) tests/bin

# --- Test targets ---
LIB_OBJECTS := $(filter-out $(OBJDIR)/src/main.o,$(OBJECTS))

$(LIBGRASS): $(LIB_OBJECTS)
	@mkdir -p $(dir $@)
	$(RM) $@
	ar rcs $@ $^

$(CI_LIBGRASS): $(CI_OBJECTS)
	@mkdir -p $(dir $@)
	$(RM) $@
	ar rcs $@ $^

.PHONY: test test-unit test-integration test-debug test-sanitizer

test: test-unit test-integration
	@echo "All tests passed."

ABSFFLAGS := -I$(CURDIR) -I$(CURDIR)/$(MODDIR) -J$(CURDIR)/$(MODDIR) \
  $(filter-out -I. -I$(MODDIR) -J$(MODDIR),$(FFLAGS))
CI_ABSFFLAGS := -I$(CURDIR) -I$(CURDIR)/$(CI_MODDIR) -J$(CURDIR)/$(CI_MODDIR) \
  $(filter-out -I. -I$(MODDIR) -J$(MODDIR),$(FFLAGS))

test-unit: $(CI_LIBGRASS)
	@$(MAKE) -C tests unit FC=$(FC) FFLAGS="$(CI_ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(CI_LIBGRASS) LIBS="$(LIBS)"

test-integration: $(CI_LIBGRASS)
	@$(MAKE) -C tests integration FC=$(FC) FFLAGS="$(CI_ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(CI_LIBGRASS) LIBS="$(LIBS)"

test-debug:
	@$(MAKE) MODE=Debug $(CI_LIBGRASS)
	@$(MAKE) -C tests sanitizer FC=$(FC) FFLAGS="$(CI_ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(CI_LIBGRASS) LIBS="$(LIBS)"

test-sanitizer:
	@$(MAKE) MODE=Sanitizer $(CI_LIBGRASS)
	@$(MAKE) -C tests unit FC=$(FC) FFLAGS="$(CI_ABSFFLAGS)" LIBGRASS=$(CURDIR)/$(CI_LIBGRASS) LIBS="$(LIBS) -fsanitize=address,undefined"
