FC ?= gfortran
FFLAGS ?=
LIBGRASS ?= ../bin/libgrass_ci.a
LIBS ?= -llapack -lblas -lfftw3
BINDIR := bin

UNIT_TESTS := test_ad test_spline test_spectral test_brent test_eos \
  test_bh_toroid_params test_bh_toroid_radial_map test_bh_toroid_green \
  test_bh_toroid_updates test_bh_toroid_solver
INTEG_TESTS := test_gr_uniform test_gr_constj test_gr_uryu test_st_uniform test_st_uniform_r07
RESTART_TEST := test_restart

.PHONY: unit integration sanitizer
.NOTPARALLEL: integration

unit: $(UNIT_TESTS)

integration: $(INTEG_TESTS) $(RESTART_TEST)

sanitizer: test_gr_uniform
	@echo "Sanitizer build passed."

$(filter-out test_spline,$(UNIT_TESTS)): %: unit/%.f90 unit/test_utils.f90 $(LIBGRASS)
	@mkdir -p $(BINDIR)
	@$(FC) $(FFLAGS) -o $(BINDIR)/$@ unit/test_utils.f90 unit/$*.f90 $(LIBGRASS) $(LIBS)
	@printf "  $@: " && cd .. && timeout 10 tests/$(BINDIR)/$@

test_spline: unit/test_spline.f90 unit/test_utils.f90 $(LIBGRASS)
	@mkdir -p $(BINDIR)
	@$(FC) $(FFLAGS) -o $(BINDIR)/$@ unit/test_utils.f90 unit/test_spline.f90 $(LIBGRASS) $(LIBS)
	@printf "  test_spline: " && cd .. && timeout 10 tests/$(BINDIR)/test_spline

$(filter-out test_gr_uryu,$(INTEG_TESTS)): %: integration/%.f90 $(LIBGRASS)
	@mkdir -p $(BINDIR)
	@$(FC) $(FFLAGS) -o $(BINDIR)/$@ $< $(LIBGRASS) $(LIBS)
	@cd .. && mkdir -p Cont Res && GRASS_TIMING=0 GRASS_E_CENTER=0.8e15 timeout 120 tests/$(BINDIR)/$@ > /dev/null 2>&1 \
	  && bash tests/check_regression.sh Cont/properties.dat tests/reference/$@.ref \
	  || (printf "  $@: FAIL\n" && exit 1)

test_gr_uryu: integration/test_gr_uryu.f90 $(LIBGRASS)
	@mkdir -p $(BINDIR)
	@$(FC) $(FFLAGS) -o $(BINDIR)/$@ $< $(LIBGRASS) $(LIBS)
	@cd .. && mkdir -p Cont Res && GRASS_TIMING=0 timeout 120 \
	  bash tests/run_regression.sh tests/$(BINDIR)/$@ Cont/properties.dat tests/reference/$@.ref \
	  || (printf "  $@: FAIL\n" && exit 1)

# Restart test: runs after all integration tests, reads whatever Res/res.rst exists.
# Only checks convergence (exit code 0), not reference values — the restart source varies.
$(RESTART_TEST): $(INTEG_TESTS) integration/test_restart.f90 $(LIBGRASS)
	@mkdir -p $(BINDIR)
	@$(FC) $(FFLAGS) -o $(BINDIR)/test_restart integration/test_restart.f90 $(LIBGRASS) $(LIBS)
	@printf "  test_restart: " && cd .. && mkdir -p Cont Res \
	  && GRASS_TIMING=0 GRASS_E_CENTER=0.8e15 timeout 120 tests/$(BINDIR)/test_restart > /dev/null 2>&1 \
	  && echo "PASS (convergence)" || (echo "FAIL" && exit 1)
