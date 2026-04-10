# Regression Test Suite — Industrial Standard

**Date:** 2026-04-10
**Branch:** `rebuild` (d2acb03 + donutization/exporter)
**Goal:** `make test` catches convergence regressions, numerical bugs, and UB — all within ~5 seconds.

---

## Lessons from the Uryu Regression

The uryu convergence bug went undetected for 5 commits because:
1. No test exercised the uryu code path
2. No test verified all 3 solver_types converge
3. No debug-mode test caught the `.and. size(unallocated)` UB
4. No EOS consistency check existed

This suite prevents all four.

---

## Test Architecture

```
make test
├── Unit tests (standalone programs, ~1s)
│   ├── test_ad          — dual-number operator correctness
│   ├── test_spline      — cubic reproduction, knot interpolation
│   ├── test_spectral    — GL/Lobatto weights, quadrature exactness
│   ├── test_brent       — root-finding convergence
│   └── test_eos         — EOS load, eval, pe_at_h consistency
│
├── Integration tests (full solver, ~3s total)
│   ├── test_gr_uniform  — GR + uniform, r_ratio=0.7 (0.03s)
│   ├── test_gr_constj   — GR + const_j, r_ratio=0.7 (0.09s)
│   ├── test_gr_uryu     — GR + uryu,    r_ratio=0.7 (0.30s)
│   └── test_st_uniform  — ST + uniform, r_ratio=1.0 (0.03s)
│
└── Sanitizer build (~2s build, instant check)
    └── make test-debug  — build with -fcheck=all, run ST+uniform
                           catches UB, uninitialized, bounds errors
```

**Total `make test` time: ~5 seconds.**

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|--------------|
| AC1 | `make test` runs all tests and reports PASS/FAIL | Exit code 0 = all pass |
| AC2 | GR+uryu test converges with axial ratio = 0.700 | grep output |
| AC3 | GR+uniform, GR+const_j, ST+uniform all converge | grep "Converged" |
| AC4 | All unit tests pass | Each test program exits 0 |
| AC5 | `make test-debug` builds with `-fcheck=all` and runs without runtime errors | Exit code 0 |
| AC6 | Introducing a known regression (e.g., changing `+` to `-` in ad_mod) is caught | Test fails |
| AC7 | Total wall time < 10 seconds | Timed |

---

## Implementation Plan

### Tier 1: Integration Tests (highest value — catches regressions like the uryu bug)

**Approach:** 4 test programs, each a minimal `main.f90` that sets config, calls the solver, and checks convergence. No pFUnit needed — just standalone Fortran programs.

**Directory structure:**
```
tests/
├── Makefile              # Test runner
├── integration/
│   ├── test_gr_uniform.f90
│   ├── test_gr_constj.f90
│   ├── test_gr_uryu.f90
│   └── test_st_uniform.f90
└── unit/
    ├── test_ad.f90
    ├── test_spline.f90
    ├── test_spectral.f90
    ├── test_brent.f90
    └── test_eos.f90
```

**Each integration test program (~30 lines):**
```fortran
program test_gr_uryu
  use para_mod
  use eos_mod, only: loadEos
  use grid_mod, only: make_grid, GridTrig
  use starting_model_mod, only: initialize_starting_model
  implicit none

  ! Override config for this test
  active_theory = THEORY_GR
  solver_type = "uryu"
  run_mode = MODE_DEFAULT
  run_task = OneModel
  eos_file = "MPA1"
  SDIV = 401; MDIV = 41
  r_ratio = 0.7_wp
  e_center_init = 0.7e15_wp   ! (set before sphere)

  call initialize_theory()
  call loadEos()
  call make_grid(); call GridTrig()
  call initialize_starting_model()
  ! If we reach here, the solver converged (it calls stop on success)
end program
```

**Key advantage:** Each test sets config at RUNTIME by overriding `para_mod` variables — no need to copy/modify source files or recompile the library. All 4 tests link against the SAME compiled library.

**Build approach:**
1. Compile all `src/` into a static library `libgrass.a`
2. Each test links against `libgrass.a`
3. Only the test `main.f90` differs — no source duplication

**Caveat:** `r_ratio` and `e_center` are set inside `initialize_starting_model`, not as config variables. Two approaches:
- **Option A:** Add thin config variables (`r_ratio_init`, `e_center_init`) to `para_mod` that `starting_model` reads instead of hardcoding
- **Option B:** Each test provides its own `starting_model` wrapper that overrides after the `merge()` call

Option A is cleaner and was partially done in the current code (the `merge()` pattern). I recommend adding `e_center_init` and `r_ratio_init` to `para_mod` so tests can set them without source modification.

### Tier 2: Unit Tests (catches numerical bugs in pure modules)

**Approach:** Standalone test programs (no pFUnit). Each program tests one module, uses simple assertion macros via a `test_utils.f90` helper.

**test_utils.f90** (~20 lines):
```fortran
module test_utils
  use precision_mod, only: wp
  implicit none
  integer :: n_pass = 0, n_fail = 0
contains
  subroutine assert_near(label, expected, actual, tol)
    character(*), intent(in) :: label
    real(wp), intent(in) :: expected, actual, tol
    if (abs(expected - actual) > tol) then
      write(*,'(A,A,2es15.7)') "FAIL: ", label, expected, actual
      n_fail = n_fail + 1
    else
      n_pass = n_pass + 1
    end if
  end subroutine
  subroutine test_summary()
    write(*,'(I0,A,I0,A)') n_pass, " passed, ", n_fail, " failed"
    if (n_fail > 0) error stop "Tests failed"
  end subroutine
end module
```

**Test programs:**
- `test_ad.f90` — 12 assertions (same as the pFUnit tests we wrote earlier)
- `test_spline.f90` — cubic reproduction, linear data, knot interpolation
- `test_spectral.f90` — GL weight sum, exactness, Lobatto endpoints
- `test_brent.f90` — root finding with known-solution callbacks
- `test_eos.f90` — load MPA1, eval at midpoint, check finite + positive

### Tier 3: Debug Sanitizer Build

**`make test-debug`** — builds the entire codebase with `-O0 -fcheck=all -finit-real=nan`, runs ST+uniform. Catches:
- Uninitialized variables (NaN propagation)
- Array bounds violations
- `.and.`/`.or.` with unallocated arrays
- Integer overflow

This is the cheapest way to catch UB — the same class of bug that caused the uryu regression.

### Tier 4: Makefile Integration

**Top-level Makefile additions:**

```makefile
.PHONY: test test-debug

TESTDIR := tests
TESTBIN := $(TESTDIR)/bin

# Build library from all sources (except main.f90)
libgrass.a: $(filter-out $(OBJDIR)/src/main.o,$(OBJECTS))
	ar rcs $@ $^

# Unit tests
test-unit: libgrass.a
	@$(MAKE) -C $(TESTDIR) unit LIBGRASS=$(CURDIR)/libgrass.a MODDIR=$(CURDIR)/$(MODDIR)

# Integration tests
test-integration: libgrass.a
	@$(MAKE) -C $(TESTDIR) integration LIBGRASS=$(CURDIR)/libgrass.a MODDIR=$(CURDIR)/$(MODDIR)

# All tests
test: test-unit test-integration
	@echo "All tests passed."

# Debug sanitizer
test-debug:
	@$(MAKE) MODE=Debug libgrass.a
	@$(MAKE) -C $(TESTDIR) sanitizer LIBGRASS=$(CURDIR)/libgrass.a MODDIR=$(CURDIR)/$(MODDIR)
```

**tests/Makefile:**
```makefile
FC = gfortran
FFLAGS = -I$(MODDIR)
LIBS = $(LIBGRASS) -llapack -lblas $(shell pkg-config --libs fftw3)

unit: test_ad test_spline test_spectral test_brent test_eos
integration: test_gr_uniform test_gr_constj test_gr_uryu test_st_uniform

test_%: unit/test_%.f90 unit/test_utils.f90
	$(FC) $(FFLAGS) -o $(TESTBIN)/$@ $^ $(LIBS) && $(TESTBIN)/$@

test_gr_% test_st_%: integration/test_$*.f90
	$(FC) $(FFLAGS) -o $(TESTBIN)/$@ $< $(LIBS) && timeout 60 $(TESTBIN)/$@
```

---

## Reference Values for Integration Tests

Each integration test checks against reference values from the known-good `d2acb03` codebase:

| Test | Axial Ratio | Mass (M_o) | Converge? | Max Time |
|------|------------|------------|-----------|----------|
| GR+uniform | 0.700 | ~1.30 | Yes | 1s |
| GR+const_j | 0.700 | ~1.30 | Yes | 1s |
| GR+uryu | 0.700 | ~1.30 | Yes | 2s |
| ST+uniform | 1.000 | ~1.8 | Yes | 1s |

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Test programs can't override `r_ratio`/`e_center` (set in starting_model) | Add `r_ratio_init`/`e_center_init` config vars to `para_mod` |
| Library approach requires relinking on any source change | `make test` depends on `libgrass.a` which depends on all objects |
| Integration tests write to `Cont/` (file collisions) | Tests run from `tests/` working directory; symlink `eos/` |
| Debug build is slow at high resolution | Use res=50 for sanitizer test (only needs to start iterating, not converge) |

---

## Execution Order

1. **Add `test_utils.f90`** — assertion helper (10 min)
2. **Add unit tests** — port from the pFUnit tests we already wrote (30 min)
3. **Add `r_ratio_init` / `e_center_init`** to `para_mod` + `starting_model` (10 min)
4. **Add integration tests** — 4 test programs (20 min)
5. **Add `tests/Makefile`** + top-level `test` target (15 min)
6. **Verify `make test`** — all pass, <10s total (5 min)
7. **Add `make test-debug`** — sanitizer build (10 min)
8. **Regression check** — introduce a deliberate bug, verify test catches it (5 min)

---

*Plan generated: 2026-04-10. Informed by the uryu regression investigation earlier in this session.*
