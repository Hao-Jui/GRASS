# GRASS Testing Framework Plan — pFUnit + CMake

**Date:** 2026-04-10
**Goal:** Bring Testing from 2/10 to 7/10 by introducing pFUnit with CMake alongside the existing Makefile.
**Scope:** Unit tests for `numerics/` (6 modules) and `eos/` (1 module), plus build integration.

---

## Requirements Summary

1. Install pFUnit and integrate with a CMake test build
2. Keep the production Makefile untouched (add only a `test` convenience target)
3. Write unit tests for all stateless numerics modules (Tier 1-2) and an integration test for EOS (Tier 3)
4. Provide a single-command test workflow: `make test`

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|--------------|
| AC1 | `cmake --build build-test` compiles all test executables without errors | Build log shows 0 errors |
| AC2 | `ctest --test-dir build-test` runs all tests, 100% pass | ctest output: "100% tests passed" |
| AC3 | `make test` from project root builds and runs all tests | Single-command workflow |
| AC4 | ad_mod tests verify all 6 operator categories + exp/log derivatives | >= 12 assertions pass |
| AC5 | spline_mod tests verify cubic polynomial reproduction to machine epsilon | Residual < 1e-14 |
| AC6 | spectral_hub_mod tests verify quadrature exactness for polynomials of correct degree | Gauss-Legendre n-point rule exact for degree 2n-1 |
| AC7 | brent_mod tests verify root-finding convergence to specified tolerance | |f(root)| < tol |
| AC8 | cheb_mod tests verify Chebyshev fit reproduces polynomials of degree <= n | Residual < 1e-12 |
| AC9 | eos_mod integration test loads a fixture table and evaluates without crash | eos_eval returns finite value matching reference |
| AC10 | Production Makefile (`make`, `make debug`) still works unchanged | Build succeeds, binary runs |

---

## Architecture Decision

**Strategy:** CMake for test build only; production Makefile stays as-is.

CMake compiles the source modules needed by each test suite and links against pFUnit.
This avoids: (a) converting the entire build to CMake, (b) fragile cross-build-system linking,
(c) `.mod` file version mismatches between separate gfortran invocations.

**Directory layout after implementation:**
```
GRASS/
├── CMakeLists.txt              # NEW — top-level CMake (test build only)
├── Makefile                    # EXISTING — add `test` target (3 lines)
├── tests/
│   ├── CMakeLists.txt          # NEW — test registration
│   ├── test_ad_mod.pf          # NEW — Tier 1
│   ├── test_toolkit_helpers.pf # NEW — Tier 1
│   ├── test_spline_mod.pf      # NEW — Tier 2
│   ├── test_cheb_mod.pf        # NEW — Tier 2
│   ├── test_spectral_hub.pf    # NEW — Tier 2
│   ├── test_brent_mod.pf       # NEW — Tier 2
│   └── test_eos_mod.pf         # NEW — Tier 3
└── tests/fixtures/
    └── test_poly.dat           # NEW — minimal EOS fixture (synthetic)
```

---

## Implementation Steps

### Phase 1: Infrastructure (files: CMakeLists.txt, tests/CMakeLists.txt, Makefile)

**Step 1.1 — Install pFUnit**
- Install via Homebrew: `brew install pfunit` (if available), otherwise:
  ```
  git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
  cd pFUnit && mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local
  make -j4 && make install
  ```
- Verify: `find $HOME/.local -name "pFUnitConfig.cmake"` returns a path

**Step 1.2 — Create top-level CMakeLists.txt**
- `CMakeLists.txt` at project root
- Minimum CMake version: 3.12 (for pFUnit 4.x compatibility)
- Project language: Fortran
- Compile source modules needed by tests (not main.f90 — tests provide their own drivers)
- `find_package(PFUNIT REQUIRED)`
- `add_subdirectory(tests)`
- Compiler flags: match debug flags from Makefile line 8

**Step 1.3 — Create tests/CMakeLists.txt**
- Register each `.pf` file as a pFUnit test via `add_pfunit_ctest()`
- Link each test against the compiled source modules it needs
- Link external libs: LAPACK, BLAS, FFTW3

**Step 1.4 — Add `test` target to Makefile**
- Append to `Makefile` (after `clean:` block):
  ```makefile
  .PHONY: test
  test:
  	@cmake -S . -B build-test -DCMAKE_BUILD_TYPE=Debug 2>&1 | tail -1
  	@cmake --build build-test 2>&1 | tail -5
  	@ctest --test-dir build-test --output-on-failure
  ```

### Phase 2: Tier 1 Tests — Pure, Zero-Dependency (files: tests/test_ad_mod.pf, tests/test_toolkit_helpers.pf)

**Step 2.1 — test_ad_mod.pf**
Module under test: `src/numerics/ad_mod.f90`
Dependencies: none (self-contained, uses `real(8)` directly)

| Test | What it verifies |
|------|-----------------|
| `test_dual_const` | `dual_const(3.0)` gives `val=3.0, der=0.0` |
| `test_dual_var` | `dual_var(3.0)` gives `val=3.0, der=1.0` |
| `test_add_dd` | `(a,a') + (b,b') = (a+b, a'+b')` |
| `test_add_dr` | `dual + real` preserves derivative |
| `test_sub_dd` | `(a,a') - (b,b') = (a-b, a'-b')` |
| `test_neg_d` | `-(a,a') = (-a, -a')` |
| `test_mul_dd` | `(a,a')*(b,b') = (ab, a'b+ab')` — product rule |
| `test_mul_dr` | `dual * real` scales derivative |
| `test_div_dd` | `(a,a')/(b,b') = (a/b, (a'b-ab')/b^2)` — quotient rule |
| `test_exp_d` | `exp(dual_var(x)).der = exp(x)` |
| `test_log_d` | `log(dual_var(x)).der = 1/x` |
| `test_chain_rule` | `exp(2*dual_var(1)).der = 2*exp(2)` — composite |

Tolerance: `1e-14` (machine epsilon for real64 is ~2.2e-16)

**Step 2.2 — test_toolkit_helpers.pf**
Module under test: `src/numerics/toolkit_mod.f90` (pure subset only)
Dependencies: `precision_mod`

| Test | What it verifies |
|------|-----------------|
| `test_same_abscissa_equal` | `same_abscissa(1.0, 1.0)` returns `.true.` |
| `test_same_abscissa_diff` | `same_abscissa(1.0, 2.0)` returns `.false.` |
| `test_same_abscissa_near` | Values within relative epsilon return `.true.` |
| `test_binary_search` | Finds correct index in sorted array |
| `test_binary_search_left` | Value below range returns index 1 |
| `test_binary_search_right` | Value above range returns index n-1 |
| `test_pow_int_real` | `x^0=1`, `x^1=x`, `x^(-1)=1/x`, `x^3=x*x*x` |
| `test_interp_polynomial` | Barycentric interpolation reproduces cubic on 5 points |

### Phase 3: Tier 2 Tests — Core Numerics (4 test files)

**Step 3.1 — test_spline_mod.pf**
Module under test: `src/numerics/spline_mod.f90`
Dependencies: `precision_mod`

| Test | What it verifies |
|------|-----------------|
| `test_cubic_reproduction` | Spline of `y = x^3` on 10 uniform knots reproduces exactly (residual < 1e-14) |
| `test_interpolation_at_knots` | Evaluating spline at knot points returns input values |
| `test_linear_data` | Spline of `y = 2x + 1` gives `b(i) = 2`, `c(i) = d(i) = 0` |
| `test_monotonicity` | Spline of monotone data doesn't introduce spurious oscillation (spot check midpoints) |

Note: `build_spline_segmented` tested separately with a piecewise dataset containing one phase transition index.

**Step 3.2 — test_cheb_mod.pf**
Module under test: `src/numerics/cheb_mod.f90`
Dependencies: `precision_mod`, `spline_mod` (for `spline_coeff` type)

| Test | What it verifies |
|------|-----------------|
| `test_cheb_fit_polynomial` | `cheb_std_base` on `y = x^2` with degree >= 2 gives residual < 1e-12 |
| `test_clenshaw_eval` | `cheb_get_val_point` on known coefficients matches direct evaluation |
| `test_cheb_diff_coeffs` | Derivative coefficients of `T_n` match analytic formula |
| `test_diff_matrix_identity` | `cheb_diff_matrix(N)` applied to `x^k` values gives `k*x^(k-1)` |

External deps: LAPACK/BLAS (for `dgels` in `cheb_std_base`)

**Step 3.3 — test_spectral_hub.pf**
Module under test: `src/numerics/spectral_hub_mod.f90`
Dependencies: `precision_mod`

| Test | What it verifies |
|------|-----------------|
| `test_gl_weights_sum` | Gauss-Legendre weights sum to 1.0 (interval [0,1]) |
| `test_gl_exactness` | n-point GL integrates `x^(2n-1)` exactly (compare to `1/(2n)`) |
| `test_lobatto_endpoints` | `gauss_lobatto(n)` has `x(1)=0`, `x(n)=1` |
| `test_lobatto_weights_sum` | Gauss-Lobatto weights sum to 1.0 |
| `test_cheby_lobatto_formula` | `chebyshev_lobatto_points(5)` matches `(1 - cos(k*pi/4))/2` |
| `test_cc_weights_sum` | Clenshaw-Curtis weights sum to 1.0 |
| `test_integrate_tabulated` | Integral of `f(x) = x^2` on [0,1] returns 1/3 within tolerance |

**Step 3.4 — test_brent_mod.pf**
Module under test: `src/numerics/brent_mod.f90`
Dependencies: `precision_mod`

| Test | What it verifies |
|------|-----------------|
| `test_sqrt2` | Root of `f(x) = x^2 - 2` near x=1.5 gives `|root - sqrt(2)| < tol` |
| `test_sin_root` | Root of `sin(x)` near x=3 gives `|root - pi| < tol` |
| `test_convergence_flag` | `ierr = 0` on success |
| `test_bad_bracket` | `ierr /= 0` when initial bracket doesn't contain root |

Note: Tests use the `brent_func` abstract interface — test module provides simple callback subroutines.

### Phase 4: Tier 3 — EOS Integration Test (files: tests/test_eos_mod.pf, tests/fixtures/test_poly.dat)

**Step 4.1 — Create synthetic EOS fixture**
- `tests/fixtures/test_poly.dat` — minimal polytropic EOS table (20 rows)
- Format matches `eos/*.dat`: columns for log(e), log(p), log(h), log(n0)
- Use polytropic relation `P = K * e^Gamma` so we have analytic reference values
- This avoids depending on real EOS tables for unit tests

**Step 4.2 — test_eos_mod.pf**
Module under test: `src/eos/eos_mod.f90`
Dependencies: `precision_mod`, `ad_mod`, `toolkit_mod`, `spline_mod`, `para_panel` (fields)

| Test | What it verifies |
|------|-----------------|
| `test_load_fixture` | `loadEos()` with fixture file completes without crash |
| `test_eval_midpoint` | `eos_eval` at table midpoint returns value within 1% of analytic reference |
| `test_eval_dual_derivative` | `eos_eval_dual` derivative matches finite-difference `(f(x+h)-f(x-h))/(2h)` |
| `test_pressure_derivative` | `pressure_derivative_n(e, 1, dp)` matches analytic `dP/de = K*Gamma*e^(Gamma-1)` |

Note: This test requires initializing `fields_mod` global state (allocate arrays, set `eos_file`, `enthalpy_min`). A setup subroutine handles this boilerplate.

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| pFUnit not available via Homebrew on macOS | Medium | Low | Fall back to source build with `cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local` |
| CMake can't find FFTW3 (needed by cheb_mod) | Low | Medium | Use `pkg-config` via `find_package(PkgConfig)` + `pkg_check_modules(FFTW3)` |
| ad_mod uses `real(8)` not `wp` — mixed precision in tests | Low | Low | Tests for ad_mod use `real(8)` to match; document the inconsistency |
| EOS integration test fragile due to global state | Medium | Medium | Isolate in its own test executable; run after numerics tests |
| `.pf` preprocessor generates unreadable intermediate `.F90` | Low | Low | Add `build-test/` to `.gitignore` |

---

## Verification Steps

1. `make test` — all tests pass (single command)
2. `make` — production build still works (no regressions)
3. `make debug` — debug build still works
4. Introduce a deliberate bug in `ad_mod.f90` (e.g., change `+` to `-` in `add_dd`) — verify test catches it
5. `build-test/` is in `.gitignore` and not tracked

---

## Execution Order

| Step | Phase | Files Created/Modified | Depends On |
|------|-------|----------------------|------------|
| 1.1 | Infra | (pFUnit install) | — |
| 1.2 | Infra | `CMakeLists.txt` | 1.1 |
| 1.3 | Infra | `tests/CMakeLists.txt` | 1.2 |
| 1.4 | Infra | `Makefile` (append test target) | 1.2 |
| 2.1 | Tier 1 | `tests/test_ad_mod.pf` | 1.3 |
| 2.2 | Tier 1 | `tests/test_toolkit_helpers.pf` | 1.3 |
| 3.1 | Tier 2 | `tests/test_spline_mod.pf` | 1.3 |
| 3.2 | Tier 2 | `tests/test_cheb_mod.pf` | 1.3 |
| 3.3 | Tier 2 | `tests/test_spectral_hub.pf` | 1.3 |
| 3.4 | Tier 2 | `tests/test_brent_mod.pf` | 1.3 |
| 4.1 | Tier 3 | `tests/fixtures/test_poly.dat` | 1.3 |
| 4.2 | Tier 3 | `tests/test_eos_mod.pf` | 4.1 |

---

*Plan generated: 2026-04-10*
