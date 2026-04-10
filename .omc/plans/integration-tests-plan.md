# Integration Tests Plan — Two Production Runs

**Date:** 2026-04-10
**Goal:** Add two full-pipeline integration tests (ST and GR+uryu) to the pFUnit/CMake test suite.

---

## Requirements Summary

Add two CTest integration tests that compile and run the full GRASS solver with:
1. **ST test** — `active_theory = THEORY_ST`, default `solver_type = "uniform"`
2. **GR+uryu test** — `active_theory = THEORY_GR`, `solver_type = "uryu"`, `r_ratio = 0.8`

Both use reduced grid resolution (`res = 100`, giving SDIV=201) to keep test runtime under 10 minutes.

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|--------------|
| AC1 | `grass_test_st` executable builds without errors | Build log: 0 errors |
| AC2 | `grass_test_gr_uryu` executable builds without errors | Build log: 0 errors |
| AC3 | ST test converges and prints "Converged" | `ctest -R integration_st`: PASS |
| AC4 | GR+uryu test converges and prints "Converged" | `ctest -R integration_gr_uryu`: PASS |
| AC5 | Both tests exit with code 0 | CTest default check |
| AC6 | Both tests complete within 600 seconds | `TIMEOUT 600` property |
| AC7 | GR+uryu test uses r_ratio=0.8 (not 0.9) | Grep starting_model_uryu.f90 line 38 |
| AC8 | Both tests run at res=100, not res=500 | Grep grid_config_test.f90 line 6 |
| AC9 | `make test` still passes (all 7 unit + 2 integration) | 9/9 tests pass |
| AC10 | Production `make` build is unaffected | `make` succeeds |

---

## Implementation Steps

### Step 1 — Create test grid config (shared)

**File:** `tests/integration/grid_config_test.f90`
**Source:** Copy of `src/macros/grid_config_mod.f90` (53 lines) with one change:
- Line 6: `integer, parameter :: res = 100` (was `500`)

Everything else identical. This gives SDIV=201 × MDIV=41 = 8,241 grid points (~25× faster than production).

### Step 2 — Create GR+uryu para_panel

**File:** `tests/integration/para_panel_gr_uryu.f90`
**Source:** Copy of `src/para_panel.f90` with two changes:
- Line 9: `integer :: active_theory = THEORY_GR` (was `THEORY_ST`)
- Line 19: `character(len=20) :: solver_type = "uryu"` (was `"uniform"`)

### Step 3 — Create GR+uryu starting model

**File:** `tests/integration/starting_model_uryu.f90`
**Source:** Copy of `src/solver/starting_model_mod.f90` (81 lines) with one change:
- Line 38: `r_ratio = 0.8e0_wp` (was `merge(0.9e0_wp, 1.e0_wp, ...)`)

### Step 4 — Add integration test targets to CMakeLists.txt

**File:** `CMakeLists.txt` (root, modify)

Add two integration test executables. Each compiles ALL source files but substitutes the test-specific config files. Use a CMake variable for the common source list (everything except grid_config_mod, para_panel, starting_model_mod, main).

```cmake
# Common sources shared by both integration tests (no grid_config, no para_panel,
# no starting_model, no main)
set(INTEGRATION_COMMON_SOURCES
  src/macros/precision_mod.f90
  src/macros/constants_mod.f90
  src/numerics/ad_mod.f90
  src/numerics/nag_compat_mod.f90
  src/numerics/brent_mod.f90
  src/numerics/cheb_mod.f90
  src/numerics/spline_mod.f90
  src/numerics/spectral_hub_mod.f90
  src/numerics/toolkit_mod.f90
  src/eos/eos_mod.f90
  src/grid/grid_mod.f90
  src/diagnostics/ope_eq_mod.f90
  src/diagnostics/constraint_mod.f90
  src/diagnostics/donutization_mod.f90
  src/diagnostics/exporter_mod.f90
  src/diagnostics/analysis_mod.f90
  src/compat/miscellaneous_mod.f90
  src/compat/set_disk.f90
  src/compat/sphere_mod.f90
  src/compat/regrid_mod.f90
  src/theory/relaxation_mod.f90
  src/theory/rotational_law_mod.f90
  src/theory/spin_derivatives_mod.f90
  src/theory/spin_updates_mod.f90
  src/theory/spin_workspace_mod.f90
  src/theory/spin_integration_mod.f90
  src/theory/spin_relaxation_mod.f90
  src/theory/rotation_solver_mod.f90
  src/solver/shoot_solver_2d_mod.f90
  src/solver/shoot_solver_1d_hc_mod.f90
  src/solver/shoot_solver_1d_r_ratio_mod.f90
  src/solver/scalar_burning_mod.f90
  src/solver/shoot_mod.f90
  src/solver/MRcurve_mod.f90
)

# -- ST integration test (THEORY_ST, uniform rotation) --
add_executable(grass_test_st
  tests/integration/grid_config_test.f90
  src/para_panel.f90                       # production defaults = THEORY_ST
  ${INTEGRATION_COMMON_SOURCES}
  src/solver/starting_model_mod.f90        # production starting model
  src/main.f90
)
target_include_directories(grass_test_st PRIVATE ${CMAKE_BINARY_DIR}/mod_st)
set_target_properties(grass_test_st PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_st)
target_link_libraries(grass_test_st ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${FFTW3_LIBRARIES})
target_link_directories(grass_test_st PRIVATE ${FFTW3_LIBRARY_DIRS})

# -- GR+Uryu integration test (THEORY_GR, uryu rotation, r_ratio=0.8) --
add_executable(grass_test_gr_uryu
  tests/integration/grid_config_test.f90
  tests/integration/para_panel_gr_uryu.f90 # THEORY_GR + solver_type="uryu"
  ${INTEGRATION_COMMON_SOURCES}
  tests/integration/starting_model_uryu.f90 # r_ratio=0.8
  src/main.f90
)
target_include_directories(grass_test_gr_uryu PRIVATE ${CMAKE_BINARY_DIR}/mod_gr)
set_target_properties(grass_test_gr_uryu PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_gr)
target_link_libraries(grass_test_gr_uryu ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${FFTW3_LIBRARIES})
target_link_directories(grass_test_gr_uryu PRIVATE ${FFTW3_LIBRARY_DIRS})
```

Key details:
- Each executable gets its own `Fortran_MODULE_DIRECTORY` (`mod_st`, `mod_gr`) to prevent `.mod` file collisions between the two builds
- The pFUnit `grass_src` library uses `${CMAKE_BINARY_DIR}/mod` — no collision
- Both link LAPACK, BLAS, FFTW3 directly (not via grass_src)

### Step 5 — Register CTest integration tests

**File:** `tests/CMakeLists.txt` (modify, append)

```cmake
# -- Integration tests (full solver runs) --
add_test(NAME integration_st
  COMMAND grass_test_st
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(integration_st PROPERTIES
  TIMEOUT 600
  PASS_REGULAR_EXPRESSION "Converged"
  LABELS "integration"
)

add_test(NAME integration_gr_uryu
  COMMAND grass_test_gr_uryu
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(integration_gr_uryu PROPERTIES
  TIMEOUT 600
  PASS_REGULAR_EXPRESSION "Converged"
  LABELS "integration"
)
```

The `LABELS "integration"` allows running only integration tests: `ctest -L integration`.

### Step 6 — Create Cont/ directory for output

Both tests write output files to `./Cont/`. This directory already exists (637 files). The tests will produce unique filenames based on their parameters. No cleanup needed — test output files are small at res=100.

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Solver doesn't converge at res=100 | Medium | High | Try res=150 or res=200; increase timeout |
| Test writes collide in Cont/ with existing data | Low | Low | Test filenames are unique (different parameters); res=100 output is tiny |
| Two full compilations slow down `make test` | Low | Medium | ~30s extra build time; integration tests labeled separately (`ctest -L unit`) |
| Copied source files drift from production | Medium | Medium | Note in plan: update test copies when production files change |

---

## Verification Steps

1. `cmake --build build-test` — 0 errors, both executables built
2. `ctest --test-dir build-test -R integration_st` — PASS
3. `ctest --test-dir build-test -R integration_gr_uryu` — PASS  
4. `ctest --test-dir build-test` — 9/9 tests pass (7 unit + 2 integration)
5. `make` — production build unaffected
6. `grep "res = 100" tests/integration/grid_config_test.f90` — confirms reduced resolution
7. `grep "0.8" tests/integration/starting_model_uryu.f90` — confirms r_ratio

---

## File Summary

| File | Action | Lines Changed |
|------|--------|--------------|
| `tests/integration/grid_config_test.f90` | New (copy) | 1 line: res=100 |
| `tests/integration/para_panel_gr_uryu.f90` | New (copy) | 2 lines: THEORY_GR, solver_type="uryu" |
| `tests/integration/starting_model_uryu.f90` | New (copy) | 1 line: r_ratio=0.8 |
| `CMakeLists.txt` | Modify | ~30 lines added (two executable targets) |
| `tests/CMakeLists.txt` | Modify | ~15 lines added (two add_test + properties) |

---

*Plan generated: 2026-04-10*
