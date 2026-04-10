# GRASS Codebase: 2026 Modernization Assessment

**Date:** 2026-04-10 (comprehensive re-assessment)
**Codebase:** 40 Fortran files, ~10,300 LOC (100% free-form `.f90`)
**Domain:** Rotating neutron stars in GR and Scalar-Tensor theories (field equations + relaxation solver)

---

## Overall Score: 62/100 (Good Foundation, Significant Gaps)

| Dimension | Score | Grade | Notes |
|-----------|-------|-------|-------|
| Language Modernization | 65/100 | C+ | Strong F2003 core undermined by ~40% `real(8)` / `1.d0` bypass |
| Parallelism & Performance | 40/100 | D | Serial only, no `do concurrent`, BLAS/FFTW used |
| Build System & Tooling | 50/100 | D+ | gfortran-only Makefile, no CMake/fpm, manual deps |
| Code Quality & Maintainability | 72/100 | B | Clean modules, good USE ONLY, extensive global state |
| Testing & CI/CD | 5/100 | F | No tests, no CI, no framework |
| Documentation | 45/100 | D+ | Sparse procedure docs, good README |

---

## 1. LANGUAGE MODERNIZATION (65/100) -- Grade C+

### What's good (F2003/2008 level)
- **Zero obsolescent patterns**: No COMMON, EQUIVALENCE, ENTRY, arithmetic IF, fixed-form, PAUSE, FORALL, numbered DO
- **IMPLICIT NONE**: 100% coverage across all 40 files (72 occurrences)
- **Module architecture**: 42 modules across 38 files with clean layered hierarchy
- **USE...ONLY**: 301 imports, all with `only:` restriction (0 bare `use` -- excellent)
- **INTENT declarations**: 649 occurrences across 32 files (extensive)
- **ALLOCATABLE preferred**: 76 occurrences; only 7 `pointer` uses (all in `eos_mod.f90`)
- **`source=` allocation**: 51 occurrences for initialized allocation (F2003+)
- **Abstract interfaces**: 6 files -- callback patterns for Brent, EOS, ODE, JFNK, shooting
- **Operator overloading**: Full AD module with `+`, `-`, `*`, `/`, `exp`, `log` (`ad_mod.f90:12-34`)
- **Type-bound procedures**: `laplacian_operator` with 9 bound methods (`ope_eq_mod.f90:15-34`)
- **Procedure pointers**: Runtime EOS dispatch (`eos_mod.f90:62-64`)
- **`error stop`**: 34 occurrences in 10 files with descriptive messages (F2008)
- **`block` constructs**: ~29 occurrences in 10 files for local scoping (F2008)
- **`associate`**: 10 occurrences in 4 files (F2003)
- **`contiguous`**: 13 occurrences in `relaxation_mod.f90` for BLAS performance (F2008)
- **`pure`**: 47 procedures; **`elemental`**: 5 procedures
- **`merge` intrinsic**: 24 occurrences for branchless conditionals
- **`newunit=`**: 16 of 19 `open` statements use it (modern)
- **`private`/`public` access control**: 22 modules with explicit defaults
- **`iso_fortran_env`**: 7 files (`real64`, `real128`, `error_unit`, `output_unit`)
- **`ieee_arithmetic`**: 4 files (`ieee_is_nan`, `ieee_is_finite`)
- **`iso_c_binding`**: FFTW interface in `cheb_mod.f90` with `bind(C)` (F2003)
- **Preprocessor eliminated**: Zero `#ifdef` in `src/`; legacy branches had them everywhere
- **`import ::`**: 14 occurrences for proper interface block type resolution

### Critical issues

| Issue | Count | Files | Impact |
|-------|------:|-------|--------|
| **`real(8)` instead of `real(wp)`** | 100+ | `ad_mod`, `rotational_law_mod`, `sphere_mod`, `miscellaneous_mod`, `set_disk`, `shoot_solver_1d_hc_mod`, `shoot_solver_1d_r_ratio_mod` | Breaks quad-precision builds; defeats `precision_mod` |
| **`1.d0`/`0.d0` literals** | 200+ | 15+ files (concentrated in `compat/`, `solver/`, `theory/`) | Hardcoded double; should be `1.0_wp` |
| **`dble()` intrinsic** | 32 | 10 files | Should be `real(..., wp)` |
| **Bare `stop` (no `error stop`)** | ~30 | `ope_eq_mod`, `spin_updates_mod`, `shoot_mod`, `sphere_mod`, `toolkit_mod` | Opaque failures, missing F2008 semantics |
| **Bare `stop` (no message)** | 2 | `shoot_mod.f90:113`, `shoot_solver_2d_mod.f90:250` | Zero diagnostic information |
| **`goto` statements** | 8 | `regrid_mod.f90:44,53,57,117,119,121,123,133` | Legacy error-handling jumps to labels 90/99 |
| **`isnan()` (non-standard)** | 4 | `eos_mod`, `rotation_solver_mod`, `spin_updates_mod` | GNU extension; should use `ieee_is_nan` |
| **Hardcoded unit numbers** | 3 | `miscellaneous_mod:190`, `sphere_mod:74,196` | Should use `newunit=` |
| **`dsqrt` (legacy intrinsic)** | 1 | `toolkit_mod.f90:376` | Should be `sqrt` |
| **Local `wp` redefinition** | 2 | `spin_updates_mod.f90:2`, `spin_derivatives_mod.f90:2` | `use iso_fortran_env, only: wp => real64` bypasses `precision_mod` |

### Not adopted (F2008/2018/2023)
| Feature | Status | Priority |
|---------|--------|----------|
| `do concurrent` | Not used (0 occurrences) | LOW |
| `submodule` | Not used | LOW |
| `findloc` | Not used | LOW |
| Coarrays | Not used (excluded by project policy) | N/A |
| `select type` / type extension | Not used | LOW |
| Parameterized derived types | Not used | LOW |
| Enumeration types (F2023) | Not used | N/A |

---

## 2. PARALLELISM & PERFORMANCE (40/100) -- Grade D

### What's present
- **BLAS/LAPACK**: `dgemm`, `dgemv`, `dposv`, `dgesv`, `dgels`, `dpbsv` used for matrix operations
- **FFTW3**: Chebyshev transforms via C interop (`cheb_mod.f90:25-42`)
- **Workspace caching**: `spin_workspace_mod.f90` pre-computes ~30 arrays to avoid recomputation
- **`cpu_time()`**: 38 call sites across 8 files for profiling

### Gaps
| Gap | Impact | Priority |
|-----|--------|----------|
| **Purely serial code** | Cannot exploit multi-core HPC | HIGH |
| **No `do concurrent`** | Missing compiler auto-parallelization hints | MEDIUM |
| **No OpenMP** | No thread-level parallelism | MEDIUM |
| **No MPI** | No distributed parallelism | LOW (single-star solver) |
| **LAPACK tied to `double precision`** | Would need S/D prefix switching for quad builds | MEDIUM |
| **`cpu_time` instead of `system_clock`** | Not wall-clock; wrong for parallel codes | LOW |
| **No SIMD directives** | Relies on compiler auto-vectorization | LOW |

---

## 3. BUILD SYSTEM & TOOLING (50/100) -- Grade D+

### Current state
- **GNU Make** with Debug/Release modes (Makefile:159 lines)
- **Debug flags**: `-Wall -Wextra -Wimplicit-interface -fcheck=all -Wconversion -Wuse-without-only -finit-real=nan`
- **Release flags**: `-O3 -march=native -fno-backtrace`
- **Explicit dependency tracking**: 50+ manual rules (Makefile:86-137)
- **Build directories**: `build/{obj,mod,bin}` keeps source clean
- **`pkg-config`** for FFTW discovery with fallback

### Gaps
| Gap | Impact | Priority |
|-----|--------|----------|
| **gfortran-only** (hardcoded `FC = gfortran`) | Cannot build with ifx, nvfortran, Cray, flang | HIGH |
| **No CMake or fpm** | No cross-platform portable builds | HIGH |
| **Manual dependency rules** | Adding `use` requires Makefile edit; drift-prone | MEDIUM |
| **macOS-specific linker flag** (`-Wl,-stack_size,...`) | Breaks on Linux/HPC | MEDIUM |
| **No Fortran standard enforcement** (no `-std=f2018`) | Non-standard extensions go undetected | MEDIUM |
| **Binary named `a.out`** | No install target | LOW |
| **No test target** | `make test` does not exist | HIGH |

---

## 4. CODE QUALITY & MAINTAINABILITY (72/100) -- Grade B

### Strengths
- **100% `USE...ONLY`**: Zero bare module imports (301 explicit imports)
- **Clean module hierarchy**: 7-directory layout with layered dependencies
- **Access control**: 22 modules with `private` default + explicit `public` exports
- **Workspace pattern**: `spin_workspace_mod` centralizes pre-computed arrays
- **Array operations**: 106 uses of `matmul`, `dot_product`, `sum`, `maxval`, `norm2`, `spread`, `where`
- **Whole-array arithmetic**: Extensive in `spin_workspace_mod.f90:69-91`
- **Assumed-shape arrays**: Used for procedure arguments in modern modules
- **BLAS for hot paths**: `dgemm`/`dgemv` chosen over `matmul` for large matrices
- **`epsilon()` for tolerances**: 37 occurrences for relative comparison
- **Consistent naming**: Mostly `snake_case`; `UPPER_CASE` for constants
- **Preprocessor replaced by runtime dispatch**: `active_theory`, `solver_type`, `run_task`

### Remaining issues
| Issue | Severity | Details |
|-------|----------|---------|
| **Global mutable state** | HIGH | ~120+ mutable module variables in `fields_mod` (~60 vars), `grid_config_mod`, `spin_workspace_mod` (~30 arrays), `eos_mod`, `rotation_law_mod`, `spin_relaxation_mod` |
| **`save` attribute proliferation** | MEDIUM | 50+ `save` declarations in 9 files; `spin_relaxation_mod:31-47` has 18 `save` variables making the module non-reentrant |
| **LAPACK interface duplication** | MEDIUM | Same BLAS routines declared in `relaxation_mod`, `shoot_solver_2d_mod`, `cheb_mod`, `spin_derivatives_mod`, `spin_integration_mod` separately |
| **PI defined in 3 places** | LOW | `constants_mod`, `cheb_mod`, `spectral_hub_mod` |
| **Typos in identifiers** | LOW | `find_omege_e` (should be `omega`), `uyru_p`/`uyru_q` (should be `uryu`) |
| **Legacy function names** | LOW | `d01gaf`, `d02pcf` (NAG), `plgndr` (Numerical Recipes) |
| **No `FINAL` subroutines** | LOW | Manual `deallocate_*` routines used instead of automatic cleanup |

---

## 5. TESTING & CI/CD (5/100) -- Grade F

### Current state
- **No test framework**: No pFUnit, FRUIT, or Vegetables
- **No test files**: No `*_test.f90` anywhere in the codebase
- **No CI/CD**: No GitHub Actions, GitLab CI, or any pipeline
- **No `make test` target**
- **No regression test infrastructure**: Output comparison is manual
- **Debug mode sanitizers**: `-fcheck=all -finit-real=nan` catches some issues at runtime (the only safety net)

### What could be tested
| Module | Test type | Effort |
|--------|-----------|--------|
| `spline_mod` | Polynomial exactness, endpoint BCs | Low |
| `brent_mod` | Known root finding problems | Low |
| `ad_mod` | Derivative correctness vs FD | Low |
| `cheb_mod` | Interpolation convergence | Low |
| `toolkit_mod` | Legendre functions, interpolation | Low |
| `eos_mod` | Polytropic round-trip, table loading | Medium |
| `grid_mod` | Grid spacing verification | Low |

---

## 6. DOCUMENTATION (45/100) -- Grade D+

### What exists
- **README.md**: Comprehensive physics and usage guide (project-level)
- **`src/theory/readme.md`**: Theory-specific notes
- **Module-level comments**: Present in `spline_mod` (13-line header), `toolkit_mod` (11-line header), `brent_mod` (algorithm phases), `ad_mod` (purpose statement)
- **Section separators**: Used in `toolkit_mod`, `brent_mod` for logical grouping

### What's missing
| Gap | Priority |
|-----|----------|
| **No FORD/Doxygen markup** (`!>`) | MEDIUM |
| **Most procedures undocumented** | MEDIUM |
| **Theory modules almost no comments** | MEDIUM |
| **No parameter file documentation** with units/valid ranges | LOW |
| **No contributing guide** | LOW |

---

## REMAINING MODERNIZATION ROADMAP

### Priority 1 -- Fix immediately (low effort, high value)

1. **Unify precision system**: Replace all `real(8)` with `real(wp)` in `ad_mod`, `rotational_law_mod`, `sphere_mod`, `miscellaneous_mod`, `set_disk`, both `shoot_solver_1d_*` modules (~100+ replacements)
2. **Fix literal constants**: Replace `1.d0`/`0.d0`/`1.d-12` with `1.0_wp`/`0.0_wp`/`1.e-12_wp` (~200+ replacements)
3. **Replace `dble()`**: Change to `real(..., wp)` (32 occurrences)
4. **Fix local `wp` redefinition**: `spin_updates_mod.f90:2` and `spin_derivatives_mod.f90:2` should import from `precision_mod`, not redefine locally
5. **Replace bare `stop`**: Change ~30 occurrences to `error stop "module: message"`, especially the 2 bare `stop` with no message
6. **Replace `isnan()`**: Change 4 occurrences to `ieee_is_nan()` for portability
7. **Replace `dsqrt`**: Change to `sqrt` at `toolkit_mod.f90:376`
8. **Replace hardcoded unit numbers**: 3 in `miscellaneous_mod`/`sphere_mod` to `newunit=`

### Priority 2 -- Short term (medium effort)

9. **Add `-std=f2018` to Makefile**: Catch non-standard extensions at compile time
10. **Consolidate LAPACK interfaces**: Single `lapack_interfaces_mod` instead of 5 duplicate declarations
11. **Replace `goto`**: Refactor 8 occurrences in `regrid_mod.f90` to structured error handling
12. **Unify PI definition**: Single source in `constants_mod`, import everywhere
13. **Fix typos**: `find_omege_e` -> `find_omega_e`, `uyru_*` -> `uryu_*`
14. **Add procedure-level docstrings**: At minimum one-sentence purpose for all public procedures

### Priority 3 -- Medium term (larger effort, high payoff)

15. **Add CMake build system**: Multi-compiler support (gfortran, ifx, nvfortran, flang), `-std=f2018` enforced, install target
16. **Add pFUnit test suite**: Cover `numerics/` (spline, Brent, AD, Chebyshev) and `eos/` at minimum
17. **Add GitHub Actions CI**: Build in Debug + Release, run tests on PR
18. **Add `do concurrent`**: Apply to safe initialization loops in `spin_workspace_mod`, `grid_mod`

### Priority 4 -- Long term (high effort, optional)

19. **Reduce global state**: Encapsulate `fields_mod` into a derived type passed by argument
20. **Submodule split**: Large modules (`eos_mod` 789L, `spin_integration_mod` 751L, `toolkit_mod` 652L)
21. **FORD documentation**: Auto-generated API docs
22. **Add `FINAL` subroutines**: For types that own allocatable members

---

## COMPARISON: GRASS vs 2026 Best Practices

| Practice | GRASS Status | 2026 Standard |
|----------|-------------|---------------|
| Free-form source | YES | YES |
| IMPLICIT NONE | YES (100%) | YES |
| Modules with USE...ONLY | YES (100%) | YES |
| ISO_FORTRAN_ENV kinds | PARTIAL (~60% of code) | YES |
| `real(wp)` everywhere | NO (~40% uses `real(8)`) | YES |
| `_wp` literal suffixes | PARTIAL (~60%) | YES |
| ERROR STOP | PARTIAL (34/64 stops) | YES |
| INTENT declarations | YES (extensive) | YES |
| Private/public access | YES (22 modules) | YES |
| Abstract interfaces | YES | YES |
| Pure/elemental | YES (52 procedures) | YES |
| DO CONCURRENT | NO | Recommended |
| Submodules | NO | Recommended |
| Type finalization | NO | Recommended |
| Zero GOTO | NO (8 remain) | YES |
| Zero preprocessor | YES | YES |
| CMake/fpm build | NO | YES |
| `-std=f2018` enforced | NO | YES |
| Multi-compiler support | NO (gfortran only) | YES |
| CI/CD pipeline | NO | YES |
| Unit testing | NO | YES |
| FORD documentation | NO | Recommended |

---

## VERDICT

GRASS has a **solid Fortran 2003/2008 foundation** -- fully free-form, modular, with no obsolescent features, universal `implicit none`, 100% `USE...ONLY`, extensive `intent` declarations, and good use of abstract interfaces, operator overloading, and type-bound procedures. The elimination of all preprocessor directives in favor of runtime dispatch is exemplary.

However, the codebase scores **62/100** due to three systemic gaps:

1. **Precision inconsistency (most critical)**: ~40% of the code bypasses the `precision_mod` kind system, using hardcoded `real(8)` and `1.d0` literals. The `ad_mod` -- a critical infrastructure module providing automatic differentiation for EOS and toolkit -- is entirely `real(8)`. This makes the `real128` quad-precision capability advertised in `precision_mod.f90:4` non-functional and creates silent precision-mixing bugs if `wp` is ever changed.

2. **No testing or CI**: Zero automated tests of any kind. For a scientific code where numerical correctness is paramount, this is the highest-risk gap. Every numerical module is a candidate for polynomial-exactness or convergence-order tests.

3. **Build portability**: gfortran-only Makefile with macOS-specific flags, manual dependency tracking, and no standard enforcement (`-std=f2018`). Cannot build on HPC systems with Intel or NVIDIA compilers without manual flag translation.

**Priority 1 items** (precision unification, `error stop`, `isnan`) are mechanical find-replace operations that would immediately raise the Language Modernization score from 65 to ~85. Combined with CMake and a basic test suite (Priority 3), the overall score would reach **~80/100** -- on par with well-maintained 2026 scientific Fortran projects.

**Bottom line**: Clean architecture, strong module design, no legacy debt -- but the precision system is only half-adopted, and the complete absence of testing and CI is the single largest risk for a code that depends on numerical accuracy.

---

*Assessment generated: 2026-04-10. Methodology: automated codebase exploration with line-level evidence.*
