# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**For physics, features, and how to use GRASS, see [README.md](README.md).**

## Code Style

- **2-space indentation** everywhere (strictly enforced)
- **`block` constructs** only when introducing local variables for scoping; never for grouping procedure calls
- **Preserve variable names** — never refactor to different naming conventions
- **Do not add comments** to `para_mod.f90` or parameter files that change frequently

## Module Organization & Dependencies

### Core hierarchy
- `para_mod.f90` — central: all parameters, grid arrays, field allocations (no dependencies on other modules except `precision_mod`)
- `eos.f90`, `grid.f90` — depend on `para_mod`; must be initialized before solver
- `constraint_eq.f90`, `analysis_v2.f90` — depend on `para_mod` and utilities; used for diagnostics

### Solver entry point
- `starting_model.f90`, `shoot_v2.f90` — initialize and run one of three tasks (OneModel, shoot, MRcurve)
- These call into the rotation solver (`rotation_solver.f90`) which orchestrates the main loop

### Theory/rotation solver chain
```
rotation_solver.f90  (outer loop on r_e)
  ├─ spin_integration.f90  (compute Green's function targets)
  ├─ spin_relaxation.f90 + relaxation_mod.f90  (4-stage: Picard → Chebyshev → Anderson → Aitken)
  └─ spin_updates.f90  (apply targets to fields)
     └─ eos.f90, rotational_law_mod.f90  (EOS lookups, rotation law)
```

### Utilities
- `toolkit_mod.f90` — interpolation, derivatives, integral accumulation (used throughout)
- `ad_mod.f90` — automatic differentiation for dP/dE (used in EOS and constraint)
- `nag_compat_mod.f90`, `cheb_mod.f90`, `spectral_hub.f90` — numerical bases and quadrature

## Key Design Patterns

**Global state in para_mod**: All grid arrays, field allocations, and parameters live here. This is by design — simplifies passing data through deeply nested solver routines without modifying interfaces.

**Workspace allocation**: `spin_workspace.f90` pre-allocates geometry caches and workspace arrays once per grid; reused across solver iterations for performance.

**Caching strategy**: Bessel function tables and geometry pre-computations cached per iteration. Derivative arrays recomputed each iteration (balances memory vs. compute).

**Log-space EOS**: All EOS quantities stored and interpolated in log-space internally. Conversions to/from physical units happen at lookup time.

**No explicit module cleanup**: Fortran automatic deallocation on program exit; no manual cleanup routines.

## Common Tasks

| Task | File | Notes |
|------|------|-------|
| Change theory (GR ↔ ST) | `para_mod.f90` line 7 | `active_theory = THEORY_GR` or `THEORY_ST` |
| Change resolution | `para_mod.f90` lines 34–37 | `res`, `MDIV`, `s_pwr` |
| Change EOS table | `para_mod.f90` line 40 | `eos_file = "APR1"` etc. |
| Add rotation law | `rotational_law_mod.f90` | Update `solver_type` selector in `rotation_solver.f90` |
| Tune relaxation | `relaxation_mod.f90`, `spin_relaxation.f90` | Adjust `PICARD_THRESH`, `CHEB_THRESH`, `w_picard`, `N_MONOTONE` |
| Add diagnostic | `analysis_v2.f90` | Use `integrate_profiles()` from `toolkit_mod.f90` |
| Modify field update | `spin_updates.f90` | Called during relaxation; updates metric/scalar fields |

## Important Caveats

- **Compilation order matters**: Makefile has explicit dependency rules (lines 72–100). If adding a new module, add it to `SOURCES` and add dependency lines.
- **Grid indices**: SDIV (radial) × MDIV (angular); Fortran column-major storage.
- **Fortran 2003 required**: No newer features; explicit interfaces used for BLAS/LAPACK.
- **No parallelization**: Single-threaded; no OpenMP or MPI. Not planned for this codebase.

## Testing & Validation

No automated test suite. Validation:
- Check Hamiltonian constraint norm (should be O(1e-10) or smaller at convergence)
- Compare M–r curves against published sequences
- Inspect diagnostic files (`Cont/Omega.dat`, `Cont/velocity.dat`) for physical bounds
