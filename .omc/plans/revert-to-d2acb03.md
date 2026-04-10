# Revert to d2acb03 and Rebuild

**Date:** 2026-04-10
**Goal:** Reset codebase to commit `d2acb03` (last version where all 3 solver_types converge), preserving donutization_mod and exporter_mod from HEAD.

---

## Context

- `d2acb03`: uniform 0.03s, const_j 0.09s, uryu 0.30s — all work
- HEAD: uniform works, uryu broken by layered regressions across 5 commits
- Debugging the regression is harder than reverting + re-adding the 2 new modules

## What d2acb03 has

```
src/
├── core/     para_mod.f90, eos.f90, grid.f90, sphere.f90, constraint_eq.f90,
│             Ope_eq.f90, analysis_v2.f90, starting_model.f90, shoot_v2.f90,
│             MRcurve.f90, miscellaneous.f90, set_disk.f90, scalar_burning.f90,
│             precision_mod.f90, restart_read.f90, newton_mod.f90, newton1d_mod.f90
├── tool/     toolkit_mod.f90, ad_mod.f90, brent.f90, cheb_mod.f90,
│             nag_compat_mod.f90, simpson_mod.f90, spectral_hub_mod.f90
└── theory/   rotation_solver.f90, rotational_law_mod.f90, relaxation_mod.f90,
              spin_derivatives.f90, spin_updates.f90, spin_relaxation.f90,
              spin_workspace.f90, spin_integration.f90
```

All parameters in one `para_mod` module.

## What HEAD adds (to preserve)

1. **donutization_mod.f90** — computes donutization number from M_b-weighted radii
2. **exporter_mod.f90** — writes solution data to Cont/ output files

Both use `fields_mod`, `constants_mod`, `grid_config_mod` (HEAD's split modules) which map to `para_mod` at d2acb03.

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|--------------|
| AC1 | `make` builds cleanly | 0 errors |
| AC2 | GR+uryu (r_ratio=0.7, e_center=.7e15, MPA1, SDIV=401) converges < 60s | grep "Converged" |
| AC3 | GR+uniform converges | grep "Converged" |
| AC4 | GR+const_j converges | grep "Converged" |
| AC5 | donutization_mod compiles and is callable | Build succeeds with use statement |
| AC6 | exporter_mod compiles and is callable | Build succeeds with use statement |

---

## Implementation Steps

### Step 1 — Save donutization and exporter from HEAD

```bash
cp src/diagnostics/donutization_mod.f90 /tmp/donutization_mod.f90
cp src/diagnostics/exporter_mod.f90 /tmp/exporter_mod.f90
```

### Step 2 — Create new branch from d2acb03

```bash
git checkout -b rebuild d2acb03
```

This gives us the known-working codebase with the old directory structure.

### Step 3 — Port donutization_mod into d2acb03 structure

Copy `/tmp/donutization_mod.f90` → `src/tool/donutization_mod.f90`

Adapt imports — replace HEAD module names with d2acb03 equivalents:
- `use fields_mod, only: ...` → `use para_mod, only: ...`
- `use constants_mod, only: ...` → `use para_mod, only: ...`
- `use grid_config_mod, only: ...` → `use para_mod, only: ...`

All constants, parameters, and field arrays live in `para_mod` at d2acb03.

### Step 4 — Port exporter_mod into d2acb03 structure

Copy `/tmp/exporter_mod.f90` → `src/tool/exporter_mod.f90`

Same import adaptation as Step 3. Also check:
- `eos_mod` references → at d2acb03 the module might be named `eos_mod` (check) or in `eos.f90`
- Verify all symbols referenced from `para_mod` exist there

### Step 5 — Update Makefile

Add both new files to the `SOURCES` list and add dependency rules:
```makefile
$(SRC_TOOL)/donutization_mod.f90
$(SRC_TOOL)/exporter_mod.f90
```

Dependency rules:
```makefile
$(OBJDIR)/$(SRC_TOOL)/donutization_mod.o: $(OBJDIR)/$(SRC_CORE)/para_mod.o
$(OBJDIR)/$(SRC_TOOL)/exporter_mod.o: $(OBJDIR)/$(SRC_CORE)/eos.o $(OBJDIR)/$(SRC_CORE)/para_mod.o
```

### Step 6 — Wire into calling code

At d2acb03, `spin_integration.f90` and `analysis_v2.f90` are the likely call sites. Check HEAD's integration points:
- Where is `donutization` called? → likely in `spin_integration` or `analysis_v2`
- Where is `exporter` called? → likely in `spin_integration` (for output_helper)

Add `use` statements and calls at the same logical points as HEAD.

### Step 7 — Build and test all 3 solver_types

```bash
make clean && make
# Test each solver_type with GR, MPA1, r_ratio=0.7, e_center=.7e15, SDIV=401
```

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| donutization/exporter reference symbols not in para_mod | Medium | Low | Check symbol list; add any missing to para_mod |
| d2acb03 has old brent.f90 without error handling | Low | Low | The old brent works for uryu — keep it |
| Losing performance optimizations from HEAD | Medium | Low | Can selectively re-add safe optimizations later |
| Git history diverges from main | Low | Medium | Use merge/rebase when ready |

---

*Plan generated: 2026-04-10*
