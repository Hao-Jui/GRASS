# Re-implement 8a7b201 in Validated Steps

**Date:** 2026-04-10
**Branch:** `rebuild` (at `d2acb03`)
**Baseline:** `Cont/properties_baseline.dat` from GR+uryu run at d2acb03
**Validation after each step:** `bash sub.sh && diff Cont/properties.dat Cont/properties_baseline.dat`

---

## What 8a7b201 did (39 files, +1773/-1322)

Grouped by risk:

### Zero-risk: Pure filename renames (no code change)
These can't affect output — the module names inside stay the same.
- `Ope_eq.f90` → `ope_eq_mod.f90` (0 lines changed)
- `sphere.f90` → `sphere_mod.f90` (0 lines changed)

### Low-risk: Renames + trivial import updates
Module name or filename changes with only `use` statement updates.
- `eos.f90` → `eos_mod.f90` (4 lines)
- `grid.f90` → `grid_mod.f90` (2 lines)
- `spin_derivatives.f90` → `spin_derivatives_mod.f90` (4 lines)
- `spin_updates.f90` → `spin_updates_mod.f90` (4 lines)
- `spin_workspace.f90` → `spin_workspace_mod.f90` (12 lines)
- `spectral_hub.f90` → `spectral_hub_mod.f90` (6 lines)
- `toolkit_mod.f90` (16 lines — import updates)
- `rotational_law_mod.f90` (2 lines)
- `relaxation_mod.f90` (12 lines)
- `constraint_eq.f90` → `constraint_mod.f90` (28 lines)
- `miscellaneous.f90` → `miscellaneous_mod.f90` (37 lines)
- `set_disk.f90` (unchanged but depends on renamed modules)
- `scalar_burning.f90` → `scalar_burning_mod.f90` (2 lines)
- `starting_model.f90` → `starting_model_mod.f90` (10 lines)
- `para_mod.f90` (214 lines — import updates + `initialize_theory` additions)
- `main.f90` (17 lines)
- `Makefile` (94 lines — source list + dependency rules)

### Medium-risk: New functionality (additive, shouldn't affect existing output)
- `donutization_mod.f90` — already on rebuild branch ✓
- `exporter_mod.f90` — already on rebuild branch ✓
- `regrid_mod.f90` — replaces `restart_read.f90` (262 lines new vs 393 deleted)
- `cheb_mod.f90` — 115 lines added (new Chebyshev features)

### HIGH-risk: Logic changes that could affect convergence
- `brent.f90` → `brent_mod.f90` — **complete rewrite** (+308/-225). THIS broke uryu.
- `spin_relaxation.f90` → `spin_relaxation_mod.f90` — **137 lines changed**. Relaxation loop modifications.
- `spin_integration.f90` → `spin_integration_mod.f90` — **249 lines changed**. Integration + output refactoring.
- `analysis_v2.f90` → `analysis_mod.f90` — **208 lines changed**. Solution properties refactoring.
- `shoot_v2.f90` → `shoot_mod.f90` — **149 lines changed**. Shooting solver refactoring.
- `shoot_solver_1d_r_ratio_mod.f90` — 36 lines changed
- `shoot_solver_2d_mod.f90` — 14 lines changed
- `shoot_solver_1d_hc_mod.f90` — 6 lines changed

---

## Execution Steps

Each step: apply change → `bash sub.sh` → `diff Cont/properties.dat Cont/properties_baseline.dat` → commit if identical.

### Step 1: Pure file renames (zero risk)

Rename files where only the filename changes, no code edits:
```
mv src/core/Ope_eq.f90      src/core/ope_eq_mod.f90
mv src/core/sphere.f90       src/core/sphere_mod.f90
```
Update Makefile references. Validate.

### Step 2: Module renames — batch 1 (core modules)

Apply filename + module-name renames for core modules that have minimal code changes:
- `eos.f90` → `eos_mod.f90`
- `grid.f90` → `grid_mod.f90`
- `constraint_eq.f90` → `constraint_mod.f90`
- `miscellaneous.f90` → `miscellaneous_mod.f90`
- `scalar_burning.f90` → `scalar_burning_mod.f90`

Cherry-pick ONLY the rename + import-update changes from 8a7b201 for these files. No logic changes. Update Makefile. Validate.

### Step 3: Module renames — batch 2 (theory modules)

- `spin_derivatives.f90` → `spin_derivatives_mod.f90`
- `spin_updates.f90` → `spin_updates_mod.f90`
- `spin_workspace.f90` → `spin_workspace_mod.f90`
- `spectral_hub.f90` → `spectral_hub_mod.f90`
- `relaxation_mod.f90` (import updates only)
- `rotational_law_mod.f90` (import updates only)

Update Makefile. Validate.

### Step 4: Solver file renames

- `starting_model.f90` → `starting_model_mod.f90`
- `shoot_v2.f90` → `shoot_mod.f90`
- `MRcurve.f90` → `MRcurve_mod.f90`

Apply ONLY rename + import changes. NO logic changes to shoot_mod or starting_model_mod. Update Makefile. Validate.

### Step 5: regrid_mod (replaces restart_read)

- Delete `restart_read.f90`
- Add `regrid_mod.f90` from 8a7b201
- Update Makefile and imports in `starting_model_mod.f90`

This changes the regrid path but NOT the default-mode path (which our test uses). Validate.

### Step 6: spin_relaxation_mod — rename only

- `spin_relaxation.f90` → `spin_relaxation_mod.f90`
- Apply ONLY module-name rename + import updates
- Do NOT apply any logic changes (no Anderson-hurt tracking, no cycle-lag changes)

Validate. **This is where the `.and.` UB previously triggered — with rename-only it should pass.**

### Step 7: spin_integration_mod — rename + safe changes

- `spin_integration.f90` → `spin_integration_mod.f90`
- Apply rename + import updates
- Apply safe structural changes (output formatting, etc.)
- Do NOT apply changes that alter the numerical integration path

Validate.

### Step 8: analysis_mod — rename + safe changes

- `analysis_v2.f90` → `analysis_mod.f90`
- Apply rename + import updates
- Apply safe refactoring (variable naming, formatting)
- Preserve numerical output computation

Validate.

### Step 9: cheb_mod additions

- Apply the 115 new lines to `cheb_mod.f90`
- These are new Chebyshev features — additive, shouldn't affect existing paths

Validate.

### Step 10: para_mod + main.f90 updates

- Apply import updates to `para_mod.f90` for new module names
- Apply changes to `main.f90`
- Apply `initialize_theory` refactoring in `para_mod.f90`

Validate.

### Step 11: .gitignore

- Add `.gitignore` from 8a7b201

No validation needed (doesn't affect build).

### Step 12 (HIGH RISK): spin_relaxation_mod logic changes

Apply the relaxation loop changes from 8a7b201:
- `ensure_allocated` restructuring
- `reset_accel_state` extraction
- Anderson-hurt tracking
- Cycle-lag changes

Validate. **If this fails, revert and apply changes one at a time.**

### Step 13 (HIGH RISK): brent_mod rewrite

Replace `brent.f90` with `brent_mod.f90` from 8a7b201.

Validate. **If this fails, keep the d2acb03 brent and skip — we know the old brent works for uryu.**

### Step 14 (MEDIUM RISK): shoot_mod logic changes

Apply the shooting solver refactoring from 8a7b201.

Validate.

### Step 15 (MEDIUM RISK): shoot_solver_* changes

Apply changes to the 1d/2d shooting solver modules.

Validate.

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|--------------|
| AC1 | After each step, `diff Cont/properties.dat Cont/properties_baseline.dat` is empty | Zero diff |
| AC2 | After all steps, `make` builds cleanly | 0 errors |
| AC3 | GR+uryu converges in < 2s at res=400 | Elapsed time in output |
| AC4 | Final properties match baseline to all digits | diff is empty |

---

## If a step fails

1. `git diff` to see what changed
2. Revert: `git checkout -- .`
3. Apply the step's changes more granularly (split into sub-steps)
4. If a specific logic change can't be reconciled with the baseline, **skip it** — the d2acb03 version is correct

---

*Plan generated: 2026-04-10*
