# Fix GR+Uryu Convergence Regression

**Date:** 2026-04-10
**Goal:** Restore GR+uryu convergence on the current dev branch.
**Test case:** r_ratio=0.7, e_center=.7e15, MPA1, SDIV=401, COLLOCATION_LEG.
**Baseline:** Commit `d2acb03` converges in 0.3s; current HEAD does not converge.

---

## Key Constraint (from user)

> `solver_type=uniform` and `const_j` **work** at commit `8a7b201`.

This means the relaxation loop, EOS module, Anderson acceleration, and general solver infrastructure are **not the cause**. The bug is isolated to the **uryu-specific code path**: the Brent root-finder (`find_omege_e` → `brent_core`) and/or the uryu rotation law function.

Only the uryu rotation law calls `find_omege_e`. Uniform and const_j compute angular velocity directly without root-finding.

---

## Root Cause Analysis

### Already fixed (not the cause of uryu failure)
1. **EOS double-log mismatch** (`eos_mod.f90:131`) — reverted to single-log. Confirmed working: spherical guess produces correct mass.
2. **`.or.` UB in `regrid_mod.f90:145`** — refactored to nested if/block.
3. **Missing `MODE_DEFAULT` import** in `para_panel.f90:3` — added.

### Prime suspect: `brent_mod.f90` rewrite

The Brent root-finder was completely rewritten across commits `6bc4d3d` and `44e4cae`. Batch test results confirm:
- `d7a3772` (before brent rewrite): **GOOD** (8.5s)
- `44e4cae` ("Corrected brent.f90"): **BAD** (bracket failure)

The rewrite changed:
1. **Bracket search algorithm**: Old code used geometric scaling with fallback to linear grid sweep. New code uses triple-nested `find_bracket` with `evaluate_point` wrappers.
2. **Scale factors**: `find_omege_e` now calls `brent_core(x_guess, 1.2, 1.1, ...)` — these control how aggressively the bracket expands. The old scale factors may have been different.
3. **Convergence check**: Changed from `abs(fb) < epsilon(fb)` (relative) to `fb == 0.0_wp` (exact). This is benign since bracket-width convergence `abs(xm) <= tol1` is the primary criterion.
4. **IEEE checks**: New code calls `ieee_is_finite(fx)` on every function evaluation. If the rotation law returns Inf/NaN for some x values, the old code would continue (possibly finding a valid bracket) while the new code rejects those points.
5. **Logfile write on failure**: `find_omege_e` writes to `./Cont/diff_rotation.dat` on bracket failure — this file can be inspected.

### Secondary suspect: `rotational_law_mod.f90`

The explore agent found only a no-op change (`+ 0.d0 * (...)`), but this should be verified by diff.

---

## Acceptance Criteria

| # | Criterion | Verification |
|---|-----------|--------------|
| AC1 | GR+uryu with r_ratio=0.7, e_center=.7e15, MPA1, res=200 prints "Converged" within 60s | `timeout 60 ./build/bin/a.out \| grep Converged` |
| AC2 | Axial ratio in output = 0.700 | grep output |
| AC3 | Debug build (`-fcheck=all`) runs without runtime errors | Exit code 0 |
| AC4 | ST+uniform still converges (no regression) | Quick test |
| AC5 | GR+const_j still converges (no regression) | Quick test |

---

## Implementation Plan

### Step 1 — Diff the brent rewrite precisely

Compare `d7a3772:src/tool/brent.f90` (last working) vs current `src/numerics/brent_mod.f90`. Focus on:
- `find_omege_e` wrapper: scale factors (`scale_up`, `scale_down`), `small_guess` parameter
- `brent_core`: bracket search loop logic (`find_bracket` subroutine)
- `evaluate_point`: IEEE finite check — does it reject valid evaluations?
- The old brent: did it have `small_guess`? Did it search differently?

### Step 2 — Inspect the bracket failure diagnostic

Run the current code and check `./Cont/diff_rotation.dat` — the Brent solver writes a log of function evaluations on bracket failure. This shows:
- The search range `[x_lo, x_hi]`
- Function values `f(x)` at sampled points
- Whether the function ever crosses zero in the search range

If the function never crosses zero, the uryu rotation law parameters are outside the valid domain for this stellar model. If it does cross zero but the bracket search misses it, the bracket algorithm is buggy.

### Step 3 — Port the old brent bracket logic

If the bracket algorithm is the issue, replace `find_bracket` in `brent_mod.f90` with the bracket search from `d7a3772:src/tool/brent.f90`. Keep the modern error handling (`ierr`, `errmsg`) and `evaluate_point` wrapper.

Specifically:
- Extract the old bracket search loop from `d7a3772:src/tool/brent.f90`
- Replace the current `find_bracket` subroutine body in `src/numerics/brent_mod.f90`
- Adapt variable names to match current code
- Keep the `ieee_is_finite` check (it's a safety improvement)

### Step 4 — Alternative: adjust bracket search parameters

If the bracket algorithm is correct but the search range is too narrow:
- In `find_omege_e` (line 251): change `scale_up=1.2, scale_down=1.1` to more aggressive values (e.g., `2.0, 2.0`)
- Or increase `MAX_ITER_BRACKET` from 200 to 500
- Or add a `small_guess` that better matches the typical omega_e scale for uryu

### Step 5 — Verify rotational_law_mod.f90

Diff `d7a3772:src/theory/rotational_law_mod.f90` vs current `src/theory/rotational_law_mod.f90` to confirm no logic changes beyond the no-op addition.

### Step 6 — Build, test, verify

1. `make clean && make` — release build
2. `make MODE=Debug` — debug build, run, no runtime errors
3. GR+uryu: converges with "Converged" within 60s
4. ST+uniform: still converges
5. GR+const_j: still converges

---

## Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Old bracket logic doesn't compile with new code structure | Low | Low | Small adaptation; the core algorithm is standalone |
| Adjusting scale factors fixes uryu but breaks other uses of brent | Low | Medium | Test all three rotation laws after change |
| The function truly has no root for this parameter set | Low | High | Inspect `diff_rotation.dat` to verify; compare with old code's function values |

---

## Execution Order

1. Step 2 first (inspect diagnostic) — cheapest information
2. Step 1 + Step 5 (diff brent + rotational_law) — understand the change
3. Step 3 or Step 4 (fix bracket) — based on findings
4. Step 6 (verify)

---

*Plan updated: 2026-04-10. Narrowed scope from full relaxation port to brent-specific fix based on user input that uniform/const_j work.*
