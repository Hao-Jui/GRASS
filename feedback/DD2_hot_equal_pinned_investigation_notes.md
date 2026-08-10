# DD2_hot_equal_pinned in GRASS: TOV star, Uryu-law bugs, and the unresolved rotating-Ham-L2 mystery

Investigation notes from building a TOV star and then a differentially-rotating
star in GRASS from `DD2_hot_equal_pinned` (PyCompOSE's hot-merger-remnant
extension of the DD2 EOS). Covers everything found and fixed, in
chronological order, plus what remains open.

## Summary

- **Phase 1 (TOV/non-rotating star): fully working**, `Ham L2≈5e-3`.
- **Two real, confirmed GRASS bugs found and fixed** in the Uryu
  differential-rotation law's implementation (independent of any EOS
  choice — both also affect MPA1 in the same regime).
- **A confirmed, large, general GRASS finding**: a rotating star's vacuum-
  surface pressure/density cutoff being far above the standard near-vacuum
  convention devastates `Ham L2` — demonstrated on MPA1 (`1e-3 → 153`,
  ~150,000×, from raising its surface alone).
- **An unresolved mystery**: `DD2_hot_equal_pinned` gives `Ham L2≈1.05` for
  *any* differential rotation law, completely insensitive to fixing its own
  surface level and every other EOS-table property tested (9 independent
  hypotheses eliminated). Root cause not found — likely something in
  GRASS's shared rotation-solver internals reacting to this star's specific
  bulk structure, not the EOS table itself.

---

## Phase 1: TOV (non-rotating) star — SOLVED

**Problem**: `sphere_mod::TOV`'s radial integrator has no iteration cap, and
`DD2_hot_equal_pinned` (built from `PyCompOSE`'s `fit_appendixA.py`, a hot
merger-remnant profile) doesn't extend anywhere near GRASS's default vacuum
surface convention (`E_SURFACE=7.8 g/cm³`, `P_SURFACE=1.01e8 dyn/cm²` —
literally iron density). The table's own low-density floor sat ~13 orders
of magnitude higher in pressure at the time, so the integrator ran
effectively forever trying to reach an unreachable target.

**Fix (at the time)**: re-pointed `E_SURFACE`/`P_SURFACE` in
`src/para_panel.f90` at the table's own low-density floor instead of the
standard convention. (Superseded later — see "Absolute surface level"
below — but this got Phase 1 unblocked.)

**Also fixed**: `r_ratio` was being unconditionally clobbered to `0.7` by
`single_model()` in `src/core/starting_model_mod.f90` regardless of what was
configured, so `OneModel` always built a rotator, never a true TOV star.
Removed the hard-coded overwrite (and an earlier `solver_type`-based guess
in the `case default` branch of `initialize_starting_model`) so `r_ratio` is
honored as configured — `1.0` for Phase 1, a genuine rotator value for
Phase 2.

**Environment note**: the build links MKL's `libmkl_rt.so.2`, which lazily
dispatches to `libmkl_intel_thread.so.2` on first real MKL call — that
threading layer needs Intel's OpenMP runtime, not present in this shell.
Run with `MKL_THREADING_LAYER=GNU` to route through `libgomp` instead (no
source change needed).

**Result**: converges cleanly, `Ham L2≈5e-3` — actually *better* than
GRASS's own MPA1 reference test achieves for a rotating model.

---

## Phase 2: Uryu-law differential rotation — two real bugs found and fixed

### Bug 1: `brent_core: failed to bracket root` (NaN in `AA_h`)

**Symptom**: crashed when spinning up from the TOV restart to
`r_ratio=0.995` with the Uryu rotation law.

**Root cause**: `AA_h`'s defining exponent, `1/(uryu_p+uryu_q) = 1/4`, is
fractional. Its bracketed base — `(F_e·Fmax_h)·tmp1/tmp2` — goes negative
for the *entire band* `0 ≤ F_e < (lambda1/lambda2)^(1/uryu_q)·Fmax_h` (not
just `F_e<0`). Fortran's `real**noninteger` on a negative base silently
returns `NaN`. `find_omega_e`'s geometric bracket search starts right in
this near-spherical, `F_e≈0` regime (exactly where a gradual spin-up seeds
it from a TOV restart) and can't escape the NaN wall.

**Fix**: `src/theory/rotational_law_mod.f90` — added a shared `safe_F_e`
function that clamps `F_e` comfortably above the whole unsafe band before
calling `AA_h`/`BB_h`. Applied at **two** call sites:
- `diff_rotation_uryu` (the original crash site).
- `cache_uryu_ab` (a second, initially-missed call site using the raw,
  unclamped `F_equator_h` global — this was poisoning `cached_aa`/`cached_bb`
  with `NaN` even after the first fix, breaking the downstream peak-search
  `zbrent_rot`/`rotation_law_uryu` calls with a different crash signature,
  `Cont/rotation_law.dat`).

**Verified universal, not EOS-specific**: MPA1 hits the *identical* crash
at the same `r_ratio`/default shape parameters — this is a general Uryu-law
implementation bug, unrelated to `DD2_hot_equal_pinned`.

### Bug 2 (not a crash, a stop): `negative F_equator_h`

**Symptom**: after fixing Bug 1, spin-up to `r_ratio=0.995` with the
*default* Uryu shape (`lambda1=1.5, lambda2=0.3`) ran cleanly through the
whole solve, then hit `STOP negative F_equator_h; L120 in uryu` (a clean
stop, not a crash) — the found `Omega_e` came out less than the local
frame-dragging value.

**Root cause**: confirmed via diagnostic instrumentation — the found root
tracks smoothly for hundreds of outer `Fmax_h`-shooting iterations, then at
a specific `Fmax_h` threshold the root **jumps discontinuously** to a
spurious branch (`Omega_e` collapsing toward `0`). Two fix attempts (a
physically-motivated `Fmax_h` reseed; a bounded-step continuation method on
the shooting loop) both failed to help — the bounded-step attempt failed
*faster*, which rules out a numerical-overshoot explanation and confirms
this is a genuine branch discontinuity in the equation for extreme
differential rotation this close to spherical.

**Resolution**: not a bug fix — milder Uryu shape parameters
(`lambda1=1.1, lambda2=0.9` instead of `1.5/0.3`) avoid the discontinuity
entirely and converge cleanly. **Confirmed universal**: MPA1 hits the same
stop at the same `r_ratio` with the default shape too — this is inherent to
extreme differential rotation near `r_ratio=1` for *any* EOS, not specific
to our table.

---

## The `Ham L2` investigation: a large real effect found, but DD2's own problem remains unexplained

With the crash fixed and milder shape parameters avoiding the stop, DD2
rotates without crashing — but `Ham L2≈1.05-1.3`, 1-3 orders of magnitude
worse than MPA1 achieves under the identical solver/shape/rotation rate
(MPA1: `≈1e-3` at `r_ratio=0.995`, `≈0.09` at `r_ratio=0.7`).

### Hypotheses tested and eliminated (in order)

1. **Vacuum-surface placement within the table's own narrow range** — moved
   the cutoff from `e=2.4e3` to `e=1.9e3 g/cm³`: no change.
2. **Table resolution** — resampled 326→6000 rows (PCHIP, same curve): no
   change.
3. **Warm-vs-cold-crust construction** — reworked `fit_appendixA.py`'s
   low-density leg to prefer staying at the table's coldest tabulated `T`
   (using DD2's own physically self-consistent NSE crust) instead of
   warming via entropy-pinning: fixed a real physics issue (the warmed
   version was unphysically soft, `Gamma1` ~0.5-0.85× MPA1's), but `Ham L2`
   unchanged.
4. **A sharp local `Gamma1` discontinuity** at the `i_min` boundary
   (`T` jumping `2.32→5.67 MeV` in one table row) — fixed with a smooth
   log-linear `T` glide window; `Gamma1` spike reduced `9.37→2.15`, but
   `Ham L2` **unchanged**, and — checked directly in the 2D constraint
   field — the exact same dominant grid shell contributed the exact same
   magnitude before and after.
5. **First-law/Euler-relation thermodynamic consistency** — checked
   directly against the table's own `Q1,Q2,Q3,Q4,Q7` at and around the
   implicated density: smooth, no violation, actually *better* than a
   random baseline point elsewhere in the table.
6. **GRASS's own grid resolution** — doubled `SDIV` (`801→1601`): `Ham L2`
   got **worse** (`1.067→1.288`), ruling out simple under-resolution (a
   discretization error should improve with resolution, not degrade).
7. **Rotation-law-specific bug** — tried `const_j` (a completely different
   formulation, no shared code with Uryu): identical `Ham L2≈1.06`,
   confirming the issue has nothing to do with the Uryu implementation.

### The one hypothesis that WAS confirmed — on MPA1, not on DD2

8. **Absolute vacuum-surface pressure level** (not just "which row within
   DD2's own table"): raised **MPA1's own** surface cutoff to match DD2's
   absolute level (`e≈2000 g/cm³`, an in-table, non-extrapolated MPA1 row).
   Result: MPA1's `Ham L2` exploded from `~1e-3` to **`153`** — a
   ~150,000× degradation. This is a large, real, confirmed effect: an
   elevated vacuum-surface cutoff, independent of which EOS underlies it,
   badly destabilizes a rotating solve.

### Fixing DD2's own surface level — still didn't help

Given finding 8, extended `DD2_hot_equal_pinned` down to the standard
near-vacuum convention using a cold degenerate polytropic append in
`fit_appendixA.py` (matched continuously at the table's own floor,
`n0=1e-12 fm⁻³`):

- First attempt used `Gamma=1.584` (matched to *MPA1's* edge slope) — hit a
  catastrophic-cancellation bug (absolute `MeV/fm³` energy-density
  arithmetic subtracting two near-equal `O(mn·n0)` terms), then, once fixed
  numerically, introduced a **new** `Gamma1` discontinuity at the seam
  (`1.58→0.95`) because DD2's own edge slope (`~0.95`, confirming the
  "warm-crust softness" finding) doesn't match MPA1's.
- Fixed by matching `Gamma` to **DD2's own edge slope** instead (via the
  same 3-point one-sided log-log derivative as PyCompOSE's
  `Table.get_polytrope(nb_idx=0)`, `Gamma≈0.9495`) — confirmed smooth
  `Gamma1` continuity at the seam (`0.946→0.951`, no jump).
- Result with the surface now at `e≈8.8 g/cm³` (matching the standard
  convention almost exactly) **and** a smooth seam: **`Ham L2` still
  `≈1.05`, completely unchanged.**

### The puzzle this leaves

MPA1 is *highly* sensitive to surface level across the tested range
(`8.9→2000 g/cm³`: `1e-3→153`). DD2 shows **zero** sensitivity across the
*same* range (`8.8→2400 g/cm³`: always `~1.05-1.08`). If surface level were
DD2's dominant problem, some sensitivity should show up somewhere in that
range — it doesn't. This means DD2 has its own, apparently much larger,
independent problem that the surface-level effect (real, proven on MPA1)
doesn't explain and fixing the surface can't touch.

**Not yet investigated**: GRASS's own shared rotation-solver internals
(`update_alpha_potential`, the Green's-function summation in
`spin_integration_mod.f90`) — the remaining candidate, since every EOS-table
property tried (9 hypotheses total) has come back either negative or
already-fixed-without-effect.

---

## Files changed

### GRASS (`/u/tlam/GRASS`)

- `src/para_panel.f90` — `r_ratio` default, `E_SURFACE`/`P_SURFACE`
  (currently `8.82 g/cm³` / matching pressure — the near-standard,
  seam-smoothed table's own row), `eos_file` (currently
  `"DD2_hot_equal_pinned"`).
- `src/core/starting_model_mod.f90` — removed the `r_ratio` hard-coded
  overwrites (case-default seed guess, `single_model`'s `0.7` clobber);
  `MODE_REGRID` branch now correctly re-imposes a *target* `r_ratio`
  captured before `regrid_read` can clobber it, instead of the old
  hard-coded `e_center=0.74e15` override.
- `src/theory/rotational_law_mod.f90` — new shared `safe_F_e` function;
  applied in `diff_rotation_uryu` and `cache_uryu_ab`.
- `src/core/constraint_mod.f90` — enabled the (previously dead-code-gated)
  2D Hamiltonian-field dump to `Cont/hamiltonian.dat`, added `r_phys` and
  `enthalpy` columns to it (diagnostic only, still enabled).

### PyCompOSE (`/sakura/ptmp/tlam/PyCompOSE/examples/DD2_hot`)

- `fit_appendixA.py`:
  - Reverse-pass low-density construction re-ordered to prefer staying at
    the table's coldest tabulated `T` (using DD2's own physical NSE crust),
    falling back to warming only as a last resort (previously the reverse:
    entropy-pin/warm first).
  - Added a smooth `T`-glide transition window (`WIN=20` rows) approaching
    the `i_min` boundary, replacing an abrupt jump to the untouched,
    original Appendix-A-fit point.
  - Added a cold degenerate polytropic low-density extension (120 rows)
    below the table's own floor, matched continuously in `Q7` (not
    absolute energy density — avoids catastrophic cancellation) with
    `Gamma` fit to the table's own edge slope, extending the GRASS-format
    output down to `n0=1e-16 fm⁻³` (`e≈0.165 g/cm³`) so a standard
    near-vacuum surface cutoff is directly tabulated, no extrapolation
    needed.
  - `write_grass_table` call now writes the full extended (446-row) table.
- New diagnostic scripts (all still present, not part of the regular
  pipeline): `resample_grass_table.py`, `compare_eos_derivatives.py`,
  `check_first_law.py`.

## Current repo state

- `/u/tlam/GRASS`: git-dirty with the four source files above; `git diff
  --stat` shows the exact changes. `eos/DD2_hot_equal_pinned.dat` (446
  rows, untracked) and `eos/DD2_hot_equal.dat` (the older, pre-investigation
  table, untracked) both present.
- `para_panel.f90` last left at: `eos_file="DD2_hot_equal_pinned"`,
  `solver_type="uryu"`, `run_mode=MODE_REGRID`, `r_ratio=0.995`,
  `lambda1=1.1`, `lambda2=0.9`, `E_SURFACE≈8.82 g/cm³`/matching `P_SURFACE`
  — i.e. mid-experiment state from the last `Ham L2` test, not reset to a
  "resting" TOV configuration.
- `Res/res.rst` holds the last-built restart (the seam-smoothed-table TOV
  star at `r_ratio=1.0`) — regenerate via `solver_type="uniform"`,
  `run_mode=MODE_DEFAULT`, `r_ratio=1.0` before any further rotating test.
- Run with `MKL_THREADING_LAYER=GNU` from `/u/tlam/GRASS` (relative
  `./eos/`, `./Cont/`, `./Res/` paths require this cwd).
