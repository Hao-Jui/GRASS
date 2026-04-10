# GRASS Codebase Reorganization Plan

## Requirements Summary

Reorganize GRASS from a 3-layer structure (`core/`, `theory/`, `tool/`) with a
monolithic `para_mod` (28 dependents, 285 lines mixing constants, grid config,
field arrays, targets, and subroutines) into a clean 8-layer DAG with no
upward dependencies.

**Goals:**
- Eliminate `para_mod` god-module by splitting into focused modules
- Enforce `numerics/` has zero physics dependencies
- Clear "where does new code go?" answer for every module
- Each phase compiles and runs before proceeding to the next

---

## Target Directory Layout

```
src/
├── main.f90
│
├── foundation/                    ← no project dependencies
│   ├── precision_mod.f90               wp = real64 (moved from core/)
│   └── constants_mod.f90               G, C, MSUN, MB, PI, KAPPA, ... (NEW)
│
├── numerics/                      ← depends only on foundation/
│   ├── toolkit_mod.f90                 interp, integrate_profiles
│   ├── ad_mod.f90                      automatic differentiation
│   ├── cheb_mod.f90                    Chebyshev polynomials
│   ├── spectral_hub_mod.f90            Gauss-Lobatto, diff matrices
│   ├── spline_mod.f90                  cubic spline
│   ├── brent_mod.f90                   root finding
│   └── nag_compat_mod.f90              adaptive quadrature
│
├── grid/                          ← depends on foundation/ + numerics/
│   ├── grid_config_mod.f90             SDIV, MDIV, LMAX, s_pwr, ... (NEW)
│   ├── grid_mod.f90                    s_gp, mu, derivative matrices
│   └── regrid_mod.f90                  adaptive remeshing
│
├── eos/                           ← depends on foundation/ + numerics/
│   └── eos_mod.f90                     EOS table + all channels
│
├── fields/                        ← depends on grid/ + eos/ + foundation/
│   └── fields_mod.f90                  all 2D/1D arrays, bulk scalars,
│                                       targets, alloc/dealloc, init (NEW)
│
├── theory/                        ← depends on fields/ + numerics/
│   ├── rotation_solver_mod.f90
│   ├── spin_integration_mod.f90
│   ├── spin_relaxation_mod.f90
│   ├── relaxation_mod.f90
│   ├── spin_updates_mod.f90
│   ├── spin_derivatives_mod.f90
│   ├── spin_workspace_mod.f90
│   └── rotational_law_mod.f90
│
├── solver/                        ← depends on theory/ + fields/
│   ├── shoot_mod.f90
│   ├── shoot_solver_2d_mod.f90
│   ├── shoot_solver_1d_hc_mod.f90
│   ├── shoot_solver_1d_r_ratio_mod.f90
│   ├── MRcurve_mod.f90
│   ├── starting_model_mod.f90
│   └── scalar_burning_mod.f90
│
├── diagnostics/                   ← depends on fields/ + numerics/
│   ├── analysis_mod.f90
│   ├── constraint_mod.f90
│   ├── ope_eq_mod.f90
│   ├── donutization_mod.f90
│   └── exporter_mod.f90
│
└── compat/                        ← legacy stubs
    ├── miscellaneous_mod.f90
    ├── set_disk.f90
    └── sphere_mod.f90
```

---

## para_mod Split Map

Current `para_mod` (src/core/para_mod.f90, 285 lines) → 3 new modules:

### constants_mod (NEW → foundation/constants_mod.f90)
Extracted from para_mod lines 155-173:
```
C, G, MSUN, MB, PI, ACCURACY, TOV_RMIN, KAPPA, KSCALE,
E_SURFACE, P_SURFACE, RHO_UNI, PRS_UNI, F_UNI, HBAR, L_UNI,
N_SAT, SCALARTON
```
Plus theory enum constants (lines 5-7): `THEORY_GR, THEORY_ST`
Plus mode/task enums (lines 11-17): `MODE_REGRID, MODE_DEFAULT, shoot, MRbuild, OneModel`
Plus shooting enums (lines 28-29): `SHOOT_FIX1_HC, SHOOT_FIX1_RP, SHOOT_2D`

### grid_config_mod (NEW → grid/grid_config_mod.f90)
Extracted from para_mod lines 34-87:
```
res, s_pwr, SDIV, MDIV, LMAX, SMAX, RDIV, DS, DM, S_E
s_gp(:), mu(:), sin_theta(:), D_mu(:,:), D_mu_t(:,:), D2_mu(:,:), w_mu(:)
angular_collocation, COLLOCATION_UNI, COLLOCATION_LEG, COLLOCATION_CHEB
```
Note: `s_gp`, `mu`, and derivative matrices are grid arrays populated by
`grid_mod.make_grid`. They move here because they define the computational
domain, not the physical solution.

### fields_mod (NEW → fields/fields_mod.f90)
Everything remaining from para_mod:
- **Config/targets** (lines 20-64): `solver_type`, `eos_file`, `M_goal`, `Mb_goal`,
  `J_goal`, `chi_goal`, `omc_goal`, `B_goal`, `mphi_goal`, `A_diff`, rotation-law
  params, `FIX1`, `FIX2`, `shooting`, `output`, `timing`, `run_mode`, `run_task`,
  `active_theory`
- **EOS state** (lines 66-74): `phase_transition`, `num_tab`, `n_PT`, `p_at_PT`,
  `log_p`, `log_e`, `log_h`, `log_n0`, `p_center`, `h_center`, `e_center`,
  `enthalpy_min`
- **Fluid arrays** (lines 98-101): `pressure`, `enthalpy`, `velocity_sq`, `energy`,
  `omg`, `F_j`, `v_plus`, `v_minus`, `V_rr_p`, `V_rr_m`, `sound_speed`
- **Metric arrays** (line 104): `gama`, `rho`, `ww`, `alpha`, `sphi`
- **Scalar field config** (lines 107-118): `has_scalar`, `B_coup`, `mphi_r`, etc.
- **Bulk properties** (lines 121-143): `mass`, `mass_0`, `chi`, `Omega_e`, etc.
- **Disk compat** (lines 88-96): `disk_present`, `edge_in`, etc.
- **Legendre weights** (line 148): `P_2n`, `P1_2n_1`, `sin_2n_1_theta`
- **Subroutines**: `initialize_theory`, `apply_gr_defaults`, `apply_st_defaults`,
  `allocate_fields`, `deallocate_fields`, `to_lower_str`

---

## Dependency DAG (enforced)

```
foundation  ──►  numerics  ──►  grid   ──►  fields  ──►  theory  ──►  solver
                            ──►  eos   ──►           ──►  diagnostics
```
No upward arrows. No lateral arrows between same-level siblings except:
- `eos/` may use `numerics/` (spline, toolkit)
- `diagnostics/` may use `eos/` (for eos_eval in analysis)
- `theory/` may use `eos/` (spin_updates needs eos_eval)
- `solver/` may use `diagnostics/` (shoot_mod calls analysis)

---

## Implementation Phases

### Phase 1: Extract constants_mod from para_mod
**Files created:** `src/core/constants_mod.f90` (stays in core/ for now)
**Files modified:** `src/core/para_mod.f90`, `Makefile`

Steps:
1. Create `constants_mod.f90` with all `parameter` physical constants (lines
   155-173) and enum constants (lines 5-7, 11-17, 28-29)
2. In `para_mod`, replace extracted parameters with `use constants_mod`
   and re-export via `public` (backward compatible — no downstream changes)
3. Add `constants_mod.f90` to Makefile `SOURCES` and dependency rules
4. **Verify:** `make clean && make -j8 && ./build/bin/a.out` with a quick
   OneModel run

**Risk:** Low — para_mod re-exports everything, so no downstream module changes.

### Phase 2: Extract grid_config_mod from para_mod
**Files created:** `src/core/grid_config_mod.f90` (stays in core/ for now)
**Files modified:** `src/core/para_mod.f90`, `Makefile`

Steps:
1. Create `grid_config_mod.f90` with grid parameters (`res`, `s_pwr`, `SDIV`,
   `MDIV`, `LMAX`, `SMAX`, `RDIV`, `DS`, `DM`, `S_E`, collocation enums)
   and grid arrays (`s_gp`, `mu`, `sin_theta`, `D_mu`, `D_mu_t`, `D2_mu`, `w_mu`)
   plus their allocate/deallocate
2. In `para_mod`, replace with `use grid_config_mod` + re-export
3. Update Makefile
4. **Verify:** `make clean && make -j8 && ./build/bin/a.out`

**Risk:** Low — same re-export pattern.

### Phase 3: Rename para_mod → fields_mod, update all consumers
**Files modified:** Every module that has `use para_mod`

This is the big phase. After Phases 1-2, `para_mod` contains only fields,
targets, EOS state, bulk properties, and subroutines. Rename it to
`fields_mod`.

Steps:
1. Rename `para_mod` → `fields_mod` in the module declaration
2. `fields_mod` keeps `use constants_mod` and `use grid_config_mod` internally
   but does NOT re-export them
3. Update every file's `use para_mod` → appropriate combination of:
   - `use constants_mod, only: ...` (for C, G, MSUN, PI, KAPPA, etc.)
   - `use grid_config_mod, only: ...` (for SDIV, MDIV, s_gp, mu, etc.)
   - `use fields_mod, only: ...` (for gama, rho, energy, mass, h_center, etc.)
4. Update Makefile dependency rules for all affected objects
5. **Verify:** `make clean && make -j8 && ./build/bin/a.out`

**Estimated changes:** ~28 files need `use` statement updates. Each file's
`only:` list already tells us exactly which symbols it needs — mechanical
but tedious.

**Risk:** Medium — many files touched. Mitigated by: each file's existing
`only:` clause makes the mapping unambiguous.

### Phase 4: Move files to new directories
**Files moved:** All source files to their target directories
**Files modified:** `Makefile` (path updates throughout)

Steps:
1. Create directories: `foundation/`, `numerics/`, `grid/`, `eos/`, `fields/`,
   `solver/`, `diagnostics/`, `compat/`
2. Move files per the layout table above (git mv for history)
3. Update all paths in Makefile (`SRC_CORE`, `SRC_THEORY`, `SRC_TOOL` →
   new per-directory variables)
4. Update any hardcoded paths in source (check `include` statements)
5. **Verify:** `make clean && make -j8 && ./build/bin/a.out`
6. Update `.codesight/` docs to reflect new layout

**Risk:** Low — no logic changes, only paths. But Makefile rewrite is the
most labor-intensive part.

---

## toolkit_mod dependency on para_mod

Currently `toolkit_mod` uses `para_mod` for `SDIV` and grid arrays in
`integrate_profiles`. After the split:
- `toolkit_mod` will `use grid_config_mod, only: SDIV` (grid config, not fields)
- This keeps `numerics/toolkit_mod` free of physics dependencies
- `grid_config_mod` is in `grid/` layer, which is above `numerics/` — this
  creates a **cross-layer dependency** (numerics → grid)

**Resolution options:**
1. **Pass SDIV as argument** to `integrate_profiles` instead of importing it
   → cleanest, but changes the call signature everywhere
2. **Move grid_config_mod to foundation/** → it has no dependencies beyond
   precision_mod, so this is valid architecturally
3. **Accept the dependency** — toolkit is a utility; knowing grid dimensions
   is reasonable

**Recommendation:** Option 2 — move `grid_config_mod` to `foundation/`.
Grid dimensions are foundational configuration, not computed grid data.
Separate `grid_mod` (which computes `s_gp`, `mu`, derivatives) stays in `grid/`.

Revised: `grid_config_mod` → `foundation/grid_config_mod.f90`
(parameters + dimension constants only; allocatable arrays for `s_gp`, `mu`
remain here since they're allocated early and used everywhere like config).

---

## Acceptance Criteria

- [ ] `make clean && make -j8` succeeds after each phase
- [ ] `./build/bin/a.out` produces identical output for a reference OneModel run
  after each phase (binary-identical `.dat` output)
- [ ] No module has `use para_mod` remaining after Phase 3
- [ ] `numerics/` modules depend only on `foundation/` (verify with `grep 'use '`)
- [ ] No circular dependencies in the `use` graph
- [ ] Every source file is in exactly one of the 8 target directories
- [ ] Makefile dependency rules match actual `use` statements
- [ ] `.codesight/CODESIGHT.md` and `graph.md` updated to reflect new layout

---

## Risks and Mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| Phase 3 breaks compilation | High | Each file's `only:` clause maps symbols unambiguously; do one file at a time |
| Makefile dependency errors | Medium | After Phase 4, run `make clean && make -j8`; missing deps show as compile errors |
| Binary output changes | High | Compare `.dat` output before/after each phase with `diff` |
| Merge conflicts with dev work | Medium | Complete reorganization on a feature branch; rebase before merge |
| `toolkit_mod` → `grid_config_mod` cross-layer | Low | Resolved by placing `grid_config_mod` in `foundation/` |

---

## Verification Steps

After each phase:
1. `make clean && make -j8` — compiles without errors
2. Run reference case: `OneModel` with MPA1 EOS, uniform rotation, Mb=1.8
3. `diff` output files against pre-reorganization baseline
4. `grep -r 'use para_mod' src/` — should show decreasing count (0 after Phase 3)
5. After Phase 4: `grep -rn 'use ' src/ | grep -v 'iso_fortran_env\|precision_mod'`
   → verify all `use` targets exist in new paths
