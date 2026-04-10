# graph.md — Import Dependency Graph

> This file maps every `USE` dependency in the codebase and identifies
> which modules are most dangerous to change (highest fan-in / most
> downstream dependents).

---

## Dependency Graph

Format: `MODULE -> [modules it USEs]`

### Foundation Layer

| Module | Uses |
|---|---|
| `precision_mod` | *(none -- root node)* |
| `constants_mod` | `precision_mod` |
| `grid_config_mod` | `precision_mod` |

### Fields Layer

| Module | Uses |
|---|---|
| `fields_mod` | `precision_mod`, `constants_mod`, `grid_config_mod` |

### Numerics Layer

| Module | Uses |
|---|---|
| `ad_mod` | `precision_mod` |
| `brent_mod` | `precision_mod` |
| `cheb_mod` | `precision_mod` |
| `spectral_hub_mod` | `precision_mod`, `cheb_mod` |
| `nag_compat_mod` | `precision_mod` |
| `spline_mod` | `precision_mod` |
| `toolkit_mod` | `precision_mod`, `ad_mod`, `spectral_hub_mod`, `grid_config_mod`, `fields_mod` |

### EOS & Grid Layer

| Module | Uses |
|---|---|
| `eos_mod` | `ad_mod`, `toolkit_mod`, `spline_mod`, `fields_mod` |
| `grid_mod` | `toolkit_mod`, `grid_config_mod`, `fields_mod` |

### Diagnostics Layer

| Module | Uses |
|---|---|
| `ope_eq_mod` | `grid_config_mod`, `fields_mod` |
| `constraint_mod` | `grid_config_mod`, `fields_mod`, `ope_eq_mod` |
| `donutization_mod` | `constants_mod` |
| `exporter_mod` | `eos_mod`, `grid_config_mod`, `fields_mod` |
| `analysis_mod` | `ad_mod`, `cheb_mod`, `nag_compat_mod`, `toolkit_mod`, `exporter_mod`, `miscellaneous_mod`, `fields_mod`, `eos_mod` |

### Compat Layer

| Module | Uses |
|---|---|
| `miscellaneous_mod` | `nag_compat_mod`, `fields_mod` |
| `set_disk` | `toolkit_mod`, `fields_mod`, `eos_mod` |
| `sphere_mod` | `toolkit_mod`, `fields_mod`, `eos_mod`, `set_disk` |
| `regrid_mod` | `fields_mod`, `grid_mod` |

### Theory Layer

| Module | Uses |
|---|---|
| `relaxation_mod` | `grid_config_mod` |
| `rotational_law_mod` | `brent_mod`, `fields_mod`, `grid_config_mod` |
| `spin_derivatives_mod` | `fields_mod`, `grid_config_mod` |
| `spin_updates_mod` | `toolkit_mod`, `fields_mod`, `rotational_law_mod`, `brent_mod`, `eos_mod` |
| `spin_workspace_mod` | `nag_compat_mod`, `fields_mod`, `grid_config_mod` |
| `spin_integration_mod` | `spin_workspace_mod`, `spin_derivatives_mod`, `toolkit_mod`, `exporter_mod`, `fields_mod`, `eos_mod`, `donutization_mod` |
| `spin_relaxation_mod` | `spin_workspace_mod`, `fields_mod`, `relaxation_mod` |
| `rotation_solver_mod` | `constants_mod`, `grid_config_mod`, `fields_mod`, `spin_workspace_mod`, `spin_integration_mod`, `spin_relaxation_mod`, `spin_updates_mod`, `analysis_mod` |

### Solver Layer

| Module | Uses |
|---|---|
| `shoot_solver_2d_mod` | `fields_mod`, `rotation_solver_mod`, `eos_mod`, `analysis_mod` |
| `shoot_solver_1d_hc_mod` | `fields_mod`, `rotation_solver_mod`, `eos_mod`, `analysis_mod` |
| `shoot_solver_1d_r_ratio_mod` | `shoot_solver_1d_hc_mod`, `fields_mod`, `rotation_solver_mod`, `eos_mod`, `analysis_mod` |
| `scalar_burning_mod` | `fields_mod`, `rotation_solver_mod` |
| `starting_model_mod` | `rotation_solver_mod`, `miscellaneous_mod`, `scalar_burning_mod`, `eos_mod`, `analysis_mod`, `regrid_mod`, `sphere_mod` |
| `shoot_mod` | `donutization_mod`, `spin_workspace_mod`, `starting_model_mod`, `eos_mod`, `analysis_mod`, `shoot_solver_2d_mod`, `shoot_solver_1d_hc_mod`, `shoot_solver_1d_r_ratio_mod` |
| `MRcurve_mod` | `donutization_mod`, `spin_workspace_mod`, `rotation_solver_mod`, `starting_model_mod`, `eos_mod`, `analysis_mod` |

### Main Program

| File | Uses |
|---|---|
| `src/main.f90` | `precision_mod`, `constants_mod`, `grid_config_mod`, `fields_mod`, `eos_mod`, `grid_mod`, `toolkit_mod`, `constraint_mod`, `shoot_mod`, `starting_model_mod`, `MRcurve_mod` |

---

## Focused Usage Tree: `ad_mod`

This is the active reachability tree for the built module in
`src/numerics/ad_mod.f90`.

```text
src/numerics/ad_mod.f90
+- src/numerics/toolkit_mod.f90
|  +- interp_linear_segment_dual
|  +- interp_dual
|  +- interp_pt_dual
+- src/eos/eos_mod.f90
|  +- module-level dual/dual_const/operator support
|  +- eos_eval_dual(eos_e2p, ...)
+- src/diagnostics/analysis_mod.f90
   +- solution_properties
   |  +- interp_dual(..., dual_var(...), ...)
   |  +- eos_eval_dual(eos_e2p, energy_dual)
   |  +- moment_inertia()
   |     +- deriv
   |        +- dual_var(s_h)
   |        +- interp_dual(...)
   +- deriv

src/main.f90
+- shoot_v2 -> shoot_mod -> solution_properties
+- MRcurve -> MRcurve_mod -> solution_properties
+- initialize_starting_model -> starting_model_mod -> solution_properties
```

Notes:

- `src/numerics/ad_mod.f90` is live in the current Makefile build.
- The AD-dependent helpers are not dead-end helpers; they are reached from
  `solution_properties` and `deriv` on active solver paths.

---

## Reverse Dependency (Fan-In) Table

Sorted by number of direct dependents -- highest fan-in = highest breakage
risk if changed.

| Rank | Module | # Direct Dependents | Dependent modules |
|---|---|---|---|
| 1 | **`fields_mod`** | ~20 | Most core, theory, solver, diagnostics, and compat modules |
| 2 | **`precision_mod`** | 15 | Nearly every module (via `wp`) |
| 3 | **`grid_config_mod`** | ~10 | grid_mod, toolkit_mod, theory modules, diagnostics, rotation_solver |
| 4 | **`constants_mod`** | ~5 | main, donutization, rotation_solver, fields_mod |
| 5 | **`toolkit_mod`** | 8 | `eos_mod`, `analysis_mod`, `constraint_mod`, `spin_integration_mod`, `spin_relaxation_mod`, `spin_derivatives_mod`, `set_disk`, `sphere_mod` |
| 6 | **`eos_mod`** | 5 | `starting_model_mod`, `sphere_mod`, `spin_updates_mod`, `scalar_burning_mod`, `analysis_mod` |
| 7 | **`rotation_solver_mod`** | 2 | `shoot_mod`, `MRcurve_mod` |
| 8 | **`analysis_mod`** | 3 | `shoot_mod`, `MRcurve_mod`, `rotation_solver_mod` |
| 9 | **`spectral_hub_mod`** | 2 | `grid_mod`, `toolkit_mod` |
| 10 | **`relaxation_mod`** | 1 | `spin_relaxation_mod` |

---

## Compilation Order (from Makefile)

The Makefile enforces this build ordering (lines 81-136):

```
precision_mod          (root)
  +- constants_mod
  +- grid_config_mod
  +- ad_mod
  +- brent_mod
  +- cheb_mod
  |   +- spectral_hub_mod
  +- nag_compat_mod
  +- spline_mod
  +- fields_mod (constants_mod, grid_config_mod)
  +- toolkit_mod (ad_mod, spectral_hub_mod, grid_config_mod, fields_mod)
  +- eos_mod (ad_mod, toolkit_mod, spline_mod, fields_mod)
  +- grid_mod (toolkit_mod, grid_config_mod, fields_mod)
  +- ope_eq_mod
  +- constraint_mod
  +- donutization_mod
  +- exporter_mod
  +- analysis_mod
  +- miscellaneous_mod
  +- set_disk
  +- sphere_mod
  +- regrid_mod
  +- relaxation_mod
  +- rotational_law_mod
  +- spin_derivatives_mod
  +- spin_updates_mod
  +- spin_workspace_mod
  +- spin_integration_mod
  +- spin_relaxation_mod
  +- rotation_solver_mod
  +- shoot_solver_2d_mod
  +- shoot_solver_1d_hc_mod
  +- shoot_solver_1d_r_ratio_mod
  +- scalar_burning_mod
  +- starting_model_mod
  +- shoot_mod
  +- MRcurve_mod
      +- main.f90  (program)
```

---

## Top-10 Highest Risk Files

Files whose modification is most likely to break other things:

| Risk | File | Why |
|---|---|---|
| :red_circle: Critical | `src/para_panel.f90` | ~20 dependents; all field arrays and bulk properties live here |
| :red_circle: Critical | `src/macros/precision_mod.f90` | `wp` type used everywhere; changing precision requires full recompile |
| :orange_circle: High | `src/macros/grid_config_mod.f90` | Grid dimensions and arrays used by ~10 modules |
| :orange_circle: High | `src/numerics/toolkit_mod.f90` | `interp`, `integrate_profiles` called throughout theory and solver layers |
| :orange_circle: High | `src/eos/eos_mod.f90` | EOS functions called inside tight loops; signature changes break multiple modules |
| :orange_circle: High | `src/theory/spin_relaxation_mod.f90` | 4-stage relaxation logic; changing thresholds affects convergence of all runs |
| :yellow_circle: Medium | `src/theory/rotation_solver_mod.f90` | Outer loop logic; used by both `shoot_mod` and `MRcurve_mod` |
| :yellow_circle: Medium | `src/diagnostics/analysis_mod.f90` | Bulk-property integrals; results feed Newton convergence checks |
| :yellow_circle: Medium | `src/numerics/spectral_hub_mod.f90` | Derivative matrices used by `grid_mod` and `toolkit_mod` |
| :green_circle: Lower | `src/solver/shoot_mod.f90` | High-level driver; Newton logic is isolated but changes affect task entry |

---

## Isolated / Safe-to-Change Modules

These modules have no dependents (leaf nodes) or very few:

| Module | Dependents | Notes |
|---|---|---|
| `donutization_mod` | 0 direct | Safe to extend without recompiling anything else |
| `brent_mod` | 1 | Only `rotational_law_mod` |
| `nag_compat_mod` | 1-2 | Quadrature wrapper; replace internals freely |
| `exporter_mod` | 2 | Output formatting; safe to extend |
| `miscellaneous_mod` | 1 | Print utilities only |
| `set_disk` | 0 active | Legacy compatibility shim |
| `scalar_burning_mod` | 1 | Only called from `starting_model_mod` |
| `relaxation_mod` | 1 | Only `spin_relaxation_mod` |
