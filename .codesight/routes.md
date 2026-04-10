# routes.md — Solver Entry Points

> GRASS has no HTTP routes. This file documents **task entry points**: the
> three high-level solver modes and all secondary entry-point subroutines,
> including their required inputs (namelist parameters), outputs, and call
> chains.

---

## Task Modes

### 1. `shoot_v2` — Newton Shooting (Single Model)
**File**: `src/solver/shoot_mod.f90` → `subroutine shoot_v2`

**Purpose**: Find a rotating neutron-star solution that hits user-specified
bulk-property targets via Newton's method on the outer loop.

**Required Inputs** (from `fields_mod` / namelist):
- One of: `Mb_goal`, `M_goal`, `J_goal`, `chi_goal`, `omc_goal`
- `eos_file` — EOS table path (string, no extension)
- `solver_type` — `"uniform"` | `"const_j"` | `"uryu"`
- `r_ratio` or `h_center` — initial seed (1D mode)

**Sub-modes** (set internally):
| Mode constant | Meaning |
|---|---|
| `SHOOT_FIX1_HC` | 1D Newton fixing `h_center` |
| `SHOOT_FIX1_RP` | 1D Newton fixing `r_ratio` |
| `SHOOT_2D` | 2D Broyden Newton on `(r_e, rho0)` |

**Outputs** (written by `exporter_mod` / `analysis_mod`):
- Console: converged `mass`, `mass_0`, `r_e`, `Omega_c`, `chi`
- `Cont/{EOS}_J{spin}_Mb{mass}_...dat` — solution file

**Call chain**:
```
shoot_v2
  ├─ initialize_starting_model()       [starting_model_mod]
  ├─ Newton loop
  │   ├─ rotation_solver()             [rotation_solver_mod]
  │   ├─ mass_radius()                 [analysis_mod]
  │   └─ broyden_update() / clamp_step()  [shoot_solver_*_mod]
  └─ write_profiles()                  [exporter_mod]
```

---

### 2. `MRcurve` — Mass–Radius Sequence
**File**: `src/solver/MRcurve_mod.f90` → `subroutine MRcurve`

**Purpose**: Sweep `Mb_goal` over a range to build a full M–R curve.

**Required Inputs**:
- `Mb_goal` — starting baryon mass
- `Mb_step` — increment (internal default)
- `eos_file`, `solver_type` (same as shoot_v2)

**Outputs**:
- `properties.dat` (or `properties_st.dat`) — table of `(Mb, M, R, chi, Love2, ...)`
  per converged solution
- One solution file per point in `Cont/`

**Call chain**:
```
MRcurve
  └─ [loop over Mb_goal]
      └─ rotation_solver()             [rotation_solver_mod]
          └─ (same as shoot_v2 inner call)
```

---

### 3. `initialize_starting_model` — Static Seed
**File**: `src/solver/starting_model_mod.f90` →
`subroutine initialize_starting_model`

**Purpose**: Generate a spherically-symmetric (TOV) initial guess before
Newton iteration begins.

**Required Inputs**:
- `h_center` — central enthalpy
- `eos_file`

**Outputs**: populates all 2D field arrays in `fields_mod` with a static seed.

---

## Secondary Entry Points

### `rotation_solver` — Inner Iteration Loop
**File**: `src/theory/rotation_solver_mod.f90` →
`subroutine rotation_solver`

**Called by**: `shoot_v2`, `MRcurve`, `initialize_starting_model`

**Purpose**: Iterate to self-consistent rotating equilibrium for a fixed
`(r_e, rho0)` pair. Not called directly by users.

**Parameters that affect behaviour**:
- `ACCURACY` (1e-5) — convergence threshold for field dif
- `angular_collocation` — grid type (1/2/3)
- Relaxation thresholds in `spin_relaxation_mod`

---

### Post-Processing Entry Points

#### `auto2D.m` (MATLAB)
**Input**: `Map/1dprofile_{EOS}_1001_B2.6*.dat`
**Output**: 2D autoencoder latent representation Z1, Z2 of scalar profiles;
correlation with morphology diagnostics; PDF figures.

#### `donut.m` (MATLAB)
**Input**: `Map/1dprofile_{EOS}_1001_B2.5*.dat`
**Output**: 1D order parameter D for donutization; correlation with
`tail_power` and `width`; PDF figure.

#### `sub.sh`
**Effect**: `make clean && make -j8 && ./build/bin/a.out`
One-shot clean build and run.

#### `sweep_b_goal.sh [START] [END] [STEP]`
**Effect**: parametric sweep of `B_goal` over `[START, END]` in steps of
`STEP`; each run logged to `logs/b_goal_sweep/B_{value:05d}.log`.
Restores original `B_goal` on exit via `trap`.
