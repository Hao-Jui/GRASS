# middleware.md — Solver Pipeline Stages

> GRASS has no HTTP middleware. This file documents the **solver pipeline**:
> the sequence of processing stages, convergence guards, diagnostic hooks,
> and error handlers that every solution must pass through — the functional
> equivalent of middleware in a web stack.

---

## Pipeline Overview

```
main.f90
  │
  ├─ [Stage 0]  initialize_theory()     — theory selection & field allocation
  ├─ [Stage 1]  loadEos()               — EOS table load & phase-transition detect
  ├─ [Stage 2]  make_grid()             — collocation grid + derivative matrices
  │
  └─ [Task dispatch]
       ├─ shoot_v2          ─┐
       ├─ MRcurve           ─┼─► [Stage 3–7] rotation_solver inner loop
       └─ initialize_model  ─┘
             │
             └─ [Stage 8]  hamiltonian()  — post-solution constraint check
```

---

## Stage 0 — Theory Initialization
**Subroutine**: `initialize_theory()` in `fields_mod`
**Role**: select GR vs. ST mode, compute grid sizing constants (DS, DM),
allocate all 2D/1D field arrays.

**Guards**:
- Calls `apply_gr_defaults()` or `apply_st_defaults()` depending on
  `active_theory`
- Scalar-field arrays (`sphi`) allocated in both modes; zeroed in GR

---

## Stage 1 — EOS Load & Validation
**Subroutine**: `loadEos()` in `eos_mod`
**Role**: read `eos/{eos_file}.dat`, convert to log-space, detect phase
transitions, build cubic spline coefficients.

**Phase-transition handling**:
- Scans table for non-monotone pressure; records jump indices in `p_at_PT`
- Always builds cubic spline tables for all 7 interpolation directions:
  - If `n_PT > 0`: uses `build_spline_segmented` (independent splines per smooth segment, linear at PT gaps)
  - Otherwise: uses `build_spline` (single continuous not-a-knot spline)
- Sets `use_spline = .true.`; all public EOS API functions switch to spline eval (Horner, ~8 FLOPs vs ~50 for barycentric)

**Error guards**:
- `STOP` if file cannot be opened
- `STOP` if `num_tab < 4` (too few points for interpolation)
- Prints warning if phase-transition count `n_PT > 0` and
  `phase_transition = .false.`

---

## Stage 2 — Grid Construction
**Subroutine**: `make_grid()` in `grid_mod`
**Role**: fill `s_gp` (radial), `mu` (angular), derivative matrices
`D_mu`, `D2_mu`, `D_mu_t`, quadrature weights `w_mu`.

**Collocation types** (set via `angular_collocation`):
| Value | Grid | Called via |
|---|---|---|
| 1 | Uniform | Direct spacing |
| 2 | Legendre Gauss-Lobatto | `gauss_lobatto(MDIV, mu, w_mu)` |
| 3 | Chebyshev extrema | `chebyshev_lobatto_points(MDIV, mu)` |

After grid construction, `barycentric_diff_matrices(mu, D_mu, D2_mu)` is
called to build the spectral derivative operators.

---

## Stage 3 — Static Seed (TOV Initialization)
**Subroutine**: `initialize_starting_model()` in `starting_model_mod`
**Role**: set all field arrays to a spherically-symmetric (TOV) initial
guess using `h_center`.

**Called by**: `shoot_v2` once before the Newton loop begins.

---

## Stage 4 — Outer Newton Loop (shoot_v2 only)
**Subroutine**: `shoot_v2` in `shoot_mod`
**Role**: iteratively adjust `(r_e, rho0)` (or `h_center` / `r_ratio` in
1D mode) until bulk-property targets are met.

**Newton modes**:
| Mode | Fixed variable | Solved for |
|---|---|---|
| `SHOOT_FIX1_HC` | `h_center` | 1D Newton on one target |
| `SHOOT_FIX1_RP` | `r_ratio` | 1D Newton on one target |
| `SHOOT_2D` | — | 2D Broyden on `(r_e, rho0)` |

**Broyden update** (`shoot_solver_2d_mod`):
- Rank-1 Broyden formula: J_new = J_old + (ΔF − J·Δx)·Δxᵀ / ‖Δx‖²
- `clamp_step` limits |Δx| to prevent divergence

**Termination**: when residual |F| < Newton tolerance (typically 1e-6).

---

## Stage 5 — Rotation Solver Inner Loop
**Subroutine**: `rotation_solver()` in `rotation_solver_mod`
**Role**: iterate fields (γ, ρ, ω, α, φ) to self-consistent equilibrium
for a given `(r_e, rho0)`.

**Loop structure** (per iteration):
1. `rescale_metric(r_e_new_sq)` — normalize metric potentials to physical scale
2. `update_equatorial_radius(...)` — update r_e estimate
3. `update_angular_velocity(...)` — compute Ω from rotation law
4. `update_eos_and_velocity(...)` — interpolate P, h, ε, v from EOS
5. `get_all_targets(...)` — compute Green's function multipole targets
   for ρ, γ, ω, φ
6. `relaxation(...)` — **4-stage adaptive relaxation** (see Stage 6)
7. `update_alpha_lapse()` — solve elliptic equation for lapse α
8. `mass_radius()` — compute bulk diagnostics
9. Check convergence: `dif < ACCURACY` (default 1e-5)

---

## Stage 6 — 4-Stage Adaptive Relaxation
**Subroutine**: `relaxation()` in `spin_relaxation_mod`
**Role**: apply one relaxation update to field arrays using the most
appropriate acceleration method for the current convergence level.

### Stage 6a — Picard Damping (`dif ≥ 0.5`)
```
x_new = (1 − W_PICARD) · x_old + W_PICARD · x_target
W_PICARD = 0.7
```
**Guard**: stall detection — if `dif` does not decrease for
`N_PICARD_STALL` consecutive iterations, reduce `W_PICARD` adaptively.

### Stage 6b — Chebyshev Acceleration (`0.1 ≤ dif < 0.5`)
Adaptive weight derived from spectral-radius estimate `ρ_spec`:
```
w_cheb = 2 / (2 − ρ_spec · (w_prev + w_prev_inv))
```
**Guards**:
- Requires `N_CHEB = 3` warm-up iterations to estimate ρ_spec
- Falls back to Picard if Chebyshev weight diverges
- `rho_spec_est` updated with exponential smoothing

### Stage 6c — Anderson Acceleration (`0.0001 ≤ dif < 0.1`)
Rank-M=3 Anderson mixing (`anderson_accel_optimized` in `relaxation_mod`):
- Maintains history of last M=3 field vectors and residuals
- Solves a small M×M least-squares problem per iteration
**Guards**:
- `N_ANDERSON = 5` warm-up iterations required
- `ANDERSON_COOLDOWN = 30` iteration cooldown after failure
- Falls back to Picard on failure

### Stage 6d — Aitken δ² Extrapolation (`dif < 1e-5`)
Quadratic extrapolation on field values using three successive iterates:
```
x_new = x_curr − (x_curr − x_prev)² / (x_curr − 2·x_prev + x_prev2)
```
**Guards**:
- Activated only after `N_AITKEN = 2` consistent-trend iterations
- `AITKEN_COOLDOWN = 3` iteration cooldown after failure
- `n_rho_locked` counter gates re-activation

---

## Stage 7 — Field Update Application
**Subroutine**: `spin_updates_mod` routines
**Role**: write updated field values back to `fields_mod` arrays after
relaxation; update rotation law profile via `rotational_law_mod`.

**Rotation law dispatch**:
| `solver_type` | Subroutine |
|---|---|
| `"uniform"` | Direct assignment Ω = Ωc |
| `"const_j"` | `apply_rotation_law` with `A_diff` |
| `"uryu"` | `find_rotation_root` via `brent_root` |

---

## Stage 8 — Post-Solution Diagnostics

### Hamiltonian Constraint Check
**Subroutine**: `hamiltonian(hamL2)` in `constrain_mod`
**Role**: compute and report L2 norm of the Hamiltonian constraint
(should be ≪ 1 for a well-converged solution).
**Output**: printed to console; not a hard stop.

**Constraint form**:
```
ham = R_ricci + KK_term − 16π ρ_H a_coup⁴ − 2(|∇φ|² + 2V_φ)
```
where the scalar term is zero in GR mode.

### Mass & Radius Diagnostics
**Subroutine**: `mass_radius()` in `analysis_mod`
**Role**: compute all bulk properties (M, Mb, J, χ, I, k₂, r_e, …) and
write them to `fields_mod` scalar outputs.

**Quadrature method**: multi-column `integrate_profiles` for angular
integrals; radial summation over DS steps.

### Solution Export
**Subroutine**: `write_profiles` in `exporter_mod`
**Role**: write 2D field arrays and 1D profiles to `Cont/` as `.dat` files.
Called after convergence; not called on failed iterations.

---

## Error Handling Summary

| Condition | Response | Location |
|---|---|---|
| EOS file not found | `STOP` with message | `eos_mod::loadEos` |
| EOS table too small | `STOP` with message | `eos_mod::loadEos` |
| Chebyshev weight divergence | Fall back to Picard | `spin_relaxation_mod` |
| Anderson failure | Cooldown + Picard fallback | `spin_relaxation_mod` |
| Aitken divergence | Cooldown + reset | `relaxation_mod::aitken_delta2` |
| Newton non-convergence | Print warning, continue | `shoot_mod::shoot_v2` |
| Hamiltonian constraint large | Print warning (hamL2) | `constrain_mod::hamiltonian` |
| BLAS/LAPACK absent | Linker error at build | Makefile |
