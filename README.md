# GRASS
## General Relativistic Axisymmetric Spacetime Solver

GRASS computes rotating neutron-star equilibria in either General Relativity (GR)
or Scalar–Tensor (ST) gravity. The code integrates the field equations on a
compactified meridional grid, reads tabulated equations of state (EOS),
and outputs diagnostic profiles that characterise the equilibrium model.

The v1 release: 6139a9f
The lastest official release: 8a7b201

---

## Features

- Uniformly rotating, axisymmetric neutron-star configurations.
- Tabulated EOS support (pressure, energy density, enthalpy tables).
- Compactified grid in meridional coordinates with spectral-type
  interpolation.
- Batch integration helper (`integrate_profiles`) to reduce repeated NAG
  calls.
- Newton shooting solvers (1D and 2D) to match target bulk properties
  such as gravitational mass, baryon mass, spin parameter, and axis ratio.

---

## 4-Stage Relaxation Strategy

The rotation solver employs an adaptive four-stage relaxation method to efficiently converge the coupled metric and scalar-field equations. Each stage activates at different convergence regimes, balancing stability and speed.

### Stage 1: Picard Damping (dif ≥ 0.1)

**Regime:** Early nonlinear phase, large residuals
**Update rule:** `x_new = (1 - w) * x + w * target`, with w = 0.7
**Purpose:** Conservative damping prevents overshoot during the initial drift phase when the solution is far from equilibrium. The fixed weight w=0.7 ensures stable convergence regardless of field coupling strength.

### Stage 2: Chebyshev Acceleration (0.01 ≤ dif < 0.1)

**Regime:** Transition zone, metrics and scalar field approaching equilibrium
**Update rule:** Chebyshev 3-term recurrence `cheb_w = 1 / (1 - (ρ²/4) * cheb_w)`
**Spectral radius:** Estimated via EMA-smoothed dif ratio: `ρ_obs = dif_n / dif_{n-1}`
**Purpose:** Adaptive acceleration using the estimated spectral radius of the linearized system. As dif contracts monotonically, `ρ_obs` becomes stable, and `cheb_w` grows toward its fixed point, giving up to 2× speedup over Picard.

**Key parameter:** `PICARD_THRESH = 1.0E-1`, `CHEB_THRESH = 1.0E-3`

### Stage 3: Anderson Acceleration (dif < 0.01, non-exponential)

**Regime:** Near-convergence, rapid residual reduction
**Method:** Rank-M Anderson (M_HIST = 3), applied per-field (rho, gama, ww, sphi)
**Purpose:** Nonlinear acceleration exploiting M=3 previous iterates. Works best when the iteration map is nearly linear and ρ is not yet stable. Provides 0.37–0.77 spectral ratio during the Anderson phase.

### Stage 4: Aitken δ² Extrapolation (Exponential regime)

**Activates when:** ρ has been constant (< 1% variance) for K ≥ 5 consecutive iterations and `dif < CHEB_THRESH`
**Method:** Quadratic extrapolation on dif monitor
```
dif_aitken = dif - (dif - dif_prev)² / (dif - 2*dif_prev + dif_prev2)
```
**Purpose:** Once linear exponential convergence is confirmed (ρ ≈ 0.84–0.87), Aitken δ² estimates the dif limit and modulates acceleration to jump toward it. Reduces exponential-phase iterations by ~20%.

**Safety bounds:**
- Aitken prediction must be 0 < dif_aitken < dif (toward zero, monotone)
- Fallback to standard Chebyshev if Aitken fails (non-monotone behavior)

### Convergence Monitoring

The solver tracks three key metrics:

1. **dif = |r_e_new / r_e_old - 1|** — equatorial-radius-based observable (outer loop convergence)
2. **Spectral radius estimate ρ_spec_est** — EMA-smoothed ratio of consecutive dif values
3. **Monotone-contraction counter n_consec_decrease** — consecutive dif reductions, gates late-phase accelerators

**Early termination criterion** (optional, looser than hard `dif < 1E-8`):
```fortran
if (dif < 1.E-5_wp .and. &
    maxval(abs(sphi - sphi_prev)) < eps .and. &
    maxval(abs(rho - prev_rho)) < eps) exit
```
Stops the Aitken tail early (~iter 35) when all fields have converged to machine precision, saving 40–50 iterations per solve.

### Typical Performance

- **Total iterations:** 40–90 (depending on warm-start state and EOS)
- **Picard phase:** 5–7 iterations
- **Chebyshev phase:** 2–3 iterations
- **Anderson phase:** 3–5 iterations
- **Aitken phase:** 30–75 iterations (scales logarithmically with target tolerance)

**Example run:**
```
iter=0-5:   Picard,    dif: 1.7E-2 → 4.5E-3
iter=6-7:   Chebyshev, dif: 4.5E-3 → 4.3E-3
iter=8-10:  Anderson,  dif: 4.3E-3 → 9.7E-4
iter=11-78: Aitken,    dif: 7.6E-4 → 9.6E-9 (67 iters at ρ ≈ 0.847)
```

### Tuning Guide

| Parameter | Default | Purpose | Adjust if… |
|-----------|---------|---------|-----------|
| `PICARD_THRESH` | 1.0E-1 | Picard→Chebyshev threshold | Divergence in dif ≥ 0.1 regime |
| `CHEB_THRESH` | 1.0E-3 | Chebyshev→Anderson threshold | Oscillation in 0.01–0.1 range |
| `w_picard` | 0.7 | Picard damping weight | Reduce to 0.5 for tighter coupling |
| `rho_spec_est` init | 0.7 | Initial spectral radius guess | Set from previous solve ρ_final |
| `N_MONOTONE` | 5 | Aitken activation gate | Reduce to 3 for earlier Aitken |

---

## Grid & Coordinates

The code uses a compactified meridional grid:
- **Radial**: `s ∈ [0, SMAX]` mapped to physical radius `r = r_e * s / (1 - s)` raised to power `s_pwr`
- **Angular**: `μ ∈ [0, 1]` (cosine of polar angle) on uniform, Legendre, or Chebyshev collocation points
- **Default resolution**: SDIV = 801, MDIV = 41 (configurable in `para_mod.f90`)
- **Interpolation**: barycentric Lagrange (4th-order stencil) for spectral-type accuracy

---

## Field Variables

**Metric** (all allocated in `para_mod.f90`):
- `gama(SDIV, MDIV)` — log conformal factor
- `rho(SDIV, MDIV)` — log metric deviation (oblate/prolate)
- `ww(SDIV, MDIV)` — angular velocity frame-dragging
- `alpha(SDIV, MDIV)` — lapse function
- `sphi(SDIV, MDIV)` — scalar field (ST theory only)

**Fluid & Thermodynamic**:
- `pressure`, `energy`, `enthalpy` — from EOS table lookup
- `velocity_sq` — fluid 3-velocity squared (must be < 1)
- `F_j` — angular-momentum flux

**Bulk Properties** (computed during each solve):
- `mass`, `ang_mom`, `chi` — ADM mass, angular momentum, spin parameter
- `r_e`, `r_ratio` — equatorial radius, polar-to-equatorial ratio
- `Omega_c`, `Omega_e` — angular velocity (center, equator)
- `I_inertia`, `Love2` — moment of inertia, tidal Love number

---

## Repository Layout

**src/core/** — central solver components
- `para_mod.f90` — global parameters, grid, field allocations, EOS/theory/rotation-law configuration
- `eos.f90` — EOS table loading and interpolation (log-space, supports phase transitions)
- `grid.f90` — compactified meridional grid construction
- `shoot_v2.f90` — Newton shooting driver; calls 1D/2D solvers in `shoot_solver_*_mod.f90`
- `analysis_v2.f90` — ADM mass, baryon mass, angular momentum, moment of inertia, Love number
- `constraint_eq.f90` — Hamiltonian constraint evaluation
- `starting_model.f90`, `sphere.f90` — static and initial-condition setup

**src/theory/** — rotation solver and relaxation
- `rotation_solver.f90` — outer loop iterating on equatorial radius
- `spin_integration.f90` — Green's function integrals (metric/scalar targets) via multipole expansion
- `spin_relaxation.f90`, `relaxation_mod.f90` — 4-stage adaptive relaxation (Picard → Chebyshev → Anderson → Aitken δ²)
- `spin_workspace.f90` — workspace allocation and geometry caching
- `spin_updates.f90`, `spin_derivatives.f90` — field updates and grid derivatives
- `rotational_law_mod.f90` — rotation law implementations (uniform, const-J, Uryu)

**src/tool/** — numerical utilities
- `toolkit_mod.f90` — barycentric Lagrange interpolation, spectral derivatives, integral accumulation
- `ad_mod.f90` — automatic differentiation (dual numbers) for dP/dE
- `cheb_mod.f90`, `spectral_hub.f90` — Chebyshev/Legendre spectral bases
- `brent.f90` — root finding for rotation laws
- `nag_compat_mod.f90` — adaptive quadrature integration

**eos/** — tabulated equation-of-state files (user-supplied)

**Cont/** — output directory: Omega.dat, velocity.dat, hamiltonian.dat

---

## Equation of State (EOS) Format

EOS tables are text files with format:
```
num_tab
log_e(1) log_p(1) log_h(1) log_n0(1)
...
log_e(num_tab) log_p(num_tab) log_h(num_tab) log_n0(num_tab)
```
where quantities are in **log-space** (e.g., `log_e = ln(e)`, energy density in CGS).
- `h = 1 + e/p + p/e` (relativistic enthalpy, dimensionless)
- `n0` is baryon number density
- Interpolation uses barycentric Lagrange in log-space with phase-transition detection

Default file: `eos/PS_G3.dat` (configurable in `para_mod.f90`).

---

## Solver Data Flow

The main computational loop (in `main.f90`):
1. **Initialize theory & grid** — set GR/ST, allocate grids in meridional coordinates
2. **Load EOS** — read tabulated table file
3. **Select task** — one of: `shoot_v2` (Newton shooting), `MRcurve` (M–R sequences), or `initialize_starting_model` (single model)
4. **Rotation solver** (in `rotation_solver.f90`):
   - Outer loop: iterate on equatorial radius `r_e`
   - Compute Green's function targets for metric and scalar field (via `spin_integration.f90`)
   - Inner loop: apply 4-stage relaxation to converge metric/scalar toward targets
   - Check convergence: `|r_e_new / r_e_old - 1| < 1e-7`
5. **Compute diagnostics** — ADM mass, baryon mass, angular momentum, Love number (in `analysis_v2.f90`)
6. **Evaluate constraint** — Hamiltonian constraint field written to `Cont/hamiltonian.dat`

---

## Build Instructions

### Prerequisites

- Fortran 2003+ compiler (tested with `gfortran`).
- BLAS and LAPACK libraries.
- The NAG D01 integration routine (`d01gaf`) or a drop-in replacement.
- EOS tables providing `log_e`, `log_p`, and `log_h` columns.

### Example Build

```bash
mkdir -p Cont
make
```

The Makefile compiles all sources in dependency order and places the
binary at `build/bin/a.out`. Override flags as needed:

```bash
make FFLAGS="-O3 -march=native" LIBS="-lopenblas -llapack"
```

Clean all build artefacts with:

```bash
make clean
```

---

## Usage

```
./build/bin/a.out
```

Theory (GR or ST) and other options are configured via parameters in
`src/core/para_mod.f90`. The program prints the active theory, loads the
EOS tables, assembles the grid, and iterates to a self-consistent
solution. Diagnostic files are
written to `./Cont`.

### Key Outputs

- `Cont/Omega.dat` — columns: circumferential radius, physical angular
  velocity, enthalpy, sound-speed squared, angular-momentum flux along
  the equatorial (`m=1`) grid line.
- `Cont/velocity.dat` — orbital velocity diagnostics (GR branch).
- Console log — convergence information, mass/radius, angular momentum.

---

## Extending GRASS

- **New diagnostics**: augment `src/core/analysis_v2.f90`; re-use
  `integrate_profiles` to accumulate multiple integrals efficiently.
- **Alternative rotation laws**: modify `para_mod.f90` defaults or add
  new parameters.
- **EOS handling**: adapt `loadEos` in `src/core/eos.f90` to support
  additional table formats.
- **Relaxation tuning**: The 4-stage relaxation strategy (Picard → Chebyshev
  → Anderson → Aitken) is implemented in `src/theory/spin_helper.f90`
  (`relaxation` subroutine, ~lines 988–1158) and `src/theory/relaxation_mod.f90`
  (Anderson and Aitken modules). Key thresholds and weights are defined as
  parameters near the start of `relaxation()`:
  - `PICARD_THRESH`, `CHEB_THRESH`: control stage transitions
  - `w_picard = 0.7`: Picard damping weight (reduce for tighter coupling)
  - `N_CHEB = 3`, `N_ANDERSON = 5`: gate activation of Chebyshev and Anderson
  - `N_RHO_LOCK = 3`, `N_AITKEN = 2`: Aitken history and ρ-lock requirements

  To modify convergence behavior: adjust thresholds, damping weights, or
  disable stages via logical guards (e.g., `if (.not. use_anderson)`).
  The Aitken δ² element-wise correction is guarded by `MAX_AITKEN_JUMP = 10`
  to prevent per-element amplification > 0.91× convergence rate.

Please keep changes Fortran 2003 compliant and document new options in
this README.

---

## Validation Tips

1. Compare global quantities (mass, radius, angular velocity) against
   published sequences or previous GRASS runs.
2. Inspect `Omega.dat` for smooth profiles and physical sound-speed
   bounds (`0 ≤ c_s^2 ≤ 1` in relativistic units).
3. Monitor Newton iteration logs for steady residual reduction.

---

## Support

If you encounter issues:

- Confirm your EOS tables match the expected format.
- Provide compiler version, OS, and build flags when reporting problems.
- Share console logs and generated `Cont/*.dat` excerpts for diagnosis.

Contributions (bug fixes, new physics modules, documentation) are
welcome—please open a pull request or issue describing the proposed
change.

Happy modelling!

