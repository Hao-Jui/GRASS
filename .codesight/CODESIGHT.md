# CODESIGHT.md — GRASS Context Map

> Read this file first. It gives you a complete mental model of GRASS before
> you touch any source file. Follow the links to specialist files for details.

---

## What Is GRASS?

**GRASS** (General Relativistic Axisymmetric Spacetime Solver) computes
self-consistent rotating neutron-star equilibria in General Relativity or
Scalar-Tensor gravity. Given an equation of state and target bulk properties
(mass, spin, ...), it iterates five 2D field arrays on a compactified
meridional grid until convergence.

**Language**: Fortran 2003+ (~50 modules, ~8 000 lines)
**Build**: `make` (gfortran, BLAS/LAPACK required)
**Run**: `./build/bin/a.out`
**Output**: solution files in `Cont/`, profile in `Cont/Omega.dat`

---

## Quick Start

```bash
# 1. Edit targets in src/para_panel.f90
#    - eos_file, Mb_goal, solver_type, active_theory

# 2. Build and run
make -j8 && ./build/bin/a.out

# 3. Parametric sweep over magnetic coupling
./sweep_b_goal.sh 3800 4800 100

# 4. Post-process in MATLAB
matlab -r "run('auto2D.m')"
```

---

## Repository Layout

```
GRASS/
├── src/
│   ├── main.f90              Program entry point
│   ├── macros/               Root types & constants
│   │   ├── precision_mod     wp = real64
│   │   ├── constants_mod     Physical constants & enums
│   │   └── grid_config_mod   Grid dimensions, arrays, alloc
│   ├── numerics/             Numerical utilities
│   │   ├── ad_mod            Forward-mode AD (dual numbers)
│   │   ├── brent_mod         Root finding
│   │   ├── cheb_mod          Chebyshev polynomials
│   │   ├── nag_compat_mod    Adaptive quadrature
│   │   ├── spectral_hub_mod  Collocation & diff matrices
│   │   ├── spline_mod        Cubic spline interpolation
│   │   └── toolkit_mod       Interp, integrate, Bessel, deriv
│   ├── fields/               Field arrays & bulk properties
│   │   └── fields_mod        (formerly para_mod)
│   ├── eos/                  Equation of state
│   │   └── eos_mod           EOS table interpolation
│   ├── grid/                 Grid construction
│   │   └── grid_mod          Compactified grid + D matrices
│   ├── theory/               Rotation solver & relaxation
│   │   ├── rotation_solver_mod
│   │   ├── spin_integration_mod
│   │   ├── spin_relaxation_mod
│   │   ├── spin_updates_mod
│   │   ├── spin_derivatives_mod
│   │   ├── spin_workspace_mod
│   │   ├── rotational_law_mod
│   │   └── relaxation_mod
│   ├── solver/               Newton shooting & drivers
│   │   ├── shoot_mod
│   │   ├── shoot_solver_{2d,1d_hc,1d_r_ratio}_mod
│   │   ├── starting_model_mod
│   │   ├── scalar_burning_mod
│   │   └── MRcurve_mod
│   ├── diagnostics/          Analysis & output
│   │   ├── analysis_mod
│   │   ├── constraint_mod
│   │   ├── ope_eq_mod
│   │   ├── donutization_mod
│   │   └── exporter_mod
│   └── compat/               Legacy compatibility
│       ├── miscellaneous_mod
│       ├── set_disk
│       ├── sphere_mod
│       └── regrid_mod
├── eos/                      EOS table files (*.dat)
├── Cont/                     Solution output files
├── tests/                    Unit tests (spline)
├── Res/, Map/                 Result & profile storage
├── Makefile                  Build system
├── sub.sh                    Clean-build-run helper
├── sweep_b_goal.sh           Parametric sweep script
├── auto2D.m, donut.m         MATLAB post-processing
└── .codesight/               <- you are here
```

---

## Context Map Files

| File | Contents |
|---|---|
| [`routes.md`](routes.md) | Every solver entry point: `shoot_v2`, `MRcurve`, `initialize_starting_model`; inputs, outputs, call chains; shell/MATLAB entry points |
| [`schema.md`](schema.md) | Every derived type (`dual`, `spline_coeff`, `laplacian_operator`, `newton_state`) and every allocatable array (metric, fluid, grid, EOS, workspace) with shape and physical meaning |
| [`components.md`](components.md) | Every module with its public subroutine/function signatures and layer (Foundation / Numerics / Fields / EOS / Grid / Theory / Solver / Diagnostics / Compat) |
| [`libs.md`](libs.md) | External libraries (BLAS, LAPACK, FFTW3), every numerics-module export with signatures, EOS file formats |
| [`config.md`](config.md) | Every configurable parameter: target quantities, EOS config, grid config, theory selection, rotation law, convergence thresholds, physical constants, build flags, output paths |
| [`middleware.md`](middleware.md) | 8-stage solver pipeline: initialization -> EOS load -> grid -> TOV seed -> Newton loop -> 4-stage relaxation -> field update -> diagnostics; all error guards and fallbacks |
| [`graph.md`](graph.md) | Full USE-dependency graph; fan-in ranking; top-10 highest-risk files; compilation order; safe-to-change leaf modules |

---

## Module Hierarchy (Quick Reference)

```
main.f90
  |
  +--- precision_mod     <- root: wp = real64
  +--- constants_mod     <- physical constants, enum parameters
  +--- grid_config_mod   <- grid dimensions, arrays, alloc/dealloc
  +--- fields_mod        <- field arrays, bulk properties, config
  |
  +--- eos_mod           <- EOS table + eos_channel containers
  +--- grid_mod          <- compactified grid + derivative matrices
  |
  +--- rotation_solver_mod    <- outer iteration loop
  |     +- spin_integration_mod   <- Green's function targets
  |     +- spin_relaxation_mod    <- 4-stage: Picard->Cheb->Anderson->Aitken
  |     |    +- relaxation_mod    <- Anderson & Aitken kernels
  |     +- spin_updates_mod       <- apply targets, EOS & velocity update
  |          +- rotational_law_mod <- uniform / const-J / Uryu
  |
  +--- shoot_mod       <- Newton shooting driver (shoot_v2)
  |     +- shoot_solver_{2d,1d_hc,1d_r_ratio}_mod
  |
  +--- MRcurve_mod     <- mass-radius sequence builder
  +--- analysis_mod    <- ADM mass, baryon mass, Love number, multipoles
  +--- constraint_mod  <- Hamiltonian constraint L2 check
  |
  +--- toolkit_mod     <- interp, integrate_profiles, Bessel, deriv
        +- ad_mod           <- forward-mode AD (dual numbers)
        +- spectral_hub_mod <- Gauss-Lobatto, Chebyshev, diff matrices
        +- cheb_mod         <- Chebyshev polynomials
        +- spline_mod       <- cubic spline (with phase-transition gaps)
        +- brent_mod        <- root finding
        +- nag_compat_mod   <- adaptive quadrature (no NAG licence)
        +- exporter_mod     <- write .dat output files
```

---

## Call Chain: Single Solution (`shoot_v2`)

```
main -> shoot_v2
  +- initialize_starting_model     [starting_model_mod]
  +- Newton outer loop
      +- rotation_solver           [rotation_solver_mod]
          +- update_eos_and_velocity
          +- get_all_targets        [Green's function multipole integrals]
          +- relaxation             [Picard -> Chebyshev -> Anderson -> Aitken]
          +- update_alpha_lapse
          +- mass_radius            [ADM mass, baryon mass, chi, Love2...]
  +- write_profiles                [output to Cont/]
  +- hamiltonian                   [constraint L2 check]
```

---

## Field Variables (2D arrays, shape SDIV x MDIV = 801 x 801)

| Variable | Physical quantity |
|---|---|
| `gama` | Conformal factor gamma |
| `rho` | Oblate distortion rho |
| `ww` | Frame-dragging omega |
| `alpha` | Lapse alpha |
| `sphi` | Scalar field phi (ST only; zero in GR) |
| `pressure` | Fluid pressure P |
| `enthalpy` | Specific enthalpy h |
| `energy` | Energy density eps |
| `omg` | Angular velocity Omega(r,mu) |

---

## Key Outputs

| Quantity | Written to | Module |
|---|---|---|
| M, Mb, chi, r_e, Love2, I | fields_mod scalars + console | `analysis_mod` |
| Omega(r) radial profile | `Cont/Omega.dat` | `exporter_mod` |
| Full 2D solution | `Cont/{EOS}_J{spin}_Mb{mass}_...dat` | `exporter_mod` |
| M-R curve | `properties.dat` | `MRcurve_mod` |
| Hamiltonian constraint | console (hamL2) | `constraint_mod` |

---

## Highest-Risk Files (from graph.md)

| Risk | File | Reason |
|---|---|---|
| :red_circle: | `src/para_panel.f90` | Many dependents; all fields live here |
| :red_circle: | `src/macros/precision_mod.f90` | `wp` type used in every module |
| :orange_circle: | `src/macros/constants_mod.f90` | Physical constants used throughout |
| :orange_circle: | `src/macros/grid_config_mod.f90` | Grid arrays used by most modules |
| :orange_circle: | `src/numerics/toolkit_mod.f90` | `interp` + `integrate_profiles` used throughout |
| :orange_circle: | `src/eos/eos_mod.f90` | EOS functions inside tight iteration loops |
| :orange_circle: | `src/theory/spin_relaxation_mod.f90` | Convergence behavior of all runs |

---

## Rotation Laws

| `solver_type` | Profile | Key params |
|---|---|---|
| `"uniform"` | Omega = const | -- |
| `"const_j"` | constant angular momentum flux | `A_diff` |
| `"uryu"` | power-law differential | `lambda1`, `lambda2`, `uyru_p`, `uyru_q` |

---

## Relaxation Stages (spin_relaxation_mod)

| Stage | Condition | Method | Key params |
|---|---|---|---|
| Picard | dif >= 0.5 | weighted average | W_PICARD = 0.7 |
| Chebyshev | 0.1 <= dif < 0.5 | spectral-radius adaptive weight | rho_spec estimated on-the-fly |
| Anderson | 0.0001 <= dif < 0.1 | rank-3 mixing | M_HIST = 3 |
| Aitken delta2 | dif < 1e-5 | quadratic extrapolation | AITKEN_COOLDOWN = 3 |

---

## How to Add Things

**New EOS**: place `{name}.dat` (format: `num_tab`, then rows of
`log_e log_p log_h log_n0`) in `eos/`; set `eos_file = "{name}"` in
`fields_mod`.

**New rotation law**: add a case to `rotational_law_mod.f90`; reuse
`brent_root` for root-finding; add a new `solver_type` string.

**New diagnostic**: add a subroutine to `analysis_mod.f90`; call it from
the post-convergence block in `rotation_solver_mod.f90`.

**New module**: add to `SOURCES` in `Makefile`; add explicit dependency
lines (see Makefile lines 81-136).
