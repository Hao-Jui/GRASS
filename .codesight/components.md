# components.md — Module Public Interfaces

> GRASS has no UI components. This file documents every **Fortran module**
> as a "component": its role, public subroutine/function signatures, and
> dependencies.

Modules are grouped by layer: **Foundation**, **Numerics**, **Fields**,
**EOS**, **Grid**, **Theory**, **Solver**, **Diagnostics**, **Compat**.

---

## Layer: Foundation

### `precision_mod` — Working Precision
**File**: `src/macros/precision_mod.f90`
**Exports**:
- `integer, parameter :: wp = real64` — global working precision

---

### `constants_mod` — Physical Constants & Enums
**File**: `src/macros/constants_mod.f90`  **Depends on**: `precision_mod`
**Exports**:
- Theory enums: `THEORY_GR`, `THEORY_ST`
- Mode enums: `MODE_REGRID`, `MODE_DEFAULT`
- Task enums: `shoot`, `MRbuild`, `OneModel`
- Shooting enums: `SHOOT_FIX1_HC`, `SHOOT_FIX1_RP`, `SHOOT_2D`
- Physical constants: `C`, `G`, `MSUN`, `MB`, `PI`, `HBAR`, `N_SAT`,
  `L_UNI`, `KAPPA`, `KSCALE`, `E_SURFACE`, `P_SURFACE`, `ACCURACY`,
  `TOV_RMIN`, `RHO_UNI`, `PRS_UNI`, `F_UNI`, `SCALARTON`

---

### `grid_config_mod` — Grid Dimensions & Arrays
**File**: `src/macros/grid_config_mod.f90`  **Depends on**: `precision_mod`
**Exports**:
- Grid dimensions: `res`, `SDIV`, `MDIV`, `s_pwr`, `SMAX`, `DS`, `DM`,
  `S_E`, `LMAX`, `RDIV`
- Collocation enums: `COLLOCATION_UNI`, `COLLOCATION_LEG`, `COLLOCATION_CHEB`,
  `angular_collocation`
- Grid arrays: `s_gp(:)`, `mu(:)`, `sin_theta(:)`, `D_mu(:,:)`,
  `D_mu_t(:,:)`, `D2_mu(:,:)`, `w_mu(:)`
- `subroutine allocate_grid()` — allocate grid arrays
- `subroutine deallocate_grid()` — free grid arrays

---

## Layer: Fields

### `fields_mod` — Field Arrays & Bulk Properties
**File**: `src/para_panel.f90`  **Depends on**: `precision_mod`,
`constants_mod`, `grid_config_mod`
**Exports** (selected):
- `subroutine initialize_theory()` — set GR/ST mode, allocate all fields
- `subroutine apply_gr_defaults()` — zero scalar coupling
- `subroutine apply_st_defaults()` — enable scalar field
- `subroutine allocate_fields()` — allocate all 2D/1D/scalar arrays
- `subroutine deallocate_fields()` — free all arrays
- `pure function to_lower_str(str) result(out)` — string lowercase
- All field arrays, bulk properties, config (see `schema.md`)

---

## Layer: Numerics

### `toolkit_mod` — Numerical Utilities
**File**: `src/numerics/toolkit_mod.f90`
**Depends on**: `precision_mod`, `ad_mod`, `spectral_hub_mod`,
`grid_config_mod`, `fields_mod`
**Exports**:
- `function interp(xp, yp, x, np, xb_local) result(y)` — 4th-order
  barycentric Lagrange interpolation
- `function interp_dual(xp, yp_dual, x, np, xb_local) result(y_dual)`
  — dual-number version
- `subroutine deriv_s_1d(f, ds)` — 1D spectral radial derivative in-place
- `subroutine integrate_profiles(mu, integrand_buffer, integral_results)`
  — angular quadrature on `mu` grid; accumulates multiple integrands at once
- `function binary_search_index(arr, n, x) result(idx)` — bracket search
- `function pow_int_real(x, n) result(y)` — x^n for integer n
- `function expm1_safe(x) result(y)` — stable exp(x)-1
- `function clip_bessel(val) result(out)` — clamp Bessel outputs
- `function bessel_j0(x)`, `bessel_j1(x)`, `bessel_jn(n, x)` — Bessel
  functions via downward recursion
**Constants**: `n_order = 4` (Lagrange stencil half-width)

---

### `ad_mod` — Automatic Differentiation
**File**: `src/numerics/ad_mod.f90`
**Exports**: `type dual`, `dual_const(x)`, `dual_var(x)`, arithmetic
operators on `dual`. (See `schema.md` for type fields.)

---

### `cheb_mod` — Chebyshev Polynomials
**File**: `src/numerics/cheb_mod.f90`
**Exports**: Chebyshev polynomial evaluation and derivative routines.

---

### `spectral_hub_mod` — Collocation & Quadrature
**File**: `src/numerics/spectral_hub_mod.f90`
**Exports**:
- `subroutine gauss_lobatto(n, x, w)` — Legendre Gauss-Lobatto nodes & weights
- `subroutine chebyshev_lobatto_points(n, x)` — Chebyshev extrema
- `subroutine clenshaw_curtis_weights(n, w)` — Clenshaw-Curtis weights
- `subroutine barycentric_diff_matrices(x, D, D2)` — first and second
  derivative matrices via barycentric Lagrange
- `subroutine legendre_sequence(x, l_vals, p_vals)` — evaluate multiple
  Legendre polynomials at a single point

---

### `brent_mod` — Root Finding
**File**: `src/numerics/brent_mod.f90`
**Exports**:
- `subroutine brent_root(f, a, b, tol, root)` — Brent's method

---

### `nag_compat_mod` — Adaptive Quadrature Wrapper
**File**: `src/numerics/nag_compat_mod.f90`
**Exports**:
- `subroutine d01gaf_compat(fun, a, b, epsabs, epsrel, result, abserr)`
  — NAG D01GAF compatibility shim for adaptive 1D quadrature (no NAG
  licence required)

---

### `spline_mod` — Cubic Spline Interpolation
**File**: `src/numerics/spline_mod.f90`
**Exports**:
- `type spline_coeff` (see `schema.md`)
- `subroutine build_spline(xp, yp, n, coeff)` — not-a-knot cubic spline
- `subroutine build_spline_segmented(xp, yp, n, pt_idx, n_pt, coeff)`
  — spline with explicit phase-transition gaps

---

## Layer: EOS

### `eos_mod` — EOS Table Interpolation
**File**: `src/eos/eos_mod.f90`  **Depends on**: `fields_mod`, `spline_mod`,
`toolkit_mod`, `ad_mod`
**Exports**:
- `subroutine loadEos` — read `eos_file.dat`, detect phase transitions,
  build spline coefficients, wire up `eos_channel` instances
- `type eos_channel` — universal EOS lookup container; holds pointers to
  log-space table arrays (`xp`, `yp`) and spline coefficients (`sc`)
- `eos_e2p`, `eos_e2n0`, `eos_h2p`, `eos_h2e`, `eos_h2n0`, `eos_p2e`,
  `eos_p2h` — pre-wired channel instances for all 7 interpolation directions
- `function eos_eval(ch, x) result(y)` — universal scalar lookup:
  `y = exp(spline(log(x)))` using the channel's table/spline pointers
- `function eos_eval_dual(ch, x) result(y)` — dual-number version for AD
- `subroutine pe_at_logh_vec(lh, pp, ee, mask, n)` — vectorized
  simultaneous P and eps lookup from pre-logged enthalpy (hot-path)
- `subroutine pressure_derivative_n(ee, n, derivative [, status, idx_hint])`
  — nth derivative of P w.r.t. eps via Fornberg weights

---

## Layer: Grid

### `grid_mod` — Grid Construction
**File**: `src/grid/grid_mod.f90`  **Depends on**: `grid_config_mod`,
`toolkit_mod`, `fields_mod`
**Exports**:
- `subroutine make_grid` — populate `s_gp`, `mu`, derivative matrices
  `D_mu`, `D2_mu`, `D_mu_t`, weights `w_mu`
- `subroutine legendre_and_deriv(n, x, p, dp)` — evaluate Pn(x) and P'n(x)

---

## Layer: Diagnostics

### `ope_eq_mod` — Differential Operators
**File**: `src/diagnostics/ope_eq_mod.f90`  **Depends on**: `grid_config_mod`,
`fields_mod`
**Exports**:
- `type laplacian_operator` with methods `init()`, `lap2(f)`,
  `laplacian(f)`, `scal(v1,v2)`, `divr(v)` (see `schema.md`)

---

### `constraint_mod` — Hamiltonian Constraint
**File**: `src/diagnostics/constraint_mod.f90`  **Depends on**: `grid_config_mod`,
`fields_mod`, `ope_eq_mod`
**Exports**:
- `subroutine hamiltonian(hamL2)` — compute 2D Hamiltonian constraint
  field, return L2 norm `hamL2`

---

### `analysis_mod` — Diagnostics & Bulk Properties
**File**: `src/diagnostics/analysis_mod.f90`  **Depends on**: `fields_mod`,
`toolkit_mod`, `eos_mod`, `ad_mod`, `cheb_mod`, `nag_compat_mod`,
`exporter_mod`, `miscellaneous_mod`
**Exports**:
- `subroutine mass_radius()` — compute ADM mass, baryon mass, ang. mom.,
  moment of inertia, Love number; writes to `fields_mod` scalar outputs
- `subroutine solution_properties()` — wrapper for full diagnostics pass
- `subroutine prepare_common_data()` — precompute metric/scalar derivatives

---

### `donutization_mod` — Donutization Diagnostics
**File**: `src/diagnostics/donutization_mod.f90`  **Depends on**: `constants_mod`
**Exports**: routines to compute "donut number" order parameter for scalar
profiles.

---

### `exporter_mod` — Output Writers
**File**: `src/diagnostics/exporter_mod.f90`  **Depends on**: `eos_mod`,
`grid_config_mod`, `fields_mod`
**Exports**:
- `subroutine write_profiles(filename, ...)` — write 1D `.dat` profile
  files (columns: r, Omega, h, cs2, F_j, ...)

---

## Layer: Compat

### `miscellaneous_mod` — Utility Subroutines
**File**: `src/compat/miscellaneous_mod.f90`  **Depends on**: `nag_compat_mod`,
`fields_mod`
**Exports**:
- `subroutine print_converged_block(...)` — formatted convergence output

---

### `set_disk` — Disk Helper (Compatibility)
**File**: `src/compat/set_disk.f90`  **Depends on**: `toolkit_mod`, `fields_mod`,
`eos_mod`
Disk equilibrium helper kept for legacy compatibility.

---

### `sphere_mod` — Static Model Setup
**File**: `src/compat/sphere_mod.f90`  **Depends on**: `toolkit_mod`, `fields_mod`,
`eos_mod`, `set_disk`
**Exports**:
- `subroutine initialize_sphere()` — static spherically-symmetric config

---

### `regrid_mod` — Adaptive Regridding
**File**: `src/compat/regrid_mod.f90`  **Depends on**: `fields_mod`, `grid_mod`
**Exports**: regridding utilities for adaptive mesh refinement during
iteration.

---

## Layer: Theory

### `rotation_solver_mod` — Main Rotation Solver
**File**: `src/theory/rotation_solver_mod.f90`
**Depends on**: `constants_mod`, `grid_config_mod`, `fields_mod`,
`spin_workspace_mod`, `spin_integration_mod`,
`spin_relaxation_mod`, `spin_updates_mod`, `analysis_mod`
**Exports**:
- `subroutine rotation_solver` — outer iteration loop until convergence;
  orchestrates Green's-function target computation, 4-stage relaxation,
  field updates, and diagnostics

---

### `spin_integration_mod` — Green's Function Targets
**File**: `src/theory/spin_integration_mod.f90`
**Depends on**: `spin_workspace_mod`, `spin_derivatives_mod`, `toolkit_mod`,
`exporter_mod`, `fields_mod`, `eos_mod`, `donutization_mod`
**Exports**:
- `subroutine get_all_targets(r_e, root_mphi_re, target_rho, target_gama,`
  `  target_ww, target_sphi)` — compute multipole-expanded 2D target arrays
  for all metric potentials and scalar field
- `subroutine compute_metric_targets(...)` — internal; solve rho, gamma, omega
- `subroutine compute_scalar_targets(...)` — internal; solve phi
- `subroutine output_helper(...)` — write diagnostic `.dat` files

---

### `spin_relaxation_mod` — 4-Stage Relaxation Orchestrator
**File**: `src/theory/spin_relaxation_mod.f90`
**Depends on**: `spin_workspace_mod`, `fields_mod`, `relaxation_mod`
**Exports**:
- `subroutine relaxation(target_rho, target_gama, target_ww, target_sphi,`
  `  root_mphi_re, n_of_it, dif)` — apply one relaxation step using the
  active stage (Picard / Chebyshev / Anderson / Aitken) based on `dif`

---

### `relaxation_mod` — Acceleration Kernels
**File**: `src/theory/relaxation_mod.f90`  **Depends on**: `grid_config_mod`
**Exports**:
- `subroutine anderson_accel_optimized(x, f, hist_x, hist_f, m, n, iter)`
  — rank-M=3 Anderson mixing update
- `subroutine aitken_delta2(x_curr, x_prev, x_prev2, x_new, converged)`
  — Aitken delta2 quadratic extrapolation
- `subroutine aitken_reset()` — clear Aitken history

---

### `spin_updates_mod` — Field Update Subroutines
**File**: `src/theory/spin_updates_mod.f90`
**Depends on**: `toolkit_mod`, `fields_mod`, `rotational_law_mod`, `brent_mod`,
`eos_mod`
**Exports**:
- `subroutine update_equatorial_radius(r_e_old, r_e_new, dif, ...)`
- `subroutine update_angular_velocity(r_e, gama_pole, rho_pole, ...)`
- `subroutine update_eos_and_velocity(r_e, ...)` — interpolate P, h, eps, v
- `subroutine reset_uryu_peak_cache()` — clear Uryu rotation law cache

---

### `spin_derivatives_mod` — Grid Derivative Operators
**File**: `src/theory/spin_derivatives_mod.f90`
**Depends on**: `fields_mod`, `grid_config_mod`
**Exports**: spatial derivative operators ds, d_mu, nabla2 on compactified
meridional grid.

---

### `spin_workspace_mod` — Workspace Allocation
**File**: `src/theory/spin_workspace_mod.f90`
**Depends on**: `nag_compat_mod`, `fields_mod`, `grid_config_mod`
**Exports**:
- `subroutine allocate_workspace()` — allocate all target + derivative
  cache arrays listed in `schema.md`
- `subroutine deallocate_workspace()` — free all workspace arrays

---

### `rotational_law_mod` — Rotation Law Implementations
**File**: `src/theory/rotational_law_mod.f90`
**Depends on**: `brent_mod`, `fields_mod`, `grid_config_mod`
**Exports**:
- Rotation law subroutines for uniform, const-J, and Uryu laws
- `subroutine cache_uryu_ab()` — cache Uryu coefficients
- `subroutine shoot_Fmax(...)` — shoot for peak angular momentum flux

**Supported laws**:
| `solver_type` | Profile |
|---|---|
| `"uniform"` | Omega = const |
| `"const_j"` | Constant angular momentum flux; parameter `A_diff` |
| `"uryu"` | Power-law; parameters `lambda1`, `lambda2`, `uyru_p`, `uyru_q` |

---

## Layer: Solver

### `shoot_mod` — Newton Shooting Driver
**File**: `src/solver/shoot_mod.f90`  **Depends on**: `donutization_mod`,
`spin_workspace_mod`, `starting_model_mod`, `eos_mod`, `analysis_mod`,
`shoot_solver_2d_mod`, `shoot_solver_1d_hc_mod`, `shoot_solver_1d_r_ratio_mod`
**Exports**:
- `subroutine shoot_v2` — outer Newton loop to reach bulk-property targets;
  selects sub-mode (`SHOOT_FIX1_HC`, `SHOOT_FIX1_RP`, `SHOOT_2D`) based on
  user parameters

---

### `shoot_solver_2d_mod` — 2D Broyden Newton Solver
**File**: `src/solver/shoot_solver_2d_mod.f90`
**Exports**:
- `type newton_state` — Jacobian, history, step (see `schema.md`)
- `subroutine init_newton_state(state, nvar)`
- `subroutine reset_newton_state(state)`
- `subroutine solve_linear(state, rhs, delta_x)`
- `subroutine broyden_update(state, F_old, F_new, delta_x)`
- `subroutine clamp_step(delta_x, x_scale)`
- `subroutine from_solver_coords(er, rho0, x)`
- `subroutine to_solver_coords(x, er, rho0)`

---

### `shoot_solver_1d_hc_mod` — 1D Newton (h_center)
**File**: `src/solver/shoot_solver_1d_hc_mod.f90`
**Exports**: 1D Newton state and driver for fixing central enthalpy
`h_center`.

---

### `shoot_solver_1d_r_ratio_mod` — 1D Newton (axis ratio)
**File**: `src/solver/shoot_solver_1d_r_ratio_mod.f90`
**Exports**: 1D Newton state and driver for fixing axis ratio `r_ratio`.

---

### `starting_model_mod` — TOV Initial Seed
**File**: `src/solver/starting_model_mod.f90`  **Depends on**:
`rotation_solver_mod`, `miscellaneous_mod`, `scalar_burning_mod`, `eos_mod`,
`analysis_mod`, `regrid_mod`, `sphere_mod`
**Exports**:
- `subroutine initialize_starting_model()` — build spherically-symmetric
  (TOV) initial guess, populate all field arrays

---

### `scalar_burning_mod` — Scalar Field Initialization
**File**: `src/solver/scalar_burning_mod.f90`
**Exports**:
- `subroutine scalar_burning()` — iterative initialization of scalar field
  phi for ST gravity runs

---

### `MRcurve_mod` — Mass-Radius Sequence
**File**: `src/solver/MRcurve_mod.f90`  **Depends on**: `donutization_mod`,
`spin_workspace_mod`, `rotation_solver_mod`, `starting_model_mod`, `eos_mod`,
`analysis_mod`
**Exports**:
- `subroutine MRcurve()` — sweep `Mb_goal`, call `rotation_solver` per
  step, write `properties.dat`
