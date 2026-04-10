# schema.md — Data Structures

> GRASS has no database. This file documents every **Fortran derived type**
> and every significant **allocatable array** in `fields_mod.f90`,
> `grid_config_mod.f90`, and the workspace modules — the closest equivalents
> to a database schema.

---

## Derived Types

### `dual` — Automatic Differentiation Dual Number
**Module**: `ad_mod` (`src/numerics/ad_mod.f90`)

| Field | Type | Meaning |
|---|---|---|
| `val` | `real(8)` | Function value |
| `der` | `real(8)` | Derivative (forward-mode) |

**Operators overloaded**: `+  -  *  /  exp  log  abs`
**Constructors**: `dual_const(x)` (zero derivative), `dual_var(x)` (unit derivative)
**Use case**: forward-mode AD for dP/dE in EOS lookups.

---

### `spline_coeff` — Cubic Spline Coefficients
**Module**: `spline_mod` (`src/numerics/spline_mod.f90`)

| Field | Type | Meaning |
|---|---|---|
| `b(:)` | `real(wp), allocatable` | Linear coefficient per segment |
| `c(:)` | `real(wp), allocatable` | Quadratic coefficient |
| `d(:)` | `real(wp), allocatable` | Cubic coefficient |

Built by `build_spline(xp, yp, n, coeff)` or
`build_spline_segmented(xp, yp, n, pt_idx, n_pt, coeff)` for phase-
transition gaps.

---

### `eos_channel` — Universal EOS Lookup Container
**Module**: `eos_mod` (`src/eos/eos_mod.f90`)

| Field | Type | Meaning |
|---|---|---|
| `xp` | `real(wp), pointer` | Log-space input table (e.g. `log_e`, `log_h`, `log_p`) |
| `yp` | `real(wp), pointer` | Log-space output table (e.g. `log_p`, `log_e`, `log_n0`) |
| `sc` | `type(spline_coeff), pointer` | Cubic spline coefficients for this direction |

**Pre-wired instances** (set by `loadEos`):
`eos_e2p`, `eos_e2n0`, `eos_h2p`, `eos_h2e`, `eos_h2n0`, `eos_p2e`, `eos_p2h`

**Usage**: `y = eos_eval(eos_h2e, h_val)` replaces the former `e_at_h(h_val)`.
For AD: `y_dual = eos_eval_dual(eos_e2p, e_dual)`.

---

### `laplacian_operator` — Differential Operator Cache
**Module**: `ope_eq_mod` (`src/diagnostics/ope_eq_mod.f90`)

| Method / Field | Meaning |
|---|---|
| `init()` | Allocate derivative-matrix cache |
| `lap2(f)` | ∇² on compactified 2D grid |
| `laplacian(f)` | Full Laplacian (alias) |
| `scal(v1, v2)` | Dot product of two vector fields |
| `divr(v)` | Divergence |

---

### `newton_state` — 2D Broyden Newton State
**Module**: `shoot_solver_2d_mod` (`src/solver/shoot_solver_2d_mod.f90`)

| Field | Type | Meaning |
|---|---|---|
| `J(:,:)` | `real(wp), allocatable` | Approximate Jacobian (2×2) |
| `J_inv(:,:)` | `real(wp), allocatable` | Inverse Jacobian |
| `F_old(:)` | `real(wp), allocatable` | Previous residual vector |
| `delta_x(:)` | `real(wp), allocatable` | Previous step |
| `n_iter` | `integer` | Iteration count |

**Subroutines**:
- `init_newton_state(state, nvar)` — allocate
- `reset_newton_state(state)` — clear history
- `solve_linear(state, rhs, delta_x)` — solve Jδ = rhs
- `broyden_update(state, F_old, F_new, delta_x)` — rank-1 Broyden update
- `clamp_step(delta_x, x_scale)` — limit step size
- `from_solver_coords(er, rho0, x)` / `to_solver_coords(x, er, rho0)` — coordinate transforms

---

## Global Allocatable Arrays (fields_mod.f90 + grid_config_mod.f90)

All 2D arrays have shape `(SDIV, MDIV)` unless noted. `SDIV = MDIV = 801`
by default.

### Metric Potentials

| Name | Shape | Type | Physical meaning |
|---|---|---|---|
| `gama` | (SDIV, MDIV) | `real(wp)` | Conformal factor γ |
| `rho` | (SDIV, MDIV) | `real(wp)` | Oblate distortion ρ |
| `ww` | (SDIV, MDIV) | `real(wp)` | Frame-dragging potential ω |
| `alpha` | (SDIV, MDIV) | `real(wp)` | Lapse function α |
| `sphi` | (SDIV, MDIV) | `real(wp)` | Scalar field φ (ST gravity; zero in GR) |

### Fluid Fields

| Name | Shape | Type | Physical meaning |
|---|---|---|---|
| `pressure` | (SDIV, MDIV) | `real(wp)` | Fluid pressure P |
| `enthalpy` | (SDIV, MDIV) | `real(wp)` | Specific enthalpy h |
| `energy` | (SDIV, MDIV) | `real(wp)` | Energy density ε |
| `velocity_sq` | (SDIV, MDIV) | `real(wp)` | v² (velocity squared) |
| `omg` | (SDIV, MDIV) | `real(wp)` | Angular velocity Ω(r,μ) |
| `F_j` | (SDIV, MDIV) | `real(wp)` | Angular momentum flux |

### 1D Velocity/Sound Profiles

| Name | Shape | Physical meaning |
|---|---|---|
| `v_plus` | (SDIV) | Co-rotating velocity at equator |
| `v_minus` | (SDIV) | Counter-rotating velocity |
| `V_rr_p` | (SDIV) | Radial velocity + |
| `V_rr_m` | (SDIV) | Radial velocity − |
| `sound_speed` | (SDIV) | Local sound speed cₛ |

### Grid Arrays (grid_config_mod.f90)

| Name | Shape | Physical meaning |
|---|---|---|
| `s_gp` | (SDIV) | Compactified radial grid points s in [0, SMAX] |
| `mu` | (MDIV) | Angular grid points mu = cos theta in [-1, 1] |
| `sin_theta` | (MDIV) | sin theta at each angular point |
| `D_mu` | (MDIV, MDIV) | First angular derivative matrix d/dmu |
| `D_mu_t` | (MDIV, MDIV) | Transposed derivative matrix |
| `D2_mu` | (MDIV, MDIV) | Second angular derivative matrix d2/dmu2 |
| `w_mu` | (MDIV) | Angular quadrature weights |

### Green's Function Tables

| Name | Shape | Physical meaning |
|---|---|---|
| `P_2n` | (MDIV, N_L) | Even-order Legendre P₂ₙ(μ) evaluated at grid |
| `P1_2n_1` | (MDIV, N_L) | Legendre first derivative P'_{2n+1}(μ) |
| `sin_2n_1_theta` | (MDIV, N_L) | sin^{2n+1}θ at angular grid |

### EOS Tables (allocated by `loadEos`)

| Name | Shape | Physical meaning |
|---|---|---|
| `log_p` | (num_tab) | ln(P · KSCALE) — natural log, CGS units scaled |
| `log_e` | (num_tab) | ln(ε · C² · KSCALE) — natural log, CGS units scaled |
| `log_h` | (num_tab) | ln(ln(h_raw)) — double-log encoding of specific enthalpy |
| `log_n0` | (num_tab) | ln(nB) — natural log baryon density |
| `p_at_PT` | (n_PT) | **Integer** indices (not pressure values) of phase-transition rows in the table |

### Workspace Arrays (spin_workspace_mod.f90)

These are allocated during `rotation_solver` and freed after:

| Name | Shape | Physical meaning |
|---|---|---|
| `target_rho` | (SDIV, MDIV) | Green's-function target for ρ |
| `target_gama` | (SDIV, MDIV) | Green's-function target for γ |
| `target_ww` | (SDIV, MDIV) | Green's-function target for ω |
| `target_sphi` | (SDIV, MDIV) | Green's-function target for φ |
| `dg_s_cache` | (SDIV, MDIV) | ∂ₛγ derivative cache |
| `dg_m_cache` | (SDIV, MDIV) | ∂_μγ derivative cache |
| `dr_s_cache` | (SDIV, MDIV) | ∂ₛρ derivative cache |
| `dr_m_cache` | (SDIV, MDIV) | ∂_μρ derivative cache |
| `dww_s_cache` | (SDIV, MDIV) | ∂ₛω derivative cache |
| `dww_m_cache` | (SDIV, MDIV) | ∂_μω derivative cache |
| `ds_s_cache` | (SDIV, MDIV) | ∂ₛφ derivative cache |
| `ds_m_cache` | (SDIV, MDIV) | ∂_μφ derivative cache |
| `d2g_ss_cache` | (SDIV, MDIV) | ∂²ₛₛγ second-derivative cache |
| `d2g_mm_cache` | (SDIV, MDIV) | ∂²_μμγ second-derivative cache |

### Scalar Bulk Outputs (fields_mod.f90)

| Name | Type | Physical meaning |
|---|---|---|
| `mass` | `real(wp)` | ADM gravitational mass M [M☉] |
| `mass_0` | `real(wp)` | Baryon mass Mb [M☉] |
| `mass_s` | `real(wp)` | Scalar-active mass component |
| `mass_p` | `real(wp)` | Pressure-weighted mass |
| `ang_mom` | `real(wp)` | Angular momentum J |
| `chi` | `real(wp)` | Dimensionless spin χ = J/M² |
| `T_kin` | `real(wp)` | Kinetic energy |
| `r_e` | `real(wp)` | Equatorial circumferential radius [km] |
| `r_ratio` | `real(wp)` | Polar-to-equatorial ratio rp/re |
| `r_circ` | `real(wp)` | Circumferential equatorial radius |
| `Omega_c` | `real(wp)` | Central angular velocity Ωc [rad/s] |
| `Omega_e` | `real(wp)` | Equatorial angular velocity |
| `Omega_K` | `real(wp)` | Keplerian (mass-shedding) frequency |
| `I_inertia` | `real(wp)` | Moment of inertia I |
| `Love2` | `real(wp)` | Tidal Love number k₂ |
| `M2, M4, M6` | `real(wp)` | Mass multipoles |
| `S3, S5` | `real(wp)` | Current multipoles |

---

## Output File Format

Each solution file in `Cont/` (`{EOS}_J{spin}_Mb{mass}_...dat`) contains
columns: `s`, `mu`, `gama`, `rho`, `ww`, `alpha`, `sphi`, `pressure`,
`enthalpy`, `energy`, `velocity_sq`, `omg`.

`Omega.dat` columns: `r [km]`, `Ω [rad/s]`, `h`, `cs²`, `F_j`
