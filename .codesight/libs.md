# libs.md — Libraries & Utility Exports

> This file covers every external library used by GRASS and every
> utility/tool module that is consumed by multiple other modules.

---

## External Libraries

### BLAS / LAPACK
**Linked via**: `LIBS ?= -llapack -lblas` in Makefile
**Used by**: matrix operations inside `shoot_solver_2d_mod` (Jacobian
solves), spectral derivative matrices in `spectral_hub_mod`
**Required**: yes — build fails without BLAS/LAPACK

---

### FFTW3
**Linked via**: `FFTW_LIBS ?= $(shell pkg-config --libs fftw3 2>/dev/null || echo -lfftw3)`
**Used by**: spectral transforms (if enabled); gracefully falls back to
`-lfftw3` if `pkg-config` is absent
**Required**: optional — code compiles without it if FFTW paths are not
referenced

---

### gfortran Runtime
**Compiler**: `gfortran` (FC = gfortran)
**Standard**: Fortran 2003+ (`real64` from `iso_fortran_env`)
**Flags (release)**: `-O3 -march=native`
**Flags (debug)**: `-O0 -g -fbacktrace -Wall -Wextra -Wimplicit-interface`
  `-fcheck=all -Wuninitialized -Wconversion -Wuse-without-only -finit-real=nan`
**Stack size**: `-Wl,-stack_size,0x4000000` (macOS; 64 MB)

---

## Internal Utility Libraries

### `toolkit_mod` — Core Numerical Utilities
**File**: `src/numerics/toolkit_mod.f90`
**Consumed by**: `eos_mod`, `analysis_mod`, `spin_integration_mod`,
`spin_relaxation_mod`, `constrain_mod`, almost every theory module.

| Export | Signature | Purpose |
|---|---|---|
| `interp` | `(xp, yp, x, np, xb_local) → real(wp)` | 4th-order barycentric Lagrange |
| `interp_dual` | `(xp, yp_dual, x, np, xb_local) → dual` | Same with AD dual numbers |
| `deriv_s_1d` | `(f(:), ds)` in-place | Spectral radial derivative |
| `integrate_profiles` | `(mu(:), integrand(:,:), result(:))` | Multi-column angular quadrature |
| `binary_search_index` | `(arr, n, x) → integer` | Interpolation bracket |
| `pow_int_real` | `(x, n) → real(wp)` | xⁿ, integer exponent |
| `expm1_safe` | `(x) → real(wp)` | Stable exp(x)−1 |
| `clip_bessel` | `(val) → real(wp)` | Clamp Bessel values |
| `bessel_j0` | `(x) → real(wp)` | J₀ Bessel function |
| `bessel_j1` | `(x) → real(wp)` | J₁ Bessel function |
| `bessel_jn` | `(n, x) → real(wp)` | Jₙ Bessel function (downward recursion) |

**Constants**:
- `n_order = 4` — Lagrange interpolation stencil half-width
- `BARY_W(-4:4)` — precomputed barycentric weights for uniform grids

---

### `ad_mod` — Forward-Mode Automatic Differentiation
**File**: `src/numerics/ad_mod.f90`
**Consumed by**: `eos_mod` (dP/dE), `toolkit_mod`

| Export | Purpose |
|---|---|
| `type dual` | Dual number: `val` + `der` fields |
| `dual_const(x)` | Dual with zero derivative |
| `dual_var(x)` | Dual with unit derivative |
| `+, -, *, /` on `dual` | Propagate derivative through arithmetic |
| `exp(dual)` | exp with chain rule |
| `log(dual)` | log with chain rule |
| `abs(dual)` | abs with sign rule |

---

### `spectral_hub_mod` — Collocation Points & Quadrature
**File**: `src/numerics/spectral_hub_mod.f90`
**Consumed by**: `grid_mod`, `toolkit_mod`, `cheb_mod`

| Export | Signature | Purpose |
|---|---|---|
| `gauss_lobatto` | `(n, x(:), w(:))` | Legendre Gauss-Lobatto nodes + weights |
| `chebyshev_lobatto_points` | `(n, x(:))` | Chebyshev extrema |
| `clenshaw_curtis_weights` | `(n, w(:))` | Clenshaw-Curtis quadrature weights |
| `barycentric_diff_matrices` | `(x(:), D(:,:), D2(:,:))` | 1st + 2nd derivative matrices |
| `legendre_sequence` | `(x, l_max, p(:))` | Pₗ(x) for l=0..l_max |

---

### `cheb_mod` — Chebyshev Spectral Basis
**File**: `src/numerics/cheb_mod.f90`
**Consumed by**: `spectral_hub_mod`, Chebyshev acceleration in
`spin_relaxation_mod`

Exports: Chebyshev polynomial evaluation Tₙ(x), first and second
derivatives T'ₙ(x), T''ₙ(x).

---

### `spline_mod` — Cubic Spline Interpolation
**File**: `src/numerics/spline_mod.f90`
**Consumed by**: `eos_mod` (EOS table splines)

| Export | Signature | Purpose |
|---|---|---|
| `build_spline` | `(xp, yp, n, coeff)` | Not-a-knot cubic spline |
| `build_spline_segmented` | `(xp, yp, n, pt_idx, n_pt, coeff)` | Spline with phase-transition gaps |

Eval is done by private kernels in `eos_mod` (`sp_eval`, `sp_eval_pair`, `sp_eval_dual`) using a hinted interval search `sp_locate`.

**Type** `spline_coeff`: fields `b(:), c(:), d(:)` (see `schema.md`).

---

### `brent_mod` — Root Finding
**File**: `src/numerics/brent_mod.f90`
**Consumed by**: `rotational_law_mod`

| Export | Signature | Purpose |
|---|---|---|
| `brent_root` | `(f, a, b, tol, root)` | Brent's bracketed root-finding |

`f` is a `real(wp)` function with a single `real(wp)` argument, passed as
a procedure argument.

---

### `nag_compat_mod` — Adaptive Quadrature Shim
**File**: `src/numerics/nag_compat_mod.f90`
**Consumed by**: `analysis_mod`, `constrain_mod`

| Export | Signature | Purpose |
|---|---|---|
| `d01gaf_compat` | `(fun, a, b, epsabs, epsrel, result, abserr)` | NAG D01GAF drop-in; adaptive 1D Gauss quadrature (no NAG licence) |

---

### `relaxation_mod` — Acceleration Kernels
**File**: `src/theory/relaxation_mod.f90`
**Consumed by**: `spin_relaxation_mod`

| Export | Signature | Purpose |
|---|---|---|
| `anderson_accel_optimized` | `(x, f, hist_x, hist_f, m, n, iter)` | Rank-M=3 Anderson mixing |
| `aitken_delta2` | `(x_curr, x_prev, x_prev2, x_new, converged)` | Aitken δ² extrapolation |
| `aitken_reset` | `()` | Clear Aitken history buffers |

---

### `donutization_mod` — Order Parameter Diagnostics
**File**: `src/diagnostics/donutization_mod.f90`
**Consumed by**: `analysis_mod`, MATLAB post-processing

Exports routines to compute "donut number" (tail-power integral over high-
density window) for scalar-field profiles — the primary observable for
donutization in ST gravity.

---

### `exporter_mod` — Output I/O
**File**: `src/diagnostics/exporter_mod.f90`
**Consumed by**: `shoot_mod`, `MRcurve_mod`, `rotation_solver_mod`

| Export | Signature | Purpose |
|---|---|---|
| `write_profiles` | `(filename, ...)` | Write 1D `.dat` profile to `Cont/` |

---

## EOS Data Files

Located in `eos/`. Format: ASCII, first line `num_tab`, then `num_tab` rows
each with:
```
log10(ε)  log10(P)  log10(h)  log10(nB)
```
All in CGS (ε in g/cm³, P in dyn/cm², nB in cm⁻³).

| File | EOS Family | Notes |
|---|---|---|
| `ALF2.dat` | Relativistic Mean Field | Stiff |
| `APR1.dat` | Akmal-Pandharipande-Ravenhall | Soft |
| `MPA1.dat` | MIT Bag (RMF) | With varying B, mphi |
| `PS.dat` | Piecewise Polytropic | Reference |
| `PS_G3.dat` | PP with GST3 | Stiffer variant |
| `PSt.dat` | PP stiff | |
| `Hadronic_{1,3,5}_S1.dat` | Hadronic phases | |
| `hybrid_{1,3}_{S1,T0}.dat` | Hybrid quark-matter | Phase transition |
| `Delta_7_S1.dat` | Delta baryon | |
| `Zdunik.dat` | Zdunik parametrization | |
| `muses.dat` | MUSES collaboration | |
| `test.dat` | Small test EOS | Quick validation |

Phase transitions detected automatically by `loadEos` when pressure is
non-monotone; gaps stored in `p_at_PT(:)`.
