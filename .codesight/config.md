# config.md — Configuration Reference

> All configuration in GRASS is done via **Fortran namelist** parameters
> defined in `src/para_panel.f90` (field/solver config) and
> `src/macros/` (constants, grid dimensions) and read at program startup.
> There are no environment variables or external config files beyond
> the EOS table path.

---

## How to Configure

Edit the parameter defaults in `src/para_panel.f90` and recompile,
**or** provide a Fortran namelist input file at runtime (if the main program
reads one — see `src/main.f90`).

For parametric sweeps, use `sweep_b_goal.sh` which patches `B_goal` in
`para_panel.f90` via `perl` substitution and recompiles.

---

## Task Selection

| Parameter | Type | Default | Options | Meaning |
|---|---|---|---|---|
| *(implicit)* | — | `shoot_v2` | `shoot_v2` / `MRcurve` / `initialize_starting_model` | Controlled by compile-time or runtime flags in `main.f90` |

---

## Target Quantities (Shooting Goals)

Set the target for the Newton shooting solver. Only the active target(s)
need be set; others are ignored depending on `SHOOT_*` mode.

| Parameter | Type | Default | Unit | Meaning |
|---|---|---|---|---|
| `M_goal` | `real(wp)` | `1.2` | M☉ | Gravitational mass target |
| `Mb_goal` | `real(wp)` | `1.8` | M☉ | Baryon mass target |
| `J_goal` | `real(wp)` | `1.6` | (dimensionless) | Angular momentum target |
| `chi_goal` | `real(wp)` | `0.1` | — | Dimensionless spin χ = J/M² target |
| `omc_goal` | `real(wp)` | `30.0` | rad/s | Central angular velocity target |
| `B_goal` | `real(wp)` | `35.0` | (ST coupling) | Scalar-tensor coupling strength target |
| `mphi_goal` | `real(wp)` | `0.2` | eV | Scalar field mass target |

---

## EOS Configuration

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `eos_file` | `character(128)` | `"test"` | EOS filename without `.dat` extension; looked up in `eos/` directory |
| `phase_transition` | `logical` | `.false.` | Set `.true.` if EOS has a first-order phase transition; enables segmented spline |
| `num_tab` | `integer` | (read from file) | Number of EOS table rows (set by `loadEos`) |
| `n_PT` | `integer` | 0 | Number of detected phase-transition pressure jumps |

**Common values for `eos_file`**: `ALF2`, `APR1`, `MPA1`, `PS`, `PS_G3`,
`PSt`, `Hadronic_1_S1`, `hybrid_1_S1`, `test`

---

## Grid Configuration

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `res` | `integer` | `400` | Resolution base; `SDIV = MDIV = 2*res + 1 = 801` |
| `SDIV` | `integer` | `801` | Radial grid points |
| `MDIV` | `integer` | `801` | Angular grid points |
| `s_pwr` | `integer, parameter` | `1` | Radial mapping power law exponent |
| `SMAX` | `real(wp)` | `1.0 - 0.1^(9/s_pwr)` | Compactification upper limit |
| `angular_collocation` | `integer` | `COLLOCATION_UNI` (1) | Angular grid type |

**Angular collocation constants**:
| Constant | Value | Grid type |
|---|---|---|
| `COLLOCATION_UNI` | 1 | Uniform spacing |
| `COLLOCATION_LEG` | 2 | Legendre Gauss-Lobatto |
| `COLLOCATION_CHEB` | 3 | Chebyshev extrema |

---

## Theory Selection

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `active_theory` | `integer` | `THEORY_GR` (0) | Active theory; 0 = GR, 1 = Scalar-Tensor |

**Theory constants**:
| Constant | Value |
|---|---|
| `THEORY_GR` | 0 |
| `THEORY_ST` | 1 |

Call `initialize_theory()` → selects `apply_gr_defaults()` or
`apply_st_defaults()` accordingly.

---

## Rotation Law

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `solver_type` | `character(20)` | `"uniform"` | Rotation profile: `"uniform"`, `"const_j"`, `"uryu"` |
| `A_diff` | `real(wp)` | `0.5` | Differential rotation parameter A (const-J law) |
| `lambda1` | `real(wp)` | `1.5` | Uryu λ₁ |
| `lambda2` | `real(wp)` | `0.3` | Uryu λ₂ |
| `uyru_p` | `integer` | `1` | Uryu power index p |
| `uyru_q` | `integer` | `3` | Uryu power index q |

---

## Convergence & Accuracy

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `ACCURACY` | `real(wp), parameter` | `1.e-5` | Field convergence threshold for iteration |
| *(relaxation stage thresholds)* | in `spin_relaxation_mod` | see below | Stage-transition dif values |

**Relaxation stage thresholds** (in `spin_relaxation_mod.f90`):
| Name | Value | Transition |
|---|---|---|
| `W_PICARD` | `0.7` | Picard damping weight |
| `PICARD_THRESH` | `0.5` | dif ≥ 0.5 → Picard |
| `CHEB_THRESH` | `0.1` | 0.1 ≤ dif < 0.5 → Chebyshev |
| `M_HIST` | `3` | Anderson history depth |
| `DIF_TAIL_THRESH` | `1.e-5` | Aitken activation |
| `ANDERSON_TAIL_THRESH` | `1.e-4` | Anderson cooldown trigger |
| `N_CHEB` | `3` | Min iterations before Chebyshev kicks in |
| `N_ANDERSON` | `5` | Min iterations before Anderson kicks in |
| `ANDERSON_COOLDOWN` | `30` | Iterations to wait after Anderson fail |
| `AITKEN_COOLDOWN` | `3` | Iterations to wait after Aitken fail |

---

## Physical Constants (compile-time, `constants_mod.f90`)

| Constant | Value | Unit | Meaning |
|---|---|---|---|
| `C` | `2.99792458e10` | cm/s | Speed of light |
| `G` | `6.67408e-8` | cm³ g⁻¹ s⁻² | Newton's constant |
| `MSUN` | `1.98847e33` | g | Solar mass |
| `MB` | `1.6749286e-24` | g | Baryon mass |
| `HBAR` | `6.582119569e-16` | eV·s | Reduced Planck |
| `N_SAT` | `2.7e14` | cm⁻³ | Nuclear saturation density |
| `L_UNI` | `1.4769994423016508` | — | Geometric length unit (G M☉/c²) |
| `KAPPA` | `1.e-15 * C²/G` | — | Pressure geometric unit |
| `KSCALE` | `KAPPA * G/C⁴` | — | Dimensionless scaling factor |
| `PI` | `acos(-1.0_wp)` | — | π |

---

## Build Configuration (Makefile)

| Variable | Default | Override | Meaning |
|---|---|---|---|
| `FC` | `gfortran` | `make FC=ifort` | Fortran compiler |
| `MODE` | `Release` | `make MODE=Debug` | Build mode |
| `FFTW_LIBS` | `pkg-config --libs fftw3` | `make FFTW_LIBS=-lfftw3` | FFTW link flags |
| `LIBS` | `-llapack -lblas $(FFTW_LIBS)` | — | All linked libraries |

**Build targets**:
| Target | Effect |
|---|---|
| `make` / `make all` | Release build (O3) |
| `make debug` | Debug build (O0 + checks) |
| `make clean` | Remove `build/` tree |
| `sub.sh` | `make clean && make -j8 && ./build/bin/a.out` |

**Output binary**: `build/bin/a.out`

---

## Output Paths

| Path | Contents |
|---|---|
| `Cont/` | Solution files: `{EOS}_J{spin}_Mb{mass}_B{field}_...dat` |
| `Cont/Omega.dat` | Radial Ω profile from last run |
| `Cont/velocity.dat` | Orbital velocity diagnostic |
| `Res/` | Research/result storage (`1dprofile.dat`, `res.rst`) |
| `Map/` | Profile mapping outputs (for MATLAB post-processing) |
| `logs/b_goal_sweep/` | Per-B_goal run logs from `sweep_b_goal.sh` |
| `build/bin/a.out` | Compiled binary |
| `build/obj/` | Object files |
| `build/mod/` | Fortran `.mod` interface files |
