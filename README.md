# GRASS
## General Relativistic Axisymmetric Spacetime Solver

GRASS computes rotating neutron-star equilibria in either General Relativity (GR)
or Scalar–Tensor (ST) gravity. The code integrates the field equations on a
compactified meridional grid, reads tabulated equations of state (EOS),
and outputs diagnostic profiles that characterise the equilibrium model.

---

## Features

- Uniformly rotating, axisymmetric neutron-star configurations.
- GR and ST theories selectable at run time (`gr` or `st`).
- Tabulated EOS support (pressure, energy density, enthalpy tables).
- Compactified grid in meridional coordinates with spectral-type
  interpolation.
- Batch integration helper (`integrate_profiles`) to reduce repeated NAG
  calls.
- Post-processing exports including angular velocity, enthalpy,
  sound-speed squared (`dP/dE`), and angular-momentum flux.

---

## Repository Layout

```
src/
  main.f90             # program entry point
  core/
    para_mod.f90       # global parameters, grid/EOS setup
    toolkit_mod.f90    # interpolation, derivatives, file writers
    shoot_v2.f90       # shooting / relaxation driver
    newton_mod.f90     # Newton iteration solver
    analysis_v2.f90    # global integrals, diagnostics, profiles
    simpson_mod.f90    # numerical integration routines
data/eos/              # expected location for EOS tables (user supplied)
Cont/                  # output directory for diagnostic profiles
matlab/, not_yet_merged_code/  # legacy utilities / archived experiments
```

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
./grass gr   # General Relativity
./grass st   # Scalar–Tensor gravity
```

The program prints the active theory, loads the EOS tables, assembles the
grid, and iterates to a self-consistent solution. Diagnostic files are
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

