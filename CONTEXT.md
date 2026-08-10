# GRASS

GRASS solves stationary axisymmetric rotating compact-object configurations — neutron stars in General Relativity (GR) or scalar-tensor theory (ST), and a black hole surrounded by a fluid toroid — using an iterative Green's-function field assembly on a Komatsu–Eriguchi–Hachisu (KEH) compactified grid.

## Language

### Theory selector

**GR**:
General Relativity branch of the field equations.
_Avoid_: Einstein, classical gravity.

**ST**:
Scalar-tensor branch (Damour–Esposito-Farèse class), adds a scalar field on top of GR.
_Avoid_: alternative gravity, modified gravity.

**Active theory**:
Compile/runtime switch (`active_theory` in `para_mod.f90`) selecting **GR** or **ST**.

### Configuration families

**Rotating star**:
Single self-gravitating fluid configuration with axisymmetric rotation. The default solver target.
_Avoid_: NS (ambiguous — covers GR-only or ST), star (too generic).

**BH-toroid**:
Black hole at the origin surrounded by a stationary self-gravitating fluid torus. Distinct solver family with a horizon boundary condition.
_Avoid_: BH+disk, accretion disk, ring.

**Single model**:
One equilibrium configuration at fixed central energy density and axis ratio.
_Avoid_: equilibrium, point.

**MR curve**:
Sequence of single models sweeping a parameter (typically central energy density) producing a mass–radius curve.
_Avoid_: sweep, scan, sequence.

**Shoot**:
Bisection-style task that adjusts one input until a target observable is hit.
_Avoid_: search, root-find.

### Solver pipeline

**Spin integration**:
Step that evaluates the Green's-function integrals on the current fields to produce **targets**.
_Avoid_: integration step, KEH step.

**Targets**:
Right-hand-side field values (`γ̂`, `ν̂`, `ω̂`, plus toroid analogues) computed from the current state. The relaxation step drives the live fields toward these.
_Avoid_: RHS, sources, source terms.

**Spin updates**:
Step that turns **targets** plus EOS / **rotation law** into updated fluid + metric fields.
_Avoid_: hydro step, field update.

**Spin relaxation**:
Iterative scheme that applies **targets** to fields with damping/acceleration. Four stages: Picard → Chebyshev → Anderson → Aitken.
_Avoid_: iteration, fixed-point loop, mixer.

**Rotation law**:
Functional form for specific angular momentum vs. radius (e.g. j-constant, KEH-style). Plug-in in `rotational_law_mod.f90`.
_Avoid_: spin profile, Ω(r) law.

**Horizon target**:
Boundary condition imposed at the inner edge of a **BH-toroid** configuration, of form `Ω_h / r_out²` (see `bh_toroid/solver_mod.f90`).
_Avoid_: inner BC, BH boundary.

**Zone classification**:
Partition of the **BH-toroid** radial domain into regions (`classify_rhat` in `bh_toroid/core_mod.f90`).
_Avoid_: region split, segmentation.

### Grid + state

**Grid**:
Two-dimensional discretization with `SDIV` radial points (compactified `s`-coordinate) by `MDIV` angular points (`μ = cos θ`). Defined in `para_mod.f90`.
_Avoid_: mesh, lattice.

**Fields**:
Set of metric (`alpha`, `gama`, `rho`, `omg`) and fluid (`energy`, `pressure`, `enthalpy`) arrays carried on the **grid**.
_Avoid_: variables, state arrays.

**r_e**:
Equatorial coordinate radius of the surface. Primary length scale for nondimensionalization.
_Avoid_: R, surface radius.

**r_ratio**:
Polar-to-equatorial radius ratio `r_p / r_e`. Controls rotational deformation. `1.0` ⇒ spherical.
_Avoid_: axis ratio, oblateness.

**Ω_e**, **Ω_c**:
Angular velocity at the equator and at the center, post-normalized by `r_e_new` per `rotation_solver_mod.f90:159–160`.
_Avoid_: equatorial spin, central spin.

**EOS**:
Tabulated equation of state interpolated in log-space internally. Selected via `eos_file` in `para_mod.f90`.
_Avoid_: equation of state (spell out only on first use), pressure law.

**Restart**:
Stream-IO binary `Res/res.rst` with magic `GRASSRST01`. Stores `header_meta` plus packed **fields** for resumption / regridding.
_Avoid_: checkpoint, snapshot.

**Header meta**:
Five-element vector `[r_e, energy(1,1), r_ratio, Ω_e, Ω_c]` at the head of a **restart** file. Stored as identity (no unit transform on either side).
_Avoid_: header, metadata block.

### Math conventions

**Code units**:
Dimensionless units used internally: lengths in units of `r_e`, energies/pressures normalized by `KSCALE`, speeds by `C`.
_Avoid_: natural units, normalized units.

**Hamiltonian constraint**:
Residual norm of the Hamiltonian constraint equation; used as primary convergence diagnostic. Target O(1e-10).
_Avoid_: constraint residual, constraint error.

## Relationships

- An **active theory** (**GR** or **ST**) is selected once per build/run and shapes the field equations every **spin integration** evaluates.
- A **rotating star** is a **single model**, an **MR curve**, or a **shoot** target. A **BH-toroid** is currently always a **single model**.
- One iteration of the solver runs **spin integration** → produces **targets** → **spin relaxation** drives **fields** toward those **targets** via **spin updates**.
- **Spin updates** consults the **EOS** and the **rotation law** to convert metric updates into fluid quantities.
- A **BH-toroid** configuration adds a **horizon target** boundary condition and a **zone classification** of the radial **grid**, on top of the standard pipeline.
- A **restart** captures **fields** plus a **header meta**; the writer (`spin_integration_mod`) and reader (`regrid_mod`) must agree on the **header meta** layout and on the unit convention (identity, in **code units**).
- **r_e**, **r_ratio**, **Ω_e**, **Ω_c** are the four scalars that together identify a **rotating star** equilibrium up to its **EOS** and **active theory**.

## Example dialogue

> **Dev**: "When a **BH-toroid** run starts from a **restart**, does the solver re-run **spin integration** before the first **spin update**?"
> **Domain**: "Yes — the **restart** restores **fields** and grid only. **Targets** are recomputed on the first iteration because they're never persisted. Same path as a fresh start, only the initial guess differs."

> **Dev**: "If I want **Ω_e** out of the **header meta**, do I divide by `r_e`?"
> **Domain**: "No. **Ω_e** in the **header meta** is already post-normalized in **code units**. Reader is identity. Multiplying by `r_e` would silently corrupt it."

## Flagged ambiguities

- "module" — overloaded. In Fortran-source context refers to a `module ... end module` unit (e.g. `spin_relaxation_mod`). In architecture-review context (`/improve-codebase-architecture`, `LANGUAGE.md`) refers to anything-with-an-interface, scale-agnostic. Both senses appear in CLAUDE.md and ADRs.
- "spin" — historically referred to neutron-star spin specifically. Now used as a prefix for the whole pipeline (`spin_integration`, `spin_updates`, `spin_relaxation`) including BH-toroid runs. Treat as a pipeline-stage prefix, not a physics term.
- "relaxation" — overloaded between **spin relaxation** (the 4-stage outer iteration) and the JFNK/GMRES inner solver in `relaxation_mod.f90`. The latter is the "relaxation kernel"; the former is the "relaxation scheme."
- "core" — appears as both `src/core/` (codebase area: para_mod, regrid, starting_model) and `bh_toroid/core_mod.f90` (validators + grid builders for the BH-toroid family). Rename one if confusion bites.
