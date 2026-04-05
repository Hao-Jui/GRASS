# Theory modules — dependency chain

```
spin_derivatives_mod    spin_updates_mod
  (finite differences)    (r_e, Omega, EOS+velocity)
        │                       │
        ▼                       │
spin_workspace_mod                  │
  (shared arrays,               │
   alloc/dealloc)               │
        │                       │
        ▼                       │
spin_integration_mod ───────────────┘
  (spectral pipeline:
   precompute → build_source → project_integrate → reconstruct,
   alpha potential, output)
        │
        ▼
spin_relaxation_mod
  (Picard / Chebyshev / Anderson mixing,
   Aitken acceleration, cycle detection)
        │
        ▼
rotation_solver
  (outer iteration loop,
   orchestrates all of the above)
```

## File summary

| File | Module | Responsibility |
|---|---|---|
| `spin_derivatives_mod.f90` | `spin_derivatives_mod` | Grid finite-difference operators (`deriv_s`, `deriv_m`) |
| `spin_updates_mod.f90` | `spin_updates_mod` | Equatorial radius balance, angular velocity solve, EOS+velocity field update |
| `spin_workspace_mod.f90` | `spin_workspace_mod` | Shared workspace arrays (derivative caches, source/coefficient/target arrays, geometry factors, quadrature weights), allocation and deallocation |
| `spin_integration_mod.f90` | `spin_integration_mod` | Spectral decomposition pipeline: `precompute` (derivatives, Bessels, exp-tables) → `build_source_terms` (Einstein sources) → `project_integrate` (tiled DGEMM angular projection fused with Green's function radial integration via `green_rho`/`green_multipole`/`green_bessel`) → `reconstruct` (sum coefficients back to 2-D grid). Also: `update_alpha_potential` (lapse), `output_helper` (file I/O). Helper predicates: `is_massive_scalar()`, `is_spherical()`. Workspace accessed directly via `use spin_workspace_mod` (no argument passing). |
| `spin_relaxation_mod.f90` | `spin_relaxation_mod` | Fixed-point relaxation iteration: tier selection (Picard/Chebyshev SOR/Anderson), Aitken delta-squared acceleration, limit-cycle detection, divergence guard, scalar field floor |
| `rotation_solver_mod.f90` | `rotation_solver_mod` | Outer driver: allocates workspace, iterates (update r_e → solve Omega → update EOS → get targets → relax → update alpha) until convergence, finalises bulk properties and output |
