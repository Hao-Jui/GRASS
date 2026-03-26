# Theory modules — dependency chain

```
spin_derivatives        spin_updates
  (finite differences)    (r_e, Omega, EOS+velocity)
        │                       │
        ▼                       │
spin_workspace                  │
  (shared arrays,               │
   alloc/dealloc)               │
        │                       │
        ▼                       │
spin_integration ───────────────┘
  (spectral pipeline:
   source → angular → radial → targets,
   alpha potential, output)
        │
        ▼
spin_relaxation
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
| `spin_derivatives.f90` | `spin_derivatives` | Grid finite-difference operators (`deriv_s`, `deriv_m`) |
| `spin_updates.f90` | `spin_updates` | Equatorial radius balance, angular velocity solve, EOS+velocity field update |
| `spin_workspace.f90` | `spin_workspace` | Shared workspace arrays (derivative caches, source/coefficient/target arrays, geometry factors, quadrature weights), allocation and deallocation |
| `spin_integration.f90` | `spin_integration` | Spectral decomposition pipeline: precompute derivatives+Bessels → build Einstein source terms → angular integration (DGEMM) → radial integration (Green's functions) → reconstruct 2D targets. Also: lapse potential (`update_alpha_potential`), file output (`output_helper`) |
| `spin_relaxation.f90` | `spin_relaxation` | Fixed-point relaxation iteration: tier selection (Picard/Chebyshev SOR/Anderson), Aitken delta-squared acceleration, limit-cycle detection, divergence guard, scalar field floor |
| `rotation_solver.f90` | `rotation_uniform` | Outer driver: allocates workspace, iterates (update r_e → solve Omega → update EOS → get targets → relax → update alpha) until convergence, finalises bulk properties and output |
