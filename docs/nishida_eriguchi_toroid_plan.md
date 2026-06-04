# Nishida & Eriguchi (1994) Black-Hole Toroid Implementation Plan

Source paper: `papers/1994ApJ...427..429N.pdf`, "A General Relativistic Toroid Around a Black Hole", Nishida & Eriguchi, ApJ 427:429-437.

## Paper Algorithm

The paper extends the Komatsu-Eriguchi-Hachisu integral-field method to stationary, axisymmetric, equatorially symmetric, asymptotically flat black-hole plus toroid equilibria.

Core variables:
- Metric: `ds^2 = -exp(2 nu) dt^2 + exp(2 alpha)(dr^2 + r^2 dtheta^2) + exp(2(gamma-nu)) r^2 sin^2(theta)(dphi - omega dt)^2`
- Horizon boundary at `r = h0`: `B = exp(gamma) = 0`, `exp(nu) = 0`, `omega = omega_h`.
- Fluid: perfect-fluid polytropic toroid, circular rotation.
- Integral equations solve for `lambda = exp(nu)`, `B`, and `omega`; alpha is solved from the same first-order constraint family as earlier KEH-style work.
- Green kernels include image terms in `h0` so the event-horizon boundary conditions are satisfied directly.
- Seven fixed model parameters define one equilibrium: `h0/rout`, `omega_h`, `rin/rout`, `pmax/epsilon_max`, `epsilon_max`, polytropic index `N`, and rotation parameter `A`.

## Fit To Current GRASS Code

The current code already has reusable pieces:
- `src/theory/rotation_solver_mod.f90`: Picard loop, relaxation, metric update sequencing.
- `src/theory/spin_integration_mod.f90`: source construction, angular projection, radial Green-function integration, output moments.
- `src/theory/spin_workspace_mod.f90`: quadrature weights, Legendre bases, target arrays.
- `src/theory/spin_updates_mod.f90`: hydrostatic update, rotation-law roots, EOS inversion.
- `src/tool/brent_mod.f90`: robust scalar root finding for angular velocity laws.
- `src/core/grid_mod.f90` and `src/para_panel.f90`: radial/angular grid, geometry arrays, global model parameters.

The implementation should be a new solver family, not a small option inside `rotation_solver`, because the radial domain starts at a horizon rather than the stellar center and the solved metric variables differ (`lambda`, `B`, `omega` rather than the current star-centered `rho`, `gama`, `ww` targets).

## Proposed Implementation Slices

1. Add black-hole toroid parameters and mode selection.
   - Extend `para_mod` with `model_family`, `bh_h0_ratio`, `bh_omega_h`, `torus_rin_ratio`, `torus_kappa_ratio`, `torus_emax`, `torus_poly_n`, and `torus_A`.
   - Keep defaults inert so current neutron-star tests are unchanged.

2. Add horizon-aware radial mapping.
   - Current radial coordinate maps a stellar center to infinity via `s/(1-s)`.
   - Add a toroid coordinate `r = rout * rhat`, with `rhat in [h0_hat, 1]` for matter support and an exterior extension for asymptotic diagnostics.
   - Do not reuse center regularity assumptions from `sphere_mod` or `rotation_solver_mod`.

3. Add Nishida-Eriguchi Green kernels.
   - Create `src/theory/bh_toroid_green_mod.f90`.
   - Implement equations equivalent to paper (3.1)-(3.5): even Legendre kernel for `lambda`, odd/sine kernel for `B`, associated-Legendre kernel for `omega`, with the `h0` image terms.
   - Unit-test limiting behavior as `h0 -> 0` against the existing regular-origin kernels where comparable.

4. Add toroid hydro/rotation update.
   - Create `src/theory/bh_toroid_updates_mod.f90`.
   - Implement the integrated Euler equation from paper (3.15), the rotation law from (3.16), and `v` from (3.17).
   - Solve constants using the three equatorial points `H(h0, pi/2)`, `S(rin, pi/2)`, and `T(rout, pi/2)`.

5. Add solver loop.
   - Create `src/theory/bh_toroid_solver_mod.f90`.
   - Reuse `spin_relaxation_mod` only after wrapping targets into compatible arrays, or factor relaxation into a metric-agnostic helper.
   - Outputs should include toroid mass/angular momentum from paper (4.1)-(4.2), plus black-hole horizon quantities.

6. Add tests.
   - Unit tests for Green kernels, rotation-law residuals, and hydro boundary conditions.
   - A low-resolution integration fixture that verifies convergence, finite fields, zero pressure at `rin/rout`, `omega = omega_h` on the horizon, and asymptotic flatness.

## Main Risks

- The paper is polytropic, while current production EOS support is tabulated. Start with a polytropic test EOS path before wiring tabulated EOS.
- Current globals assume star-centered quantities such as `h_center`, `r_ratio`, and `r_e`; black-hole toroids need separate semantics to avoid corrupting neutron-star workflows.
- Existing moment extraction assumes exterior stellar multipoles on the current compactified grid; black-hole horizon surface terms need separate accounting.
