# BH Toroid Steps 1-3 Implementation Plan

## Decision

Create an inert black-hole toroid subsystem under `src/bh_toroid/` for Nishida & Eriguchi (1994) steps 1-3:
1. Parameters and mode selection.
2. Horizon-aware radial mapping.
3. Nishida-Eriguchi Green kernels.

Do not route production execution to this subsystem yet.

## Drivers

- Preserve existing neutron-star behavior by default.
- Keep horizon-domain assumptions out of star-centered modules.
- Make equations (3.1)-(3.5) independently unit-testable before solver integration.
- Let implementation proceed in parallel with clear file ownership.

## Alternatives

- `src/bh_toroid/`: chosen. Clear ownership boundary and low risk to current solver.
- `src/theory/`: rejected. Current theory modules assume a regular center, `s_gp`, and star-centered metric targets.
- Documentation/prototype only: rejected. The next step needs implementation-ready Fortran and test design.

## Shared Validation Contract

Prefer `src/bh_toroid/bh_toroid_validation_mod.f90` owned by Stream A.

```fortran
type :: validation_result
  integer :: status
  character(len=160) :: message
end type
```

Status constants:
- `VALID_OK = 0`
- `VALID_BAD_MODEL_FAMILY = 1`
- `VALID_BAD_HORIZON = 2`
- `VALID_BAD_RADIAL_ORDER = 3`
- `VALID_BAD_SCALE = 4`
- `VALID_BAD_POLYTROPE = 5`
- `VALID_BAD_ROTATION = 6`
- `VALID_BAD_GRID_SIZE = 7`
- `VALID_BAD_GREEN_ARGS = 8`

Tests assert status codes; messages are diagnostics only. Invalid-state APIs must not use `error stop`.

## Step 1: Parameters And Mode Selection

Files:
- `src/bh_toroid/bh_toroid_validation_mod.f90`
- `src/bh_toroid/bh_toroid_params_mod.f90`
- `src/para_panel.f90`
- `CMakeLists.txt`
- `tests/unit/test_bh_toroid_params.f90`

Plan:
- Add inert constants to `para_mod`: `MODEL_NS = 1`, `MODEL_BH_TOROID = 2`, `model_family = MODEL_NS`.
- Add `type(bh_toroid_params)` with `h0_hat`, `omega_h`, `rin_hat`, `kappa_ratio`, `emax`, `poly_n`, `rotation_A`, `rout_scale`.
- Add `default_bh_toroid_params()`, `validate_bh_toroid_params()`, and `validate_model_family()`.
- Validate `0 < h0_hat < rin_hat < 1`, positive `rout_scale`, `emax`, `kappa_ratio`, `poly_n`, and `rotation_A`.
- Add root CMake source discovery for `src/bh_toroid/*.f90`.

Tests:
- Defaults are valid and leave `model_family == MODEL_NS`.
- A paper-like sample parameter set is valid.
- Each invalid field returns the expected status code.
- No production solver routing is added.

## Step 2: Horizon-Aware Radial Mapping

Files:
- `src/bh_toroid/bh_toroid_radial_map_mod.f90`
- `tests/unit/test_bh_toroid_radial_map.f90`

Plan:
- Do not mutate `s_gp`, `SDIV`, `DS`, or existing grid globals.
- Add:
  - `validate_radial_domain(h0_hat, rin_hat)`
  - `build_rhat_trapezoid_grid(n, h0_hat, rhat, weights)`
  - `classify_rhat(rhat, h0_hat, rin_hat)`
  - `radius_from_rhat(rout, rhat)`
  - `radial_jacobian(rout)`
- Zone constants: `ZONE_HORIZON`, `ZONE_VACUUM_GAP`, `ZONE_TORUS`, `ZONE_OUTSIDE`.
- Trapezoid grid on `[h0_hat, 1]`: `dr = (1 - h0_hat)/(n - 1)`, endpoint weights `dr/2`, interior weights `dr`, total `1 - h0_hat`.

Tests:
- Domain validation rejects bad horizon/order values.
- Grid endpoints are exactly `h0_hat` and `1`.
- Grid is monotone and never below horizon.
- Weights are positive where expected and sum to `1 - h0_hat`.
- `radius_from_rhat` maps horizon and outer edge correctly.
- Classification handles horizon, vacuum gap, torus, and outside zones.

## Step 3: Green Kernels

Files:
- `src/bh_toroid/bh_toroid_green_mod.f90`
- `tests/unit/test_bh_toroid_green.f90`

Plan:
- Implement pure scalar kernels:
  - `validate_green_args(n, r, rp, h0)`
  - `ne_f1_kernel(n, r, rp, h0)` for paper eq. (3.4)
  - `ne_f2_kernel(n, r, rp, h0)` for paper eq. (3.5)
  - `lambda_radial_kernel(...)`
  - `b_radial_kernel(...)`
  - `omega_radial_kernel(...)`
- Transcribe exact formulas from `papers/1994ApJ...427..429N.pdf`, pages 431-432, into equation-numbered comments.
- Specify branch behavior for `r < rp`, `r = rp`, and `r > rp`; equality uses the stable common value.
- No source integration arrays yet.

Tests:
- Include a comment table with formula, inputs, expected value, and arithmetic source before assertions.
- Cover `n = 0, 1, 2`, `r < rp`, `r = rp`, `r > rp`, `h0 > 0`.
- Cover comparable `h0 -> 0` limits.
- Assert branch continuity at `r = rp`.
- Assert finite values near the horizon for `r > h0`.
- Assert horizon-image cancellation only where explicitly derived from the formula.
- Ambiguous transcription blocks implementation until reported upward.

## Parallel Streams

Stream A, bootstrap/params, executor:
- Owns `CMakeLists.txt`, `src/para_panel.f90`, shared validation module, params module/test.
- Creates `src/bh_toroid/` and `BH_TOROID_SOURCES`.
- Registers `test_bh_toroid_params`.

Stream B, radial map, executor:
- Owns radial map module/test.
- Starts after A bootstrap or in a forked workspace against planned paths.

Stream C, Green kernels, executor or architect/executor, high reasoning:
- Owns Green module/test.
- Must include paper equation comments and oracle tables.

Stream D, verifier/test engineer, sequential after merge:
- Owns final `tests/CMakeLists.txt` consistency.
- Checks whether legacy `Makefile` / `tests/Makefile` need entries for `make test-unit`.
- Runs verification and no-op behavior checks.

## Verification

Baseline:
- `cmake --preset debug`
- `cmake --build --preset debug`
- `ctest --preset debug --output-on-failure`

If the full suite fails because of known local timing/profile configuration, D records that and runs:
- `ctest --test-dir build-debug -R 'test_(ad|spline|spectral|brent|eos)$' --output-on-failure`

New tests:
- `ctest --test-dir build-debug -R 'test_bh_toroid' --output-on-failure`

No-op check:
- Confirm no imports/calls from `src/main.f90`, `src/core/grid_mod.f90`, `src/theory/rotation_solver_mod.f90`, or `src/theory/spin_workspace_mod.f90` to BH/toroid modules.

## Acceptance Criteria

- Existing neutron-star default remains `MODEL_NS`.
- No BH solver routing exists yet.
- Existing unit baseline is unchanged.
- New BH/toroid tests pass.
- Invalid inputs are tested via validation status, not process termination.
- APIs document dimensionless conventions and paper equation references.
- No new dependencies.
- No behavioral changes to `main`, `make_grid`, `rotation_solver`, or `spin_workspace`.
