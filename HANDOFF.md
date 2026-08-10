# Development handoff

## 2026-08-09 — DD2/Uryu feedback audit

Scope was investigation only. No production source or test was changed.
`src/para_panel.f90` remains a tracked local-configuration change and was
excluded from this work: comparing the compiler's symbol trees for `HEAD` and
the working copy produced `schema_symbol_diff_exit=0`.

### Confirmed: Uryu coefficient-domain failure

`src/theory/rotational_law_mod.f90` evaluates the raw equatorial momentum in
both `cache_uryu_ab` (lines 12-16) and `diff_rotation_uryu` (lines 39-59).
`AA_h` line 98 then raises a negative base to `1/(p+q)`.

For the repository defaults (`lambda1=1.5`, `lambda2=0.3`, `p=1`, `q=3`), a
direct call through the current module gave:

```text
threshold=  1.709976E+00
F_e=  1.0101E-01 denominator=  7.5051E-01 base= -2.0184E-01 finite=F
F_e=  2.2222E+00 denominator=  1.8111E+00 base=  2.1990E+00 finite=T
```

This isolates the cause from a zero denominator and from Brent itself. The
published real-domain condition is
`F_e > (lambda1/lambda2)^(1/q) F_max` (Iosif & Stergioulas 2021, equations
23-25): https://academic.oup.com/mnras/article/503/1/850/6133452

The existing integration executable also reproduced the downstream symptom:

```text
$ GRASS_TIMING=0 GRASS_E_CENTER=0.8e15 ./build/tests/test_gr_uryu
Note: The following floating-point exceptions are signalling: IEEE_INVALID_FLAG
ERROR STOP brent_core: failed to bracket root
check  ./Cont/diff_rotation.dat
```

The resulting bracket log had 200 rows: 141 finite and 59 NaN.

The feedback's proposed `safe_F_e` clamp is not accepted as a verified fix.
Clamping only the coefficient inversion substitutes a different state and can
hide an inadmissible root. A fix needs domain-aware bracketing/parameterisation,
explicit rejection of invalid trials, and post-root validation before caching.

Falsifier checked: a finite `AA_h` below the published threshold. It was false.

### CANNOT REPRODUCE: claimed intrinsic branch discontinuity

The current checkout stops on the confirmed domain failure first. The feedback's
diagnostic instrumentation and continuation experiments are not present, and no
raw trace was supplied. Remaining candidates are multiple-root switching,
crossing the admissibility boundary, fixed-step `Fmax_h` overshoot/non-convergence,
and a physical branch endpoint. The uncapped fixed-step loop is at
`src/theory/spin_updates_mod.f90:166-219`; negative `F_equator_h` is checked only
after root selection at line 174.

Required falsifier was not checked: enumerate all admissible roots versus
`Fmax_h` and demonstrate that the continued physical branch itself terminates
discontinuously. A smaller-step attempt failing sooner does not distinguish the
live candidates.

### Feedback/check-out inconsistencies

- No `safe_F_e` implementation exists.
- The `r_ratio=0.7` overwrite in `single_model` is removed locally, but
  `MODE_DEFAULT` still overwrites `r_ratio` in
  `src/core/starting_model_mod.f90:34`, and restart loading still overwrites it
  in `src/core/regrid_mod.f90:69`.
- The Uryu integration test does not set `r_ratio`, `lambda1`, or `lambda2`; its
  reference still expects the removed hard-coded ratio `0.7`.
- The feedback's DD2 Hamiltonian-L2 numbers and constraint-field instrumentation
  have no preserved outputs in this checkout and were not independently verified.

### Adjacent confirmed issue from the same feedback

The DD2 TOV run can stall without a progress/iteration cap in
`src/core/sphere_mod.f90:142`. An isolated run produced 51,780,960 repeated
non-progress rows (2.38 GB) before it was stopped; the temporary reproduction
was removed. This was not an Uryu-law verdict and no fix was made.

## 2026-08-09 — Uryu admissible-domain fix and rotation smoke tests

Scope was the confirmed Uryu coefficient-domain failure and smoke coverage for
the small Uryu/constant-J rotation helpers. `src/para_panel.f90` was not edited
or included: the coefficient/root-search change does not alter its schema.

### Implementation

- `uryu_coefficients` now checks the strict published boundary
  `F_e > F_max*(lambda1/lambda2)^(1/q)`, parameter ranges, denominators, power
  bases, and finite results before returning `A` and `B`.
- `cache_uryu_ab` refuses an invalid final state instead of retaining or
  clamping it. Residual callbacks return quiet NaN for invalid trial states.
- Uryu equatorial roots use `find_omega_e_admissible`, which searches only the
  physical velocity interval, does not bracket across a non-finite gap, and
  reports no-root/non-convergence through an optional status for tests.
- Uryu local roots use the bounded, guess-centred
  `zbrent_rot_admissible`. Constant-J retains the original `find_omega_e` and
  `zbrent_rot` implementations exactly.
- `test_rotation_law` exercises 43 success, boundary, invalid/refusal, and
  disconnected-domain checks across the constant-J and Uryu coefficient,
  cache, residual, integral, context, equatorial-root, and local-root helpers.

The first implementation put both laws on the bounded root routines. That was
refuted by `test_gr_constj`: the patched path reached the 2,000-iteration stop,
while the same checkout linked against the `HEAD` Brent implementation
converged in 0.1418 s. The final split restores the legacy constant-J path and
uses the admissible routines only at the six Uryu call sites.

### Verification evidence

```text
$ ./build/tests/test_rotation_law
rotation_law_mod: 43 passed, 0 failed

$ ctest --test-dir build --output-on-failure -R '^(test_rotation_law|test_brent)$'
100% tests passed, 0 tests failed out of 2

$ ctest --test-dir build-debug --output-on-failure -R '^(test_rotation_law|test_brent)$'
100% tests passed, 0 tests failed out of 2
```

Direct integration checks on the final split both exited zero and printed
`Converged`:

```text
$ GRASS_TIMING=0 GRASS_E_CENTER=0.8e15 ./build/tests/test_gr_constj
Elapsed time [s]: 0.1605
STOP One model solved!

$ GRASS_TIMING=0 GRASS_E_CENTER=0.8e15 ./build/tests/test_gr_uryu
Fmax   1.128592172E-01   Omega_e   5.880283896E-02
Elapsed time [s]: 21.9531
STOP One model solved!
```

The direct integrations still report IEEE flags already seen in this dirty
checkout. The CTest regression wrappers fail their stored-output comparisons:
both references record the removed hard-coded `r_ratio=0.7`, while the current
dirty configuration produces `0.95` for constant-J and `0.9` for Uryu. This is
separate from the root-domain fix. The final wrapper run reported two failures;
the constant-J comparison was isolated as
`current_axial_ratio=9.500000000E-01`,
`reference_axial_ratio=7.000000000E-01`.

### Performance evidence

An `-O3 -march=native` temporary benchmark (10,000,000 coefficient evaluations;
source removed after the run) compared the former two coefficient functions to
the checked consolidated calculation:

```text
legacy_coefficient_seconds=  0.162750
checked_coefficient_seconds=  0.210142
checked_over_legacy=  1.2912
checked_nanoseconds_per_call=    21.014
```

Thus the isolated kernel is about 29% slower, but the absolute increment is
4.739 ns per coefficient pair. The same benchmark measured the 512-interval
Uryu equatorial search at 87.881 microseconds/root (540 residual evaluations).
These costs are negligible relative to the measured 21.9531 s integration, and
constant-J pays neither bounded-search cost.

## 2026-08-09 — Uryu full-model gate reinforcement

The existing `test_gr_uryu` inherited its central energy, axial ratio, and Uryu
parameters from mutable production defaults while comparing against a fixed
`r_ratio=0.7` reference. Its wrapper also retained stale `Cont/properties.dat`,
hid solver output, and allowed missing reference labels to pass. Labels that
contained digits (for example `M2/M^3`) were parsed from the label rather than
the value, so those diagnostics compared `2` with `2` instead of comparing the
reported multipoles.

The in-progress strengthened gate explicitly selects MPA1, 401x41,
`e_center=0.8e15`, `r_ratio=0.7`, `lambda1=1.5`, `lambda2=0.3`, and
`p/q=1/3`; it bypasses `initialize_starting_model`, removes stale output before
the solve, exposes solver output, serialises the CTest case, and requires every
numeric reference label to be present and within tolerance. Refusal probes
confirmed that an empty properties file and a zero-output status-0 executable
fail. Mutating only `M2/M^3` from `0.697368191` to `9.697368191` produced
`rdiff=1.29e+01` and failed.

The stronger historical-parity gate exposed a dev convergence regression. On
detached local main `13372cd`, with the exact current
`src/para_panel.f90` (SHA-256
`ce92d32815759585eaae94d4b36a6e0bb157fbeffbaeecfd817641177f13b7b2`), the
model converged to the stored solution in 0.85 s wall time (GRASS CPU timer
3.3871 s). The dirty dev solve timed out after 240 s and oscillated through
iteration 300; a clean detached `8498af8` solve also missed 130 s, with
`dif=2.901E-06` at iteration 100. Therefore the 120 s gate failure is not
coefficient-check overhead and could not be resolved honestly by increasing
the timeout.

Detached A/B builds isolated the slowdown to the new equatorial
`find_omega_e_admissible` selection in `spin_updates_mod.f90`: tolerance-scale
root stair-stepping perturbed the repeated `Fmax_h` update. Restoring the
`rotational_law_mod.f90` and `spin_updates_mod.f90` pair to `8ab4661` converged
in 0.48 s; restoring only `spin_updates_mod.f90` failed because legacy Brent
then encountered the checked NaN residuals. The authorised shared-tree
restoration therefore put the complete three-file Uryu root path back at
`8ab4661`:

- `src/theory/rotational_law_mod.f90`
- `src/theory/spin_updates_mod.f90`
- `src/tool/brent_mod.f90`

The removed fix is archived as `feedback/uryu_fix_8498af8.patch` (with a second
copy in `/tmp/grass-uryu-fix.L4eKn6/`); the fix patch SHA-256 is
`53d3ba5b756867bd38863d665854b69bdf29d484c8c4e78e3c99ecfff0d32aab`.
The helper unit test introduced with that API was removed from the build. This
also means the active restored Uryu implementation has no coefficient-domain
or admissible-root checks; it again evaluates raw `AA_h`/`BB_h` expressions.

The strengthened full-model driver converged in 0.77 s on the dirty shared
tree. Strict comparison initially exposed `M4/M^5=-0.125280283` versus reference
`-0.125294142` (`rdiff=1.11E-04`), a label the old checker never compared. A
documented per-field tolerance annotation now permits `1E-3` for Uryu `M4/M^5`,
supported by the `8.99E-4` variation observed across clean restored builds; the
other 15 fields retain `1E-4`. Final CTest verification passed `test_brent` and
`test_gr_uryu` in 0.58 s. A deliberate `5E-4` perturbation to default-tolerance
`M2/M^3` failed, and a malformed tolerance annotation was refused.

## 2026-08-10 — Native-F Uryu attempts refuted by the runtime gate

The requested follow-up was a strict admissibility fix with no more than a 5%
full-model runtime regression. Before editing, the restored legacy executable
was copied to `/tmp/grass-uryu-baseline.Qh0oa9/test_gr_uryu_legacy`; its SHA-256
is `08b99d6aae04fa8f57833106bcf197825fcf9c8c377bb104d5e3cdc9f21e03d2`.
Fifteen runs after three warm-ups had a 0.350 s median (0.350--0.370 s), making
the acceptance ceiling 0.3675 s under the same harness.

Two native-momentum implementations were tried. Both retained the legacy
constant-J path, checked the strict Uryu coefficient domain, inverted
`F -> Omega` analytically, and solved Uryu equatorial/local residuals in `F`.
The inversion identities and strict coefficient boundary were independently
checked by `wolfram/verify_uryu_momentum_inversion.wls`:

```text
$ /Applications/Wolfram.app/Contents/MacOS/wolframscript -file wolfram/verify_uryu_momentum_inversion.wls
URYU_MOMENTUM_INVERSION=PASS
```

Attempt 1 used opposite geometrically expanded endpoints. Its focused tests
passed (`rotation_law_mod: 19 passed, 0 failed`; `brent_mod: 4 passed, 0
failed`), but the full gate refused an apparently empty equatorial domain in
0.28 s. A failing-context trace showed the real cause one iteration earlier:
at `Fmax=1.8582348002331E-4`, the local solve jumped to a distant root with
`Omega_max=4.3990526263928`, driving the next `Fmax` to
`0.3636530649922072`. The next admissible lower bound was
`F=0.621837994071932` while the continuation guess mapped to
`F=0.145506287197738`; a 2,000-point admissible scan found no sign change.

Attempt 2 instead selected the nearest adjacent sign change on either
geometric ray from the previous root. Its 44 focused coefficient, mapping,
boundary, refusal, constant-J, Uryu, and multiple-root branch-selection checks
all passed, as did the four existing Brent checks. Nevertheless, the `.7`
full-model executable was still running after 52 s, over 148 times the legacy
median and already decisively outside the 5% ceiling, so the run was stopped.
This refutes the assumption that mathematically equivalent native-F root
selection preserves the coupled fixed-point/Fmax iteration map.

Per the two-attempt stop rule, both implementations were removed. The active
`rotational_law_mod.f90`, `spin_updates_mod.f90`, and `brent_mod.f90` again
match `8ab4661`, and the rebuilt Uryu executable is byte-identical to the
frozen legacy executable (same SHA-256 above). Thus no runtime or convergence
regression remains in the working implementation, but the original Uryu
admissibility bug also remains open; it must not be reported as fixed.

## 2026-08-10 — EOS regression-data retention and local history rewrite

All seven automated test programs that load an EOS explicitly select `MPA1`:
`test_eos`, `test_gr_uniform`, `test_gr_constj`, `test_gr_uryu`,
`test_st_uniform`, `test_st_uniform_r07`, and `test_restart`. Accordingly,
`.gitignore` now ignores every `*.dat` file except `/eos/MPA1.dat`, and that is
the only tracked EOS table. Eighteen formerly tracked tables were removed from
the index without deleting their working-tree copies. The strongest
counterexample is that legacy manual defaults still name `PS` or `H4`; those
tables are intentionally local-only because the requested retention criterion
was the automated regression suite.

Before rewriting history, tracked staged and unstaged state was saved as binary
patches and all untracked/ignored files were archived under
`/Users/horay/ptmp/GRASS-eos-recovery.WkgHZx`. A complete verified bundle at
`GRASS-before-eos-filter.bundle` has SHA-256
`e9c360d4b49d15182c265b96deddca745bc4b2f1f352d9f51db145ef14aeced2`.
It contains the old local branches, all five pre-rewrite stash commits, and the
Codex checkpoint tree ref that otherwise caused `git-filter-repo` to skip 78
EOS paths. The checkpoint ref was then CAS-deleted before filtering.

The purge manifest contains 80 unwanted paths (79 historically reachable plus
the untracked-only `DD2_hot_equal_pinned.dat`) and has SHA-256
`880a6feaf228d28a4959529dbe4d1ce941ce4347d465c6709d73d9b28a5f8518`.
`git filter-repo --force --paths-from-file ... --invert-paths` rewrote all 177
local commits and the stash reflog. Final reachability checks reported:

```text
REACHABLE_EOS_PATHS_AFTER
eos/MPA1.dat
UNWANTED_REV_OBJECTS_AFTER=0
CHECKPOINT_REF_COUNT=0
```

The rewrite removed `origin`; it was re-added with the original URL but was not
fetched or pushed, so no old remote-tracking refs were made reachable again.
The remote repository itself still has the old history until the owner
explicitly force-pushes the rewritten branches. Four pre-existing user stashes
remain after the temporary rewrite stash was applied and dropped.

The user worktree was restored exactly: SHA-256 values for both the staged
patch (`ec95491027f3b9a8427fa53d81dcf270cf4f261f32580f5b7437f6fa05dd42c1`)
and unstaged patch
(`3774e02202249d54d85e65076616adc54708e2f5d1ae59e4ae16727b80ad354b`)
matched before and after filtering. All 79 physical top-level EOS tables remain
present locally and ignored except for tracked `MPA1.dat`. The retained table
passed its unit test:

```text
$ cmake --build build --target test_eos -j2 && ctest --test-dir build --output-on-failure -R '^test_eos$'
1/1 Test #5: test_eos ... Passed 0.28 sec
100% tests passed, 0 tests failed out of 1
```
