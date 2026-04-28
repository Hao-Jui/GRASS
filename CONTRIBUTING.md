# Contributing to GRASS

Thanks for your interest. GRASS is a small, research-driven Fortran
codebase; the workflow is light. Read this once before opening your
first PR.

## Getting started

```bash
git clone https://github.com/Hao-Jui/GRASS.git
cd GRASS
cmake --preset release
cmake --build build -j
ctest --test-dir build
```

If the test suite is green on your machine, you are set up correctly.
See [TESTING.md](TESTING.md) for the full harness inventory.

## Branch conventions

- `main` — released, citable, reproducible. Direct pushes are
  prohibited; merges land via PR from `dev`.
- `dev` — active development. Default branch for new work.
- `feat/<short-name>`, `fix/<short-name>`, `perf/<short-name>` — topic
  branches off `dev`. Keep them short-lived (< 1 week ideally) and
  rebase onto `dev` before opening a PR.

Never force-push `main`. Force-pushing your own topic branch before it
is merged is fine.

## Commit conventions

Follow Conventional Commits, lowercase type:

```
<type>(<scope>): <short summary in imperative mood>

<body — why, not what; reference issues / commits as needed>
```

Types we use: `feat`, `fix`, `perf`, `refactor`, `chore`, `docs`,
`test`, `build`, `ci`. Subject ≤ 72 characters, body wrap at 72.
Example:

```
perf(eos): switch interpolation to PCHIP, add test-config harness

PCHIP gives O(1) lookup after binary search, exact at table nodes,
and naturally degrades to linear at flat phase-transition steps.
Net -29 lines vs the prior file. tests/run_with_test_config.sh
keeps integration regression alive without reverting local config.
```

Do **not** mix unrelated changes in one commit. Local run-config edits
(`e_center`, `eos_file`, `THEORY_*`, output paths) are not commits;
keep them as uncommitted working-tree changes or stash them.

## Code style

- 2-space indentation; tabs are forbidden.
- Fortran 2003 only. No newer-than-2003 features (no submodules
  requiring 2008+, no coarrays, no `do concurrent` in committed code).
- `block` constructs only when introducing local scope; never as
  visual grouping for procedure calls.
- Preserve existing variable names; do not rename across patches.
- Do not add comments to `src/para_panel.f90` or other parameter files
  that change between runs.
- Default to writing no comments; add one only when the *why* is
  non-obvious. Never narrate *what* the code does.
- Public API of a module goes in a single `public ::` block at the top.

## Tests

- Every behaviour change must come with a passing test, an updated
  reference output, or a documented justification for neither.
- Unit tests live under `tests/unit/` and use `tests/unit/test_utils.f90`
  for assertions.
- Integration tests live under `tests/integration/`, set their own
  configuration in-program, and are gated by
  `tests/check_regression.sh` (1e-4 relative tolerance) against
  `tests/reference/test_<name>.ref`.
- If you change the solver and the reference output legitimately
  shifts, regenerate it deterministically (release build, default
  thread count) and include the regen rationale in the PR body.

Run before opening a PR:

```bash
bash tests/run_with_test_config.sh
```

This wraps `ctest` with the source-baked baseline config so
integration tests pass even when your working tree has local sweep
config.

## Pull-request process

1. Fork or branch from `dev`.
2. Keep the PR focused — one feature, one fix, one perf win.
3. Update `README.md` / `TESTING.md` / module docstrings if you change
   user-visible behaviour.
4. CI must be green (unit + integration + regression + sanitizer).
5. Maintainer reviews on a one-pass cadence; expect rebase requests if
   `dev` has moved.
6. Squash-merge to `dev`. `dev` → `main` happens at release time.

## Reporting bugs

Open a GitHub issue with:

- compiler version + OS + build flags,
- minimal reproducer (config edits + a one-line invocation),
- the `Cont/properties.dat` excerpt or console log,
- the EOS file used (or its first 5 rows).

## Releases

Releases are tagged on `main`. CI must be green on the merge commit;
no separate release branch.

---

Questions that aren't bug reports are welcome as GitHub Discussions or
direct email to the maintainer. Keep PRs for code.
