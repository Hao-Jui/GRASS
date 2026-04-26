# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**For physics, features, and how to use GRASS, see [README.md](README.md).**

## Code Style

- **2-space indentation** everywhere (strictly enforced)
- **`block` constructs** only when introducing local variables for scoping; never for grouping procedure calls
- **Preserve variable names** — never refactor to different naming conventions
- **Do not add comments** to `para_mod.f90` or parameter files that change frequently

## Module Organization & Dependencies

### Core hierarchy
- `para_mod.f90` — central: all parameters, grid arrays, field allocations (no dependencies on other modules except `precision_mod`)
- `eos.f90`, `grid.f90` — depend on `para_mod`; must be initialized before solver
- `constraint_eq.f90`, `analysis_v2.f90` — depend on `para_mod` and utilities; used for diagnostics

### Solver entry point
- `starting_model.f90`, `shoot_v2.f90` — initialize and run one of three tasks (OneModel, shoot, MRcurve)
- These call into the rotation solver (`rotation_solver.f90`) which orchestrates the main loop

### Theory/rotation solver chain
```
rotation_solver.f90  (outer loop on r_e)
  ├─ spin_integration.f90  (compute Green's function targets)
  ├─ spin_relaxation.f90 + relaxation_mod.f90  (4-stage: Picard → Chebyshev → Anderson → Aitken)
  └─ spin_updates.f90  (apply targets to fields)
     └─ eos.f90, rotational_law_mod.f90  (EOS lookups, rotation law)
```

### Utilities
- `toolkit_mod.f90` — interpolation, derivatives, integral accumulation (used throughout)
- `ad_mod.f90` — automatic differentiation for dP/dE (used in EOS and constraint)
- `nag_compat_mod.f90`, `cheb_mod.f90`, `spectral_hub.f90` — numerical bases and quadrature

## Key Design Patterns

**Global state in para_mod**: All grid arrays, field allocations, and parameters live here. This is by design — simplifies passing data through deeply nested solver routines without modifying interfaces.

**Workspace allocation**: `spin_workspace.f90` pre-allocates geometry caches and workspace arrays once per grid; reused across solver iterations for performance.

**Caching strategy**: Bessel function tables and geometry pre-computations cached per iteration. Derivative arrays recomputed each iteration (balances memory vs. compute).

**Log-space EOS**: All EOS quantities stored and interpolated in log-space internally. Conversions to/from physical units happen at lookup time.

**No explicit module cleanup**: Fortran automatic deallocation on program exit; no manual cleanup routines.

## Common Tasks

| Task | File | Notes |
|------|------|-------|
| Change theory (GR ↔ ST) | `para_mod.f90` line 7 | `active_theory = THEORY_GR` or `THEORY_ST` |
| Change resolution | `para_mod.f90` lines 34–37 | `res`, `MDIV`, `s_pwr` |
| Change EOS table | `para_mod.f90` line 40 | `eos_file = "APR1"` etc. |
| Add rotation law | `rotational_law_mod.f90` | Update `solver_type` selector in `rotation_solver.f90` |
| Tune relaxation | `relaxation_mod.f90`, `spin_relaxation.f90` | Adjust `PICARD_THRESH`, `CHEB_THRESH`, `w_picard`, `N_MONOTONE` |
| Add diagnostic | `analysis_v2.f90` | Use `integrate_profiles()` from `toolkit_mod.f90` |
| Modify field update | `spin_updates.f90` | Called during relaxation; updates metric/scalar fields |

## Important Caveats

- **Adding modules**: Add new source files to the appropriate `src/` subdirectory; CMake picks them up via `file(GLOB)`. Re-run `cmake -B build` after adding files.
- **Grid indices**: SDIV (radial) × MDIV (angular); Fortran column-major storage.
- **Fortran 2003 required**: No newer features; explicit interfaces used for BLAS/LAPACK.
- **No parallelization**: Single-threaded; no OpenMP or MPI. Not planned for this codebase.

## Testing & Validation

CMake-driven test suite under `tests/` (unit + integration). Run via `cmake --build build --target <test>` then `./build/tests/<test>`. Reference outputs in `tests/reference/*.ref`. CI runs them on push.

Manual validation:
- Hamiltonian constraint norm (should be O(1e-10) or smaller at convergence)
- M–r curves against published sequences
- Diagnostic files (`Cont/Omega.dat`, `Cont/velocity.dat`) for physical bounds
- `Res/res.rst` round-trip via `test_restart` (writer ↔ reader unit consistency)

## Restart Binary (`Res/res.rst`)

Stream-IO binary, magic `GRASSRST01`, version 1. Writer: `spin_integration_mod.write_restart_file`. Reader: `regrid_mod.read_binary_restart`.

**Unit convention (do not break):** `header_meta = [r_e, energy(1,1), r_ratio, Omega_e, Omega_c]` stored as-is in code units. `Omega_e/Omega_c` are post-normalization (already `/r_e_new` per `rotation_solver_mod.f90:159-160`). Reader reads identity — NO `*r_e`, NO `/(C*C*KSCALE)`. Asymmetric transforms = silent corruption of rotating restarts.

`MODE_REGRID` currently restores fields + grid only; `starting_model_mod` overrides `e_center` and `single_model` resets `r_ratio = 1.0_wp`. Full state-restore not yet supported.

## MCP Usage Policy

Three MCPs registered. Route by question type, not habit.

### codebase-memory-mcp (CMM) — DEFAULT for any code question

Use INSTEAD OF `grep`/`Read` for code discovery. Index is graph (583 nodes, 626 edges at last build).

| Question | Tool |
|---|---|
| "where is X defined / who calls X" | `search_graph(name_pattern=".*X.*")` |
| natural-language discovery | `search_graph(query="...")` (BM25 + structural boost) |
| read full module body | `get_code_snippet(qualified_name)` — NOT `Read` |
| text-in-code | `search_code(pattern, path_filter="^src/", mode="compact")` — filters legacy out |
| call chain / impact | `trace_path(function_name, mode=calls\|data_flow)` |
| multi-hop / aggregations | `query_graph(cypher)` |
| project map | `get_architecture(aspects=['all'])` |
| re-index after edits | `detect_changes` (delta) or `index_repository` (full) |
| persist invariants | `manage_adr(mode='update', content=...)` |

Always pass `path_filter="^src/"` when searching active code — `Legacy_branches/` and `*_backup_*` produce noise (e.g. `restart` matches in 4 legacy branches).

### context-mode (ctx) — only when raw output would flood context

Use for: build/test logs >20 lines, web docs, large bash dumps, data analysis. NOT for code discovery (CMM owns that).

- `ctx_batch_execute(commands, queries)` — multi-step build+test+grep, FTS5-indexed.
- `ctx_execute(language="js"|"shell", code)` — analyze/filter/diff data without piping into context.
- `ctx_fetch_and_index(url)` — Fortran/CMake/BLAS docs.

### claude-mem — cross-session memory

- `mem-search "<topic>"` at session start before re-investigating.
- `get_observations` for prior decisions on this codebase.
- Save observations with `file:line` anchors (e.g. "regrid_mod.f90:130 reader unit transform").

### Order of operations (bug fix / investigation)

1. `mem-search` — prior fix?
2. `detect_changes` since last index — re-sync graph if stale.
3. CMM (`search_graph` → `get_code_snippet` → `trace_path`) — locate + understand.
4. Edit → `cmake --build build --target <test>` → run.
5. ctx for log analysis if output floods.
6. `manage_adr(mode='update')` if a new invariant emerged.

### Anti-patterns

- `grep -r` across repo — use `search_code` with `path_filter`.
- `Read` of file just to scan — use `get_code_snippet` (returns metadata + source).
- `cat Legacy_branches/...` — graph filtering already excludes; re-introducing legacy noise wastes context.
- Skipping `manage_adr` after fixing a unit/convention bug — next session will repeat the investigation.
