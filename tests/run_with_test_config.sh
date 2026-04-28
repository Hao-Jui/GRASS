#!/bin/bash
# Temporarily rewrite source-baked config knobs to the values used to
# generate tests/reference/*.ref, run a command (default: ctest), then
# restore the developer's local config.
#
# Reference baseline (recovered from tests/reference/test_gr_uniform.ref):
#   src/core/starting_model_mod.f90 :: e_center = .8e15_wp
#
# All other knobs (active_theory, solver_type, eos_file, run_task, ...)
# are overridden inside each tests/integration/test_*.f90, so we do not
# need to touch para_panel.f90 or MRcurve_mod.f90 here.
#
# Usage:
#   tests/run_with_test_config.sh [ctest-args...]
#   tests/run_with_test_config.sh -- <command...>
#
# Examples:
#   tests/run_with_test_config.sh
#   tests/run_with_test_config.sh -R test_gr_uryu --output-on-failure
#   tests/run_with_test_config.sh -- ./build/tests/test_gr_uniform

set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

STARTING="src/core/starting_model_mod.f90"
PARA="src/para_panel.f90"
PATCHED_FILES=("$STARTING" "$PARA")
BACKUP_DIR="$(mktemp -d -t grass-test-config.XXXXXX)"
RESTORED=0

restore() {
  if [ "$RESTORED" -eq 1 ]; then return; fi
  RESTORED=1
  for f in "${PATCHED_FILES[@]}"; do
    b="$BACKUP_DIR/$(basename "$f")"
    if [ -f "$b" ]; then
      cp "$b" "$f"
      echo "[run_with_test_config] restored $f"
    fi
  done
  rm -rf "$BACKUP_DIR"
}
trap restore EXIT INT TERM

for f in "${PATCHED_FILES[@]}"; do
  cp "$f" "$BACKUP_DIR/$(basename "$f")"
done

# Rewrite e_center in the `case default` branch (line ~37) to the
# reference baseline. The MODE_REGRID branch (line ~31) is left alone
# because no integration test uses MODE_REGRID. Tolerant of any current
# literal value the developer has hand-edited the file to.
python3 - "$STARTING" <<'PY'
import re, sys, pathlib
p = pathlib.Path(sys.argv[1])
src = p.read_text()
matches = list(re.finditer(r"e_center\s*=\s*[^\s!]+_wp", src))
if len(matches) < 2:
    sys.exit(f"run_with_test_config: expected >=2 e_center literals in {p}, found {len(matches)}")
m = matches[1]
new = src[:m.start()] + "e_center = .8e15_wp" + src[m.end():]
if new == src:
    sys.exit("run_with_test_config: rewrite produced no change")
p.write_text(new)
PY
echo "[run_with_test_config] applied test baseline (e_center = .8e15_wp)"

# Force timing=.false. so the profiling early-exit doesn't short-circuit
# convergence (the developer's local config flips this on for benchmarking).
python3 - "$PARA" <<'PY'
import re, sys, pathlib
p = pathlib.Path(sys.argv[1])
src = p.read_text()
new, n = re.subn(r"(timing\s*=\s*)\.true\.", r"\1.false.", src, count=1)
if n == 0:
    # already .false. -- nothing to do
    sys.exit(0)
p.write_text(new)
PY
echo "[run_with_test_config] applied test baseline (timing = .false.)"

cmake --build build >/dev/null 2>&1 || {
  echo "[run_with_test_config] build failed under test config" >&2
  exit 2
}

if [ "${1-}" = "--" ]; then
  shift
  "$@"
  status=$?
else
  ( cd build && ctest "$@" )
  status=$?
fi

exit "$status"
