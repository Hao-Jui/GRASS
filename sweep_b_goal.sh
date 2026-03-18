#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARA_FILE="$ROOT_DIR/src/core/para_mod.f90"
BUILD_MODE="${MODE:-Release}"
START_VALUE="${1:-3800}"
END_VALUE="${2:-4800}"
STEP_VALUE="${3:-50}"
LOG_DIR="${ROOT_DIR}/logs/b_goal_sweep"

mkdir -p "$LOG_DIR"

original_line="$(grep -E '^[[:space:]]*real\(wp\)[[:space:]]*::[[:space:]]*B_goal[[:space:]]*=' "$PARA_FILE")"
if [[ -z "$original_line" ]]; then
  echo "Could not find B_goal declaration in $PARA_FILE" >&2
  exit 1
fi

restore_original() {
  perl -0pi -e 's@^[ \t]*real\(wp\)[ \t]*::[ \t]*B_goal[ \t]*=.*$@'"$original_line"'@m' "$PARA_FILE"
}
trap restore_original EXIT

for ((b_value=START_VALUE; b_value<=END_VALUE; b_value+=STEP_VALUE)); do
  b_literal="${b_value}.e0_wp"
  run_tag="$(printf 'B_%05d' "$b_value")"
  run_log="${LOG_DIR}/${run_tag}.log"

  echo "=== Running B_goal = ${b_literal} ==="
  perl -0pi -e 's@^[ \t]*real\(wp\)[ \t]*::[ \t]*B_goal[ \t]*=.*$@  real(wp) :: B_goal   = '"$b_literal"'@m' "$PARA_FILE"

  make clean >"$run_log" 2>&1
  make MODE="$BUILD_MODE" -j8 >>"$run_log" 2>&1
  "$ROOT_DIR/build/bin/a.out" >>"$run_log" 2>&1
done

