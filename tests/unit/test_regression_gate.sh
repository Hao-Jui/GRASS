#!/bin/bash
set -u

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
props=$(mktemp /tmp/grass-regression-props.XXXXXX) || exit 1
ref=$(mktemp /tmp/grass-regression-ref.XXXXXX) || exit 1
trap 'rm -f "$props" "$ref"' EXIT

passed=0
failed=0

expect_status() {
  expected=$1
  name=$2
  if bash "$root/tests/check_regression.sh" "$props" "$ref" >/dev/null 2>&1; then
    actual=0
  else
    actual=1
  fi
  if [ "$actual" -eq "$expected" ]; then
    passed=$((passed + 1))
  else
    echo "FAIL: $name expected status $expected, got $actual"
    failed=$((failed + 1))
  fi
}

printf 'M4/M^5 = -1.00015\n' > "$props"
printf 'M4/M^5 = -1.0 tol=2e-4\n' > "$ref"
expect_status 0 "field override accepts value inside boundary"

printf 'M4/M^5 = -1.00021\n' > "$props"
expect_status 1 "field override rejects value outside boundary"

printf 'ADM Mass = 1.00015\n' > "$props"
printf 'ADM Mass = 1.0\n' > "$ref"
expect_status 1 "default tolerance remains strict"

printf 'M4/M^5 = -1.0 tol=0\n' > "$ref"
expect_status 1 "non-positive override is refused"

echo "regression_gate: $passed passed, $failed failed"
test "$failed" -eq 0
