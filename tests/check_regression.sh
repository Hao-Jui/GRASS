#!/bin/bash
# Reference-value regression test.
# Compares Cont/properties.dat against tests/reference/*.ref
# with relative tolerance 1e-4 (4 significant digits).
# Usage: bash tests/check_regression.sh <properties_file> <reference_file>
set -e

PROPS="${1:-Cont/properties.dat}"
REF="${2:-tests/reference/gr_uryu.ref}"
TOL="1e-4"

if [ ! -f "$PROPS" ]; then
  echo "FAIL: $PROPS not found"
  exit 1
fi
if [ ! -f "$REF" ]; then
  echo "FAIL: $REF not found"
  exit 1
fi

# Extract numbers from both files and compare
FAIL=0
while IFS= read -r ref_line; do
  # Extract the label (everything before =) and value (first number after =)
  label=$(echo "$ref_line" | sed 's/=.*//' | xargs)
  ref_val=$(echo "$ref_line" | grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' | head -1)
  [ -z "$ref_val" ] && continue

  # Find matching line in properties file
  prop_line=$(grep "$label" "$PROPS" 2>/dev/null | head -1)
  [ -z "$prop_line" ] && continue
  prop_val=$(echo "$prop_line" | grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' | head -1)
  [ -z "$prop_val" ] && continue

  # Compute relative difference using awk
  result=$(awk -v a="$ref_val" -v b="$prop_val" -v tol="$TOL" '
    BEGIN {
      diff = (a == 0 && b == 0) ? 0 : (a == 0) ? b : (b - a) / a
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "FAIL: %-20s ref=%-15s got=%-15s rdiff=%.2e\n", "'"$label"'", a, b, diff
        exit 1
      }
    }')
  if [ $? -ne 0 ]; then
    echo "$result"
    FAIL=1
  fi
done < "$REF"

if [ $FAIL -eq 0 ]; then
  echo "  regression: PASS (all values within $TOL relative tolerance)"
  exit 0
else
  echo "  regression: FAIL"
  exit 1
fi
