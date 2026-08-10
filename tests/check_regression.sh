#!/bin/bash
# Reference-value regression test.
# Compares Cont/properties.dat against tests/reference/*.ref
# with default relative tolerance 1e-4. A reference line may override this
# with a trailing `tol=<positive value>` annotation.
# Usage: bash tests/check_regression.sh <properties_file> <reference_file>
set -u

PROPS="${1:-Cont/properties.dat}"
REF="${2:-tests/reference/gr_uryu.ref}"
DEFAULT_TOL="1e-4"

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
COMPARED=0
OVERRIDES=0
while IFS= read -r ref_line || [ -n "$ref_line" ]; do
  case "$ref_line" in
    *"="*) ;;
    *) continue ;;
  esac
  # Extract the label (everything before =) and value (first number after =)
  label=$(echo "$ref_line" | sed 's/=.*//' | xargs)
  ref_val=$(echo "${ref_line#*=}" | grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' | head -1)
  if [ -z "$label" ] || [ -z "$ref_val" ]; then
    echo "FAIL: malformed numeric reference line: $ref_line"
    FAIL=1
    continue
  fi

  line_tol="$DEFAULT_TOL"
  if [[ "$ref_line" == *"tol="* ]]; then
    line_tol=$(echo "$ref_line" | sed -nE 's/.*[[:space:]]tol=([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?).*/\1/p')
    if [ -z "$line_tol" ] || ! awk -v tol="$line_tol" 'BEGIN { exit !(tol > 0) }'; then
      echo "FAIL: malformed tolerance annotation: $ref_line"
      FAIL=1
      continue
    fi
    OVERRIDES=$((OVERRIDES + 1))
  fi

  # Find matching line in properties file
  prop_line=$(grep -F "$label" "$PROPS" 2>/dev/null | head -1)
  if [ -z "$prop_line" ]; then
    echo "FAIL: missing reference label in properties: $label"
    FAIL=1
    continue
  fi
  prop_val=$(echo "${prop_line#*=}" | grep -oE '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?' | head -1)
  if [ -z "$prop_val" ]; then
    echo "FAIL: reference label has no numeric property value: $label"
    FAIL=1
    continue
  fi

  # Compute relative difference using awk
  if ! result=$(awk -v a="$ref_val" -v b="$prop_val" -v tol="$line_tol" '
    BEGIN {
      diff = (a == 0 && b == 0) ? 0 : (a == 0) ? b : (b - a) / a
      if (diff < 0) diff = -diff
      if (diff > tol) {
        printf "FAIL: %-20s ref=%-15s got=%-15s rdiff=%.2e tol=%.2e\n", "'"$label"'", a, b, diff, tol
        exit 1
      }
    }'); then
    echo "$result"
    FAIL=1
  fi
  COMPARED=$((COMPARED + 1))
done < "$REF"

if [ "$COMPARED" -eq 0 ]; then
  echo "FAIL: no numeric reference values were compared"
  FAIL=1
fi

if [ "$FAIL" -eq 0 ]; then
  echo "  regression: PASS ($COMPARED values; default tol=$DEFAULT_TOL, overrides=$OVERRIDES)"
  exit 0
else
  echo "  regression: FAIL"
  exit 1
fi
