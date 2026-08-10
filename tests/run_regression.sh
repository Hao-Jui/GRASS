#!/bin/bash
# Run one integration model and compare only output created by that run.
set -eu

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <executable> <properties_file> <reference_file>"
  exit 2
fi

EXE="$1"
PROPS="$2"
REF="$3"
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

# Some Fortran STOP paths return status zero. Removing stale output ensures
# those paths cannot pass by comparing a previous successful model.
rm -f -- "$PROPS"
"$EXE"
bash "$SCRIPT_DIR/check_regression.sh" "$PROPS" "$REF"
