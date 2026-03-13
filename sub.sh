#!/bin/bash
set -e

make clean
make -j4 #> build.log 2>&1

./build/bin/a.out
