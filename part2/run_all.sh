#!/bin/sh
# Reproduces every computation of Part II (a few minutes). Requires python3, python-flint (0.9), mpmath.
set -e
cd "$(dirname "$0")/src"
python3 x10r2.py     > ../logs/x10r2.log      # exact rank-two series for X10 (mock part + modular part)
python3 x8r2.py      > ../logs/x8r2.log       # exact rank-two series for X8
python3 published2.py > ../logs/published2.log # agreement with the 88 published coefficients
python3 smatrix2.py  > ../logs/smatrix2.log   # S-matrices of the modular parts (ranks one and two)
python3 rigorous2.py > ../logs/rigorous2.log  # constants: rank-one majorant K, class-number constant K_G, merger gap
python3 signs2.py    > /dev/null              # Theorem C: uniform threshold (writes ../logs/signs2.log)
python3 signs2c.py   > /dev/null              # Theorem C: pointwise thresholds (writes ../logs/signs2c.log)
python3 growth2.py   > /dev/null              # Theorems D and E: growth and log-concavity (writes ../logs/growth2.log)
python3 turan2.py    > /dev/null              # higher-order Turan inequalities (writes ../logs/turan2.log)
python3 tables2.py   > ../logs/tables2.log    # tables and figure data
echo "all done"
