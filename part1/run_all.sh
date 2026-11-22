#!/bin/sh
# Reproduces every computation of Part I (about one minute). Requires python3, python-flint (0.9), mpmath.
set -e
cd "$(dirname "$0")/src"
python3 vvmf.py          > ../logs/construction.log   # Theorem A: exact forms, agreement with the 161 published coefficients
python3 smatrix.py       > ../logs/smatrix.log        # S-matrix check (Proposition 5.1)
python3 rad.py           > ../logs/rademacher_c1.log  # leading Rademacher term against exact coefficients
python3 rigorous.py      > ../logs/rigorous.log       # Theorems C and D: certified thresholds, exact checks below them
python3 logconcave.py    > /dev/null                  # Theorem E: writes ../logs/logconcave.log
python3 certify_error.py > ../logs/certify_error.log  # |c/L-1| < R at 2097 pairs
python3 extended.py 200  > ../logs/extended_200.log   # integrality, positivity, growth, log-concavity for n <= 200
python3 turan.py         > ../logs/turan.log          # higher-order Turan inequalities for n <= 200
python3 tables.py        > ../logs/tables.log
python3 tables2.py       > ../logs/tables2.log
python3 examples.py      > ../logs/examples.log
python3 figdata.py       > ../logs/figdata.log
echo "all done"
