#!/bin/sh
# Reproduces every computation of Part III (a few minutes). Requires python3, python-flint (0.9), mpmath, numpy.
set -e
cd "$(dirname "$0")/src"
python3 closed3.py     > /dev/null              # Theorem A: closed forms, Wronskian free basis (writes ../logs/closed3.log)
python3 check3.py      > ../logs/check3.log     # Theorem A: dimensions, polar terms, agreement with the printed coefficients
python3 express3.py    > ../logs/express3.log   # generator expressions of the one-polar forms U_{kappa,c2}
python3 rc3.py         > ../logs/rc3.log        # Rankin-Cohen expressions of the numerators
python3 validate3.py   > /dev/null              # free-basis method validated on the hypergeometric models (writes ../logs/validate3.log)
python3 lattice3.py    > /dev/null              # lattices of the multiparameter models (writes ../logs/lattice3.log)
python3 smatrix3.py    > ../logs/smatrix3.log   # S-matrices, numerical check
python3 major3.py      > /dev/null              # Theorem B: |c/L - 1| < R certified (writes ../logs/major3.log)
python3 run_signs3.py  > /dev/null              # Theorems C and G: thresholds and exact checks (writes ../logs/signs3.log)
python3 classify3.py   > /dev/null              # Theorem D: invariants with at most two polar terms, sign windows (writes ../logs/classify3.log)
python3 window3.py     > /dev/null              # bicubic quotient: exact window for t (writes ../logs/window3.log)
python3 extended3.py   > /dev/null              # integrality to n <= 400, third-order Turan (writes ../logs/extended3.log)
python3 tables3.py     > ../logs/tables3.log    # tables and figure data (tab_*.txt, *.dat)
echo "all done"
