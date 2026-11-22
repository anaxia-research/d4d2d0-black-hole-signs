"""Validation of the free-basis computation of dim V on the hypergeometric models (kappa <= 7):
the number of obstructions should be 1 for X_{4,4} and 0 for the others."""
from models3 import MODELS, virtual
from forms3 import space, polar_positions, rref_rows, offsets
out = open('../logs/validate3.log', 'w')
for nm, k, c2 in [('X5', 5, 50), ('X6', 3, 42), ('X8', 2, 44), ('X10', 1, 34), ('X4,3', 6, 48), ('X6,4', 2, 32), ('X4,4', 4, 40), ('X6,6', 1, 22), ('X6,2', 4, 52)]:
    d = virtual(k, c2)
    mus, off, V = space(d, 4)
    P = polar_positions(d, off)
    img = rref_rows([[h[mu][n] for (mu, n) in P] for h in V])
    line = "%-5s kappa=%d c2=%d: polar terms %d, dim V %d, obstructions %d" % (nm, k, c2, len(P), len(img), len(P) - len(img))
    print(line); print(line, file=out)
