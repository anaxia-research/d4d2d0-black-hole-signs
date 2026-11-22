"""Express the one-polar forms U_{kappa,c2} in the generators E4^a E6^b D^l theta / eta^N
(minimal-l solution), to obtain closed formulas."""
from fractions import Fraction as Fr
from flint import fmpq_mat, fmpq
from models3 import MODELS, derived
from forms3 import basis, polar_positions
def F2q(x): return fmpq(x.numerator, x.denominator)
for (k, c2) in [(1, 10), (1, 22), (2, 8), (3, 6), (4, 4), (5, 2)]:
    d = MODELS['V%d_%d' % (k, c2)]
    mus, off, cols, labels, pre = basis(d, 8)
    # restrict to generators with l <= d-1 (free basis)
    dd = k//2 + 1
    idx = [j for j, (l, a, b) in enumerate(labels) if l <= dd - 1]
    P = polar_positions(d, off)
    rows, rhs = [], []
    for mu in mus:
        for i in range(pre[mu]): rows.append([cols[j][mu][i] for j in idx]); rhs.append(Fr(0))
    for (mu, n) in P: rows.append([cols[j][mu][pre[mu] + n] for j in idx]); rhs.append(Fr(1) if (mu, n) == (0, 0) else Fr(0))
    A = fmpq_mat(len(rows), len(idx) + 1, [F2q(x) for r, v in zip(rows, rhs) for x in r + [v]])
    R, rank = A.rref()
    piv = []
    for i in range(rank):
        for j in range(len(idx) + 1):
            if R[i, j] != 0: piv.append(j); break
    assert len(idx) not in piv
    x = {}
    for i, j in enumerate(piv): x[labels[idx[j]]] = Fr(int(R[i, len(idx)].p), int(R[i, len(idx)].q))
    free = len(idx) - rank
    chi, chiD, N, w = derived(d)
    print("(%d,%d) N=%d generators (l<=%d): %d, rank %d, free %d:" % (k, c2, N, dd - 1, len(idx), rank, free),
          " + ".join("%s*E4^%d E6^%d D^%d theta" % (v, a, b, l) for (l, a, b), v in sorted(x.items()) if v != 0))
