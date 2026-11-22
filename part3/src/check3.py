from fractions import Fraction as Fr
from models3 import MODELS
from forms3 import space, polar_positions, rref_rows, solve
M = 12
tot = 0
for name, d in MODELS.items():
    if name.startswith('V'): continue
    if name in ('G7', 'Rey', 'B3'): continue
    h, off = solve(name, M)
    ok = True
    for mu in h:
        P = [n for n in range(M) if off[mu] + n < 0]
        start = len(P)
        got = h[mu][:start + len(d['check'][mu])]
        want = d['polar'][mu] + d['check'][mu]
        if got != [Fr(x) for x in want]: ok = False; print(name, mu, got, want)
        tot += len(d['check'][mu])
    print(name, "reproduces McGovern:", ok, {mu: [str(x) for x in h[mu][:6]] for mu in h})
print("printed non-polar coefficients checked:", tot)
# bicubic: two-parameter family
for t in [486]:
    h, off = solve('B3', M, polar={0: [-18], 1: [t]})
    d = MODELS['B3']
    print("B3 t=%d" % t, {mu: [str(x) for x in h[mu][:6]] for mu in h})
    for mu in h: print("  match", h[mu][1:1+len(d['check'][mu])] == [Fr(x) for x in d['check'][mu]])
# obstruction functionals
for name in ('G7', 'Rey'):
    d = MODELS[name]
    mus, off, V = space(d, M)
    P = polar_positions(d, off)
    img = rref_rows([[h[mu][n] for (mu, n) in P] for h in V])
    # annihilator: vector a with a . v = 0 for all v in img
    from flint import fmpq_mat, fmpq
    A = fmpq_mat(len(img), len(P), [fmpq(x.numerator, x.denominator) for r in img for x in r])
    # solve nullspace
    R, rank = A.rref()
    piv = []
    for i in range(rank):
        for j in range(len(P)):
            if R[i, j] != 0: piv.append(j); break
    free = [j for j in range(len(P)) if j not in piv]
    for f in free:
        a = [Fr(0)]*len(P); a[f] = Fr(1)
        for i, j in enumerate(piv): a[j] = -Fr(int(R[i, f].p), int(R[i, f].q))
        print(name, "polar positions", [(mu, str(off[mu] + n)) for (mu, n) in P])
        print("   image (rref):", [[str(x) for x in r] for r in img])
        print("   obstruction functional:", [str(x) for x in a])
        v = [d['polar'][mu][n] for (mu, n) in P]
        print("   value on McGovern (3.55):", sum(x*y for x, y in zip(a, v)))
