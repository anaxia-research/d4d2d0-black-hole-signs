"""The bicubic quotient P5[3,3]/Z3: h^t = -18 F1 + t F2, t = Omega_{1,1}(5/24) unknown.
Find the exact set of t for which every non-polar coefficient is positive, and certify it."""
from fractions import Fraction as Fr
from models3 import MODELS, derived
from forms3 import solve
from rig3 import Engine
M = 80
out = open('../logs/window3.log', 'w')
def say(*a):
    print(*a); print(*a, file=out); out.flush()
F1, off = solve('B3', M, {0: [1], 1: [0]})
F2, _ = solve('B3', M, {0: [0], 1: [1]})
say("F1:", {mu: [str(x) for x in F1[mu][:6]] for mu in F1})
say("F2:", {mu: [str(x) for x in F2[mu][:6]] for mu in F2})
lo, hi = None, None; lo_at = hi_at = None
for mu in F1:
    for n in range(M + 1):
        if off[mu] + n < 0: continue
        f1, f2 = F1[mu][n], F2[mu][n]
        if f2 == 0:
            assert -18*f1 > 0; continue
        b = Fr(18)*f1/f2
        if f2 > 0:
            if lo is None or b > lo: lo, lo_at = b, (mu, n)
        else:
            if hi is None or b < hi: hi, hi_at = b, (mu, n)
say("signs of F2 on non-polar terms:", {mu: ''.join('+' if F2[mu][n] > 0 else '-' for n in range(M + 1) if off[mu] + n >= 0)[:40] for mu in F2})
say("window from n<=%d: %s < t < %s, binding at %s and %s" % (M, lo, hi, lo_at, hi_at))
chi, chiD, N, w = derived(MODELS['B3'])
for t in (lo, hi):
    h = {mu: [-18*F1[mu][n] + t*F2[mu][n] for n in range(M + 1)] for mu in F1}
    eng = Engine(3, w, off, h, -1)
    res = eng.analyse(M)
    say("t =", t, {mu: (r['first'], r['n0'], r['nonpos'], [str(h[mu][m]) for m in r['nonpos']]) for mu, r in res.items()})
# integrality / divisibility
say("F1, F2 integral up to n=%d:" % M, all(x.denominator == 1 for G in (F1, F2) for mu in G for x in G[mu]))
say("F2 = 2 theta/eta^6 check: first coefficients", F2[1][:4])
