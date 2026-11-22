"""Direct check of Theorem B: |c_mu(n)/L(E) - 1| < R(E) in ball arithmetic, for all components of all
forms and all n <= M with E = off_mu + n > 0 (an independent test of the normalisation of the Rademacher
expansion: the sign eps, the factor K^{-1/2} and the Bessel index beta = 1 - w)."""
from fractions import Fraction as Fr
from flint import arb
from models3 import MODELS, derived
from forms3 import solve
from rig3 import Engine
M = 80
out = open('../logs/major3.log', 'w')
def say(*a): print(*a); print(*a, file=out); out.flush()
cases = []
for name, d in MODELS.items():
    if name in ('G7', 'Rey') or name.startswith('V'): continue
    pol = {0: [-18], 1: [d['t_mcg']]} if name == 'B3' else None
    h, off = solve(name, M, pol); chi, chiD, N, w = derived(d)
    cases.append((d['label'], d, h, off, (-1)**(int(chiD) + 1)))
for (k, c2) in [(1, 10), (1, 22), (2, 8), (3, 6), (4, 4), (5, 2)]:
    name = 'V%d_%d' % (k, c2); d = MODELS[name]; chi, chiD, N, w = derived(d); eps = (-1)**(int(chiD) + 1)
    h, off = solve(name, M, {mu: ([eps] if mu == 0 else []) for mu in range(k//2 + 1)})
    cases.append(('U(%d,%d)' % (k, c2), d, h, off, eps))
tot = 0
for lab, d, h, off, eps in cases:
    chi, chiD, N, w = derived(d)
    eng = Engine(d['K'], w, off, h, eps)
    cnt = 0; ok = True
    for mu in off:
        for n in range(M + 1):
            E = off[mu] + n
            if E <= 0: continue
            r, x = eng.R(E); L = eng.L(E)
            rel = abs(arb(h[mu][n].numerator)/h[mu][n].denominator/L - 1)
            if not (rel < r): ok = False; say("  FAIL", lab, mu, n)
            cnt += 1
    tot += cnt
    say("%-10s pairs checked: %3d  |c/L - 1| < R certified: %s" % (lab, cnt, ok))
say("total pairs:", tot)
