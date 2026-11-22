"""Certified check |c_mu(n)/L(E) - 1| < R(E) at every pair (mu, n) with E > 0 and n < 70, all models:
an independent test of the normalisation of the Rademacher expansion (sign eps, factor kappa^{-1/2})."""
from flint import arb, ctx
from vvmf import solve
from models import MODELS
from rigorous import polar, terms, R, q2a
from logconcave import L
ctx.prec = 300
bad = tot = 0
for name in MODELS:
    d, mus, off, h, xs, lab = solve(name, 82)
    P = polar(d, mus, off, h); D0, c00, T = terms(d, P); d['_c00'] = c00
    for mu in mus:
        for n in range(70):
            E = q2a(off[mu] + n)
            if not E > 0: continue
            l, _ = L(d, D0, E); r, _ = R(T, D0, E)
            ok = abs(arb(int(h[mu][n].p))/l - 1) < r
            tot += 1; bad += (not ok)
print(f"certified |c/L - 1| < R at {tot} pairs (n < 70, E > 0, all models); failures: {bad}")
