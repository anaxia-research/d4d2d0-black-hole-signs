"""Data files for the figures: relative error |c/L - 1| against the certified majorant R(E)."""
from flint import arb, ctx
from vvmf import solve
from models import MODELS
from rigorous import polar, terms, R, q2a
from logconcave import L
import math
ctx.prec = 256
M = 70
for name, tag in [('X5', 'x5'), ('X10', 'x10'), ('X6,2', 'x62')]:
    d, mus, off, h, xs, lab = solve(name, M + 12)
    P = polar(d, mus, off, h); D0, c00, T = terms(d, P); d['_c00'] = c00
    first = next(n for n in range(M) if off[0] + n > 0)
    with open(f'../logs/err_{tag}.dat', 'w') as f:
        f.write("E err maj\n")
        for n in range(first, M + 1):
            E = q2a(off[0] + n); l, x = L(d, D0, E); r, _ = R(T, D0, E)
            err = abs(arb(int(h[0][n].p))/l - 1)
            f.write(f"{float(E.mid()):.4f} {math.log10(float(err.mid())):.4f} {math.log10(float(r.mid())):.4f}\n")
with open('../logs/maj_all.dat', 'w') as f:
    names = list(MODELS)
    f.write("E " + " ".join(n.replace(',', '') for n in names) + "\n")
    data = {}
    for name in names:
        d, mus, off, h, xs, lab = solve(name, 20)
        P = polar(d, mus, off, h); D0, c00, T = terms(d, P)
        data[name] = (D0, T)
    for i in range(1, 61):
        E = arb(i)/4
        f.write(f"{i/4:.2f} " + " ".join(f"{math.log10(float(R(data[nm][1], data[nm][0], E)[0].mid())):.4f}" for nm in names) + "\n")
print("ok")
