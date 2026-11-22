"""Worked-example tables: exact c_0(n), leading term L(E), relative error and majorant, for X_{6,6}, X_{10}, X_5;
and the entropy constants of Theorem F for all models."""
from flint import arb, ctx
from vvmf import solve
from models import MODELS
from rigorous import polar, terms, R, q2a, PI
from logconcave import L
from fractions import Fraction as Fr
import math
ctx.prec = 256
def fr(x):
    x = Fr(x)
    if x.denominator == 1: return str(x.numerator)
    return f"{'-' if x<0 else ''}\\tfrac{{{abs(x.numerator)}}}{{{x.denominator}}}"
def sci(v):
    v = float(v)
    if v == 0: return "0"
    e = math.floor(math.log10(abs(v))); m = v/10**e
    return f"{m:.3f}\\cdot10^{{{e}}}" if e != 0 else f"{m:.3f}"
rows = []
for name, nmax in [('X6,6', 8), ('X10', 8), ('X5', 8)]:
    d, mus, off, h, xs, lab = solve(name, 20)
    P = polar(d, mus, off, h); D0, c00, T = terms(d, P); d['_c00'] = c00
    first = next(n for n in range(20) if off[0] + n > 0)
    for n in range(first, first + nmax):
        E = q2a(off[0] + n); l, x = L(d, D0, E); r, _ = R(T, D0, E)
        c = int(h[0][n].p)
        rows.append(f"$X_{{{name[1:]}}}$ & {n} & ${fr(off[0]+n)}$ & ${c}$ & ${sci(l.mid())}$ & ${sci((arb(c)/l - 1).mid())}$ & ${sci(r.mid())}$ \\\\")
open('../logs/tab_examples.txt', 'w').write("\n".join(rows))
rows = []
for name in MODELS:
    d, mus, off, h, xs, lab = solve(name, 20)
    P = polar(d, mus, off, h); D0, c00, T = terms(d, P)
    k = d['kappa']; chiD = k + d['c2']; chi_D = abs(int(c00))
    C = arb(chi_D*chiD)/(24*arb(2*k).sqrt())
    b = max(bb for w, bb in T)
    rate = (1 - b)*2*PI*(arb(chiD)/6).sqrt()
    c01 = int(h[0][1].p) if off[0] + 1 < 0 else None
    pred = (-1)**(chi_D + 1)*d['chiY']*(chi_D - 1)
    last = f"${c01}$ & ${pred}$" if c01 is not None else "-- & --"
    rows.append(f"$X_{{{name[1:]}}}$ & {chiD} & {float(C.mid()):.4f} & {float(b.mid()):.4f} & {float(rate.mid()):.4f} & {last} \\\\")
open('../logs/tab_const.txt', 'w').write("\n".join(rows))
print(open('../logs/tab_examples.txt').read()); print(open('../logs/tab_const.txt').read())
