"""Data tables (model data, polar terms, first non-polar coefficients) in LaTeX tabular row format."""
from fractions import Fraction as Fr
from flint import arb, ctx
from vvmf import solve, build
from models import MODELS
from rigorous import polar, terms, R, q2a, PI
import math
ctx.prec = 256
def tex(name): return name.replace('X', 'X_{').replace(',', ',') + '}'
def fr(x):
    x = Fr(x)
    if x.denominator == 1: return str(x.numerator)
    s = '-' if x < 0 else ''
    return f"{s}\\tfrac{{{abs(x.numerator)}}}{{{x.denominator}}}"
rows1, rows2, rows3, app, dens = [], [], [], [], []
for name in MODELS:
    d, mus, off, h, xsol, labels = solve(name, 40)
    k, c2 = d['kappa'], d['c2']; chiD = k + c2; chi_D = Fr(k, 6) + Fr(c2, 12)
    P = polar(d, mus, off, h); D0, c00, T = terms(d, P)
    bmax = max(b for w, b in T)
    Ds = sorted({D for nu, D, c in P}, reverse=True)
    D1 = Ds[1] if len(Ds) > 1 else None
    npol = sum(1 for nu in mus for n in range(10) if off[nu] + n < 0)
    rows1.append(f"${tex(name)}$ & {k} & {c2} & {d['chiY']} & {int(chi_D)} & {chiD} & {len(labels)} & {npol} & " + ", ".join(f"${fr(off[mu])}$" for mu in mus) + r" \\")
    pol = "; ".join(f"$\\mu{{=}}{mu}$: " + ",\\,".join(str(v) for v in d['polar'][mu]) for mu in mus)
    rows2.append(f"${tex(name)}$ & {pol}" + r" \\")
    # thresholds
    first = {mu: next(n for n in range(40) if off[mu] + n >= 0) for mu in mus}
    rows3.append((name, fr(D0), fr(D1) if D1 else '--', float(bmax.mid())))
    for mu in mus:
        cs = [int(h[mu][n].p) for n in range(first[mu], first[mu] + 4)]
        app.append(f"${tex(name)}$ & {mu} & ${fr(off[mu] + first[mu])}$ & " + " & ".join(f"${c}$" for c in cs[:3]) + r" \\")
    dx = 1
    for v in xsol: dx = dx*int(v.q)//math.gcd(dx, int(v.q))
    dens.append((name, dx))
open('../logs/tab_models.txt', 'w').write("\n".join(rows1))
open('../logs/tab_polar.txt', 'w').write("\n".join(rows2))
open('../logs/tab_first.txt', 'w').write("\n".join(app))
print(rows3); print(dens)
