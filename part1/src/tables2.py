"""Threshold and exception tables from the certified runs."""
from fractions import Fraction as Fr
from vvmf import solve
from models import MODELS
from rigorous import run
import io, re
def fr(x):
    x = Fr(x)
    if x.denominator == 1: return str(x.numerator)
    s = '-' if x < 0 else ''
    return f"{s}\\tfrac{{{abs(x.numerator)}}}{{{x.denominator}}}"
def tex(name): return '$X_{' + name[1:] + '}$'
lc = {}
for line in open('../logs/logconcave.log'):
    pass
cur = None
for line in open('../logs/logconcave.log'):
    if not line.startswith(' '): cur = line.strip()
    else:
        m = re.search(r'mu=(\d+): analytic from n=(\d+)', line); lc[(cur, int(m.group(1)))] = int(m.group(2))
thr, exc = [], []
for name in MODELS:
    buf = io.StringIO(); res = run(name, out=buf)
    d, mus, off, h, x, l = solve(name, 40)
    from rigorous import polar, terms
    P = polar(d, mus, off, h); D0, c00, T = terms(d, P); b = max(bb for w, bb in T)
    cells = []
    for mu in mus:
        first, n0, neg, n2, nonmono, r = res[mu]
        cells.append(f"{mu}: {first}/{n0}/{n2}/{lc[(name, mu)]}")
        for (m, v) in neg:
            exc.append(f"{tex(name)} & {mu} & {m} & ${fr(-(off[mu] + m))}$ & ${v}$ \\\\")
    thr.append(f"{tex(name)} & {float(b.mid()):.4f} & " + "; ".join(cells) + r" \\")
open('../logs/tab_thr.txt', 'w').write("\n".join(thr))
open('../logs/tab_exc.txt', 'w').write("\n".join(exc))
print(len(exc)); print("\n".join(thr)); print("\n".join(exc))
