import sys
from fractions import Fraction as Fr
from models3 import MODELS, derived
from forms3 import solve
from rig3 import Engine
M = 80
out = open('../logs/signs3.log', 'w')
def say(*a):
    print(*a); print(*a, file=out); out.flush()
for name, d in MODELS.items():
    if name.startswith('V'): continue
    if name in ('G7', 'Rey'): continue
    pol = {0: [-18], 1: [d['t_mcg']]} if name == 'B3' else None
    h, off = solve(name, M, pol)
    chi, chiD, N, w = derived(d)
    eng = Engine(d['K'], w, off, h, (-1)**(int(chiD) + 1))
    integral = all(x.denominator == 1 for mu in h for x in h[mu])
    say(f"{name} {d['label']}: K={d['K']} w={w} nu={1-w} chi(D)={chi} D0={eng.D0} c00={eng.c00} b={eng.bmax.mid().str(5)} integral(n<={M})={integral}")
    res = eng.analyse(M)
    for mu, r in res.items():
        say(f"   mu={mu} off={off[mu]} first={r['first']} n0={r['n0']} nD={r['nD']} nE={r['nE']} nonpositive={r['nonpos']} growth failures={r['nongrowth']} log-concavity failures={r['nonlc']}")
