"""Extended exact checks: integrality to n <= 400 and the third-order Turan inequality
4(a_n^2 - a_{n-1}a_{n+1})(a_{n+1}^2 - a_n a_{n+2}) > (a_n a_{n+1} - a_{n-1} a_{n+2})^2 for n* + 1 <= n <= 200,
for the nine quotient models (bicubic at t = 486), the six one-polar forms and the ten two-polar basis forms."""
from fractions import Fraction as Fr
from models3 import MODELS, derived
from forms3 import solve, offsets, polar_positions
M = 400
out = open('../logs/extended3.log', 'w')
def say(*a):
    print(*a); print(*a, file=out); out.flush()
jobs = []
for name, d in MODELS.items():
    if name in ('G7', 'Rey') or name.startswith('V'): continue
    pol = {0: [-18], 1: [486]} if name == 'B3' else None
    jobs.append((name + ' ' + d['label'], name, pol, True))
for (k, c2) in [(1, 10), (1, 22), (2, 8), (3, 6), (4, 4), (5, 2)]:
    jobs.append(('U(%d,%d)' % (k, c2), 'V%d_%d' % (k, c2), {mu: ([1] if mu == 0 else []) for mu in range(k//2 + 1)}, True))
for (k, c2) in [(1, 34), (1, 46), (2, 20), (3, 18), (11, 2)]:
    d = MODELS['V%d_%d' % (k, c2)]; off = offsets(d); P = polar_positions(d, off)
    for i, lab in ((0, 'G1'), (1, 'G2')):
        pol = {mu: [0]*sum(1 for (m, n) in P if m == mu) for mu in off}
        m, n = P[i]; pol[m][n] = 1
        jobs.append(('%s(%d,%d)' % (lab, k, c2), 'V%d_%d' % (k, c2), pol, False))
for lab, name, pol, turan in jobs:
    h, off = solve(name, M, pol)
    d = MODELS[name]; chi, chiD, N, w = derived(d); eps = (-1)**(int(chiD) + 1)
    integral = all(x.denominator == 1 for mu in h for x in h[mu])
    msg = "%-14s integral for n <= %d: %s" % (lab, M, integral)
    if turan:
        fails = []
        for mu in h:
            a = [eps*x for x in h[mu]] if lab.startswith('U') else h[mu]
            first = next(n for n in range(M) if off[mu] + n >= 0)
            for n in range(first + 1, 201):
                l = 4*(a[n]**2 - a[n-1]*a[n+1])*(a[n+1]**2 - a[n]*a[n+2]); r = (a[n]*a[n+1] - a[n-1]*a[n+2])**2
                if not l > r: fails.append((mu, n))
        msg += "; third-order Turan for n*+1 <= n <= 200: %s" % ("no failure" if not fails else fails)
    say(msg)
