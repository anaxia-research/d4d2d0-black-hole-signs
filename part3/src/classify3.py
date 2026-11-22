"""Classification of b2 = 1 invariants with at most two polar terms, the forms with one polar term,
and the sign windows for the forms with two polar terms."""
from fractions import Fraction as Fr
from models3 import MODELS, derived
from forms3 import solve, offsets, polar_positions, space, rref_rows
from rig3 import Engine
out = open('../logs/classify3.log', 'w')
def say(*a):
    print(*a); print(*a, file=out); out.flush()
M = 90
say("Invariants (kappa, c2) with chi_D = kappa/6 + c2/12 in Z_{>=1}, c2 >= 1, and at most two polar terms:")
cases = []
for k in range(1, 200):
    for c2 in range(1, 200):
        chiD = Fr(k, 6) + Fr(c2, 12)
        if chiD.denominator != 1 or chiD < 1: continue
        chi = k + c2
        if chi > 48: continue            # more than two polar terms in mu = 0 alone otherwise
        P = []
        for mu in range(k//2 + 1):
            off = Fr(-chi, 24) + (Fr(mu*mu, 2*k) + Fr(mu, 2)) % 1
            n = 0
            while off + n < 0: P.append((mu, str(off + n))); n += 1
        if len(P) <= 2: cases.append((k, c2)); say("  kappa=%d c2=%d chi(D)=%d chi_D=%d polar exponents %s" % (k, c2, chi, chiD, P))
# completeness of the search: chi(D) <= 48 is forced (mu = 0 has floor(chi/24)+1 polar terms), and kappa < chi(D)
say("search complete: chi(D) = kappa + c2 <= 48 is necessary, so kappa, c2 < 48")
for (k, c2) in cases:
    name = 'V%d_%d' % (k, c2); d = MODELS[name]
    chi, chiD, N, w = derived(d)
    off = offsets(d)
    P = polar_positions(d, off)
    mus, off, V = space(d, 6)
    img = rref_rows([[h[mu][n] for (mu, n) in P] for h in V])
    say("\n(%d,%d): N=%d, polar positions %s, dim V = %d" % (k, c2, N, P, len(img)))
    eps = (-1)**(int(chiD) + 1)
    if len(P) == 1:
        h, off = solve(name, M, {mu: ([1] if mu == 0 else []) for mu in off})
        he = {mu: [eps*x for x in h[mu]] for mu in h}
        eng = Engine(k, w, off, he, eps)
        res = eng.analyse(M)
        sgn = all(eps*h[mu][n] > 0 for mu in h for n in range(M + 1) if off[mu] + n >= 0)
        say("  U: ", {mu: [str(x) for x in h[mu][:5]] for mu in h})
        say("  eps = %d; eps*c(n) > 0 for all non-polar n <= %d: %s; integral: %s" % (eps, M, sgn, all(x.denominator == 1 for mu in h for x in h[mu])))
        for mu, r in res.items():
            say("   mu=%d first=%d n0=%d nD=%d nE=%d failures: pos %s growth %s lc %s" % (mu, r['first'], r['n0'], r['nD'], r['nE'], r['nonpos'], r['nongrowth'], r['nonlc']))
    else:
        (m1, n1), (m2, n2) = P
        def pol(a, b):
            pd = {mu: [0]*sum(1 for (m, n) in P if m == mu) for mu in offsets(d)}
            pd[m1][n1] = a; pd[m2][n2 if m2 != m1 else n2] = b
            return pd
        G1, off = solve(name, M, pol(1, 0)); G2, _ = solve(name, M, pol(0, 1))
        # eps * (G1 + t G2) > 0 on non-polar terms
        lo = hi = None; lo_at = hi_at = None
        for mu in G1:
            for n in range(M + 1):
                if off[mu] + n < 0: continue
                f1, f2 = eps*G1[mu][n], eps*G2[mu][n]
                if f2 == 0:
                    assert f1 > 0; continue
                b = -f1/f2
                if f2 > 0:
                    if lo is None or b > lo: lo, lo_at = b, (mu, n)
                else:
                    if hi is None or b < hi: hi, hi_at = b, (mu, n)
        say("  G1:", {mu: [str(x) for x in G1[mu][:5]] for mu in G1}, " G2:", {mu: [str(x) for x in G2[mu][:5]] for mu in G2})
        say("  window (from n <= %d): %s < t < %s   binding at %s, %s" % (M, lo, hi, lo_at, hi_at))
        if hi is None:
            # then eps*G2 >= 0 on all non-polar terms up to M; certify for all n
            he = {mu: [eps*x for x in G2[mu]] for mu in G2}
            o2 = {mu: off[mu] for mu in off}
            # G2 has its polar term at (m2, n2); as a form it is the one-polar form of a smaller model
            say("  eps*G2 non-polar coefficients all >= 0 for n <= %d:" % M, all(he[mu][n] >= 0 for mu in he for n in range(M + 1) if off[mu] + n >= 0))
        for t in (lo, hi):
            if t is None: continue
            h = {mu: [eps*(G1[mu][n] + t*G2[mu][n]) for n in range(M + 1)] for mu in G1}
            eng = Engine(k, w, off, h, eps)
            res = eng.analyse(M)
            say("  endpoint t=%s: zero/negative non-polar coefficients" % t, {mu: [(m, str(h[mu][m])) for m in r['nonpos']] for mu, r in res.items()},
                "R<1 from", {mu: r['n0'] for mu, r in res.items()})
