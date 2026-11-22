"""Rank-2 sign theorem: exact coefficients below the analytic threshold, certified majorants above it.
Also checks the majorants against the exact coefficients (sanity)."""
from fractions import Fraction as Fr
from flint import arb, fmpq
from rigorous2 import setup, Rtot, Lhat, Mstar, Istar, q2a, PI, R
from r2common import MODELS2
import sys

def F(c): return Fr(int(c.p), int(c.q))

def intseries(S, mu, M):
    """coefficient list (exponent E, rational c2, integer cint) for the full h_{2,mu}"""
    out = S['out']; o, p = out[mu][2]
    ic = MODELS2[S['name']]['intcorr']; r1 = S['r1']
    res = []
    for n in range(M + 1):
        E = o + n; c = F(p[n]); ci = c
        if mu in ic:
            nu = ic[mu]; o1, p1 = S['h1'][nu]; e1 = E/2 - o1
            if e1.denominator == 1 and e1 >= 0:
                assert e1 <= p1.degree(), "rank-1 series too short"
                ci = c - F(p1[int(e1)])/4
        res.append((E, c, ci))
    return res

def run(name, out=sys.stdout):
    S = setup(name, M=160)
    D1 = S['r1']['D1']
    Es = max(2*q2a(D1)/q2a(S['y0'])**2, (6/S['delta'])**2, arb(1))
    print(f"{name}: L=2kappa={S['L']} D0=chi(2D)/24={S['D0']} chi_2D={S['chi']} D1={D1} delta={S['delta'].mid().str(6)} "
          f"b={S['bmax'].mid().str(6)} K_G={S['Kg'].mid().str(6)} K={[v.mid().str(6) for v in S['r1']['K'].values()]} "
          f"sign threshold E>={Es.mid().str(5)}", file=out)
    res = {}
    for mu in S['mus']:
        ser = intseries(S, mu, S['M'])
        first = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
        # analytic threshold
        i = first
        while True:
            E = ser[i][0]
            if q2a(E) >= Es:
                rh, mr, tot, x = Rtot(S, mu, q2a(E), integer=True)
                if tot < 1: break
            i += 1
        nstar = i
        negr = [(i - first, ser[i][0], ser[i][1]) for i in range(first, nstar) if not ser[i][1] > 0]
        negi = [(i - first, ser[i][0], ser[i][2]) for i in range(first, nstar) if not ser[i][2] > 0]
        integral = all(ser[i][2].denominator == 1 for i in range(len(ser)))
        # sanity: exact vs majorants for E in [5, 150]
        worst = [0, 0, 0]
        hat = S['out'][mu][0]; mock = S['out'][mu][1]
        for i in range(first, len(ser)):
            E = ser[i][0]
            if E < 3: continue
            ch = F(hat[1][int(E - hat[0])]); cm = F(mock[1][int(E - mock[0])])
            Ea = q2a(E); lh = Lhat(S['L'], q2a(S['chi']), S['D0'], Ea)
            rh, _ = R(S['T'], S['D0'], Ea)
            e1 = abs(q2a(ch)/lh - 1)/rh; e2 = abs(q2a(cm))/Mstar(S, mu, Ea) if Ea >= 2*q2a(D1)/q2a(S['y0'])**2 else arb(0)
            ici = abs(q2a(ser[i][2] - ser[i][1]))/Istar(S, mu, Ea) if (Ea >= 2*q2a(D1)/q2a(S['y0'])**2 and mu in MODELS2[name]['intcorr']) else arb(0)
            worst = [max(worst[0], float(e1.upper())), max(worst[1], float(e2.upper())), max(worst[2], float(ici.upper()))]
        res[mu] = dict(first=first, nstar=nstar - first, Estar=ser[nstar][0], negr=negr, negi=negi, tot=tot)
        print(f"   mu={mu}: Delta={-ser[0][0]} first non-polar E={ser[first][0]}; Rhat+(M+I)/Lhat < 1 certified at E={ser[nstar][0]} (n={nstar-first}, value {float(tot.upper()):.2e})", file=out)
        print(f"      rational invariants, non-positive non-polar: {[(n, str(E), str(c)) for n, E, c in negr] or 'none'}", file=out)
        print(f"      integer invariants,  non-positive non-polar: {[(n, str(E), str(c)) for n, E, c in negi] or 'none'}; integral up to n={len(ser)-1}: {integral}", file=out)
        print(f"      sanity (max over 3<=E<=~160): |hatc/Lhat-1|/Rhat <= {worst[0]:.3f}, |m|/Mstar <= {worst[1]:.3e}, |int corr|/Istar <= {worst[2]:.3e}", file=out)
    return S, res

if __name__ == "__main__":
    import os
    os.makedirs('../logs', exist_ok=True)
    with open('../logs/signs2.log', 'w') as f:
        for name in (sys.argv[1:] or ['X10', 'X8']):
            run(name, f); run(name)
