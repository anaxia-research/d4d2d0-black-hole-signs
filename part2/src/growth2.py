"""Rank-2 growth c(n+1) > c(n) and log-concavity c(n)^2 > c(n-1) c(n+1), rational and integer invariants.
Analytic part with Rtot = Rhat + (Mstar + Istar)/Lhat (a majorant for both kinds of invariant):
  growth:   2 Rtot/(1 - Rtot) < a/(4 sqrt(E+1)) at E2, where (1-b)x >= 8, x >= 14, sqrt(E) >= 7/delta, y <= y0;
  log-conc: 32 Rtot(E')(E'+2)^2/sqrt(E') < a (1 - Rtot(E')) at E' = E3 - 1, (1-b)x >= 10, x >= 14, sqrt(E') >= 9/delta.
Exact comparison of the coefficients below the thresholds."""
from fractions import Fraction as Fr
from flint import arb
from rigorous2 import setup, Rtot, q2a, PI
from signs2 import intseries
import sys

def run(name, out=sys.stdout, M=420):
    S = setup(name, M=M)
    D0 = q2a(S['D0']); a = 4*PI*D0.sqrt(); b = S['bmax']; dl = S['delta']
    D1 = q2a(S['r1']['D1']); ymin = 2*D1/q2a(S['y0'])**2
    print(f"{name}: a=4pi sqrt(D0)={a.mid().str(6)}", file=out)
    res = {}
    for mu in S['mus']:
        ser = intseries(S, mu, S['M'])
        first = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
        i = first
        while True:
            E = q2a(ser[i][0]); x = a*E.sqrt()
            if E >= ymin and (1 - b)*x >= 8 and x >= 14 and E.sqrt() >= 7/dl:
                r = Rtot(S, mu, E, integer=True)[2]
                if r < 1 and 2*r/(1 - r) < a/(4*(E + 1).sqrt()): break
            i += 1
        n2 = i
        i = first + 1
        while True:
            Ep = q2a(ser[i - 1][0]); x = a*Ep.sqrt()
            if Ep >= ymin and (1 - b)*x >= 10 and x >= 14 and Ep.sqrt() >= 9/dl:
                r = Rtot(S, mu, Ep, integer=True)[2]
                if r < 1 and 32*r*(Ep + 2)**2/Ep.sqrt() < a*(1 - r): break
            i += 1
        n3 = i
        assert n3 + 1 < len(ser)
        R = {}
        for kind, j in (('rational', 1), ('integer', 2)):
            c = [ser[m][j] for m in range(len(ser))]
            mono = [m - first for m in range(first, n2 + 1) if not c[m + 1] > c[m]]
            lc = [m - first for m in range(first + 1, n3 + 1) if not c[m]**2 > c[m - 1]*c[m + 1]]
            R[kind] = (mono, lc)
            print(f"   mu={mu} {kind:8s}: growth analytic from n={n2-first} (E={ser[n2][0]}), exact failures of c(n+1)>c(n): {mono or 'none'};"
                  f" log-concavity analytic from n={n3-first} (E={ser[n3][0]}), exact failures: {lc or 'none'}", file=out)
        integral = all(ser[m][2].denominator == 1 for m in range(len(ser)))
        print(f"   mu={mu} integer invariants integral for n <= {len(ser)-1-first}: {integral}", file=out)
        res[mu] = dict(n2=n2 - first, n3=n3 - first, R=R, integral=integral, N=len(ser) - 1 - first)
    return S, res

if __name__ == "__main__":
    with open('../logs/growth2.log', 'w') as f:
        for name in (sys.argv[1:] or ['X10', 'X8']):
            run(name, f)
    print(open('../logs/growth2.log').read())
