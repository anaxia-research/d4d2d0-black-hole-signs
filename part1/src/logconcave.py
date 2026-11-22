"""Log-concavity c(n)^2 > c(n-1) c(n+1) for n >= first + 1, and the closed-form leading term.

Analytic part (E - 1 >= E3, x = 4 pi sqrt(Delta_0 t) >= 14 on [E-1, E+1]):
   (log L)''(t) = (6 - x + x^2 psi_xx - x psi_x)/(4 t^2) <= -(x - 7)/(4 t^2) <= -x/(8 t^2),
so the second difference of log L is <= -a sqrt(E-1)/(8 (E+1)^2), while the error contributes at most
4 R(E-1)/(1 - R(E-1)).  G(E') = R(E') (E'+2)^2/sqrt(E') is decreasing once (1 - b_max) x >= 10, so it suffices
to certify 32 R(E') (E'+2)^2 / sqrt(E') < a (1 - R(E')) at E' = E3.
Exact check for first + 1 <= n below the analytic threshold."""
from flint import arb, ctx
from vvmf import solve
from models import MODELS
from rigorous import polar, terms, R, q2a, PI, NU
import sys
ctx.prec = 256

def L(d, D0, E):
    k = d['kappa']; chiD = q2a(D0)*24
    chi_D = abs(int(d['_c00']))
    x = 4*PI*(q2a(D0)*E).sqrt()
    return 2*PI*chi_D/arb(k).sqrt()*(q2a(D0)/E)**(arb(5)/4)*x.bessel_i(NU), x

def closed(d, D0, E):
    k = d['kappa']; chi_D = abs(int(d['_c00']))
    x = 4*PI*(q2a(D0)*E).sqrt()
    F = 1 - 3/x + 3/x**2 - (-2*x).exp()*(1 + 3/x + 3/x**2)
    return chi_D*q2a(D0)/arb(2*k).sqrt()*E**(-arb(3)/2)*x.exp()*F

if __name__ == "__main__":
    M = 220
    out = open('../logs/logconcave.log', 'w')
    for name in (sys.argv[1:] or MODELS):
        d, mus, off, h, xsol, labels = solve(name, M + 12)
        P = polar(d, mus, off, h); D0, c00, T = terms(d, P); d['_c00'] = c00
        bmax = max(b for w, b in T); a = 4*PI*q2a(D0).sqrt()
        line = [name]
        for mu in mus:
            first = next(n for n in range(M) if off[mu] + n >= 0)
            # analytic threshold for log-concavity at n (uses E' = E - 1)
            n = first + 1
            while True:
                Ep = q2a(off[mu] + n - 1)
                if Ep > 0:
                    x = a*Ep.sqrt()
                    if x >= 14 and (1 - bmax)*x >= 10:
                        r, _ = R(T, D0, Ep)
                        if r < 1 and 32*r*(Ep + 2)**2/Ep.sqrt() < a*(1 - r): break
                n += 1
            n3 = n
            assert n3 + 1 <= M
            c = [int(h[mu][m].p) for m in range(M + 1)]
            fails = [m for m in range(first + 1, n3 + 1) if not c[m]**2 > c[m-1]*c[m+1]]
            # closed form = L check at n = 20
            E = q2a(off[mu] + 20); l, _ = L(d, D0, E); cf = closed(d, D0, E)
            line.append(f"mu={mu}: analytic from n={n3}, exact failures on [{first+1},{n3}]: {fails if fails else 'none'}; |closed/L-1| at n=20 < {float(abs(cf/l - 1).upper()):.1e}; c/L-1 at n=20 = {float(((arb(c[20])/l) - 1).mid()):.3e}")
        print("\n   ".join(line)); print("\n   ".join(line), file=out)
