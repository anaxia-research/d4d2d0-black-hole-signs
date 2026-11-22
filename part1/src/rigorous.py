"""Rigorous (ball-arithmetic) versions of the sign and monotonicity theorems.

For each model and each mu:
  * R_mu(E): explicit majorant of |c_mu(n)/L(E) - 1|, E = Delta_mu + n, built from the Rademacher
    expansion (|Kl_c| <= phi(c) * mult_nu) and I_nu(x/c) <= (3/c)^nu I_nu(x/3) for c >= 3.
    Each term is w * I(b x)/I(x) with b < 1, x = 4 pi sqrt(D0 E), hence decreasing in E.
  * n0 = first index with R < 1 (certified in arb); then c_mu(n) > 0 for all n >= n0.
  * Exact rational coefficients for n < n0 give the complete list of non-positive values.
  * Monotonicity: E2 with (1 - b_max) x >= 8 and x >= 14; there R(E) sqrt(E+1) is decreasing and
    2 R/(1-R) < a/(4 sqrt(E+1)) is certified at E2; below E2 exact integers are compared.
"""
from flint import arb, ctx, fmpq
from vvmf import solve
from models import MODELS
import sys
ctx.prec = 256
PI = arb.pi()
NU = arb(5)/2
K3 = arb(3)**NU*(NU.zeta() - 1 - arb(2)**(-NU))

def q2a(x): return arb(int(x.p))/int(x.q) if hasattr(x, 'p') else arb(x.numerator)/x.denominator

def polar(d, mus, off, h):
    P = []
    for nu in mus:
        n = 0
        while off[nu] + n < 0:
            c = h[nu][n]
            if c != 0: P.append((nu, -(off[nu] + n), c))
            n += 1
    return P

def terms(d, P):
    """list of (weight, b) with R(E) = sum weight * I(b x)/I(x)"""
    k = d['kappa']
    D0 = max(D for nu, D, c in P); c00 = [c for nu, D, c in P if D == D0][0]
    T = []
    for nu, D, c in P:
        mult = 1 if nu == 0 or 2*nu == k else 2
        w = q2a(abs(c))/q2a(abs(c00))*(q2a(D)/q2a(D0))**(arb(5)/4)*mult
        s = (q2a(D)/q2a(D0)).sqrt()
        if D != D0: T.append((w, s))
        T.append((w*arb(k).sqrt(), s/2))
        T.append((w*arb(k).sqrt()*K3, s/3))
    return D0, c00, T

def R(T, D0, E):
    x = 4*PI*(q2a(D0)*E).sqrt()
    Ix = x.bessel_i(NU)
    return sum((w*(b*x).bessel_i(NU)/Ix for w, b in T), arb(0)), x

def run(name, M=48, out=sys.stdout):
    d, mus, off, h, xsol, labels = solve(name, M)
    P = polar(d, mus, off, h)
    D0, c00, T = terms(d, P)
    k = d['kappa']; chiD = k + d['c2']
    bmax = max(b for w, b in T)
    a = 4*PI*q2a(D0).sqrt()
    integral = all(h[mu][n].q == 1 for mu in mus for n in range(M + 1))
    print(f"{name}: kappa={k} chi(D)={chiD} Delta_0={D0} c_0(0)={c00} terms={len(T)} b_max={bmax.mid().str(6)} integral(n<={M})={integral}", file=out)
    res = {}
    for mu in mus:
        first = next(n for n in range(M) if off[mu] + n >= 0)
        n = first
        while True:
            E = q2a(off[mu] + n)
            if E > 0:
                r, x = R(T, D0, E)
                if r < 1: break
            n += 1
        n0 = n
        neg = [(m, int(h[mu][m].p)) for m in range(first, n0) if not h[mu][m] > 0]
        # monotonicity
        n = first
        while True:
            E = q2a(off[mu] + n)
            if E > 0:
                x = a*E.sqrt()
                if (1 - bmax)*x >= 8 and x >= 14:
                    r, _ = R(T, D0, E)
                    if r < 1 and 2*r/(1 - r) < a/(4*(E + 1).sqrt()): break
            n += 1
        n2 = n
        assert n2 + 1 <= M, "increase M"
        nonmono = [m for m in range(first, n2 + 1) if not h[mu][m + 1] > h[mu][m]]
        res[mu] = (first, n0, neg, n2, nonmono, r)
        print(f"   mu={mu}: Delta_mu={off[mu]}  first non-polar n={first}  R<1 certified from n={n0}"
              f"  non-positive coefficients: {neg if neg else 'none'}  monotone (analytic) from n={n2}"
              f"  exact failures of c(n+1)>c(n) on [{first},{n2}]: {nonmono if nonmono else 'none'}", file=out)
    return res

if __name__ == "__main__":
    for name in (sys.argv[1:] or MODELS):
        run(name)
