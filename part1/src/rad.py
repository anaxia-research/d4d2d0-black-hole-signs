"""Leading Rademacher term T1 versus exact coefficients, and the trivial bound for c >= 2."""
from mpmath import mp, mpf, pi, sqrt, cos, besseli, gamma, cosh, exp, zeta
from vvmf import solve
from models import MODELS
mp.dps = 40
def fq(x): return mpf(int(x.p))/int(x.q)
def comp(k, nu): nu %= k; return min(nu, k - nu)
def polar_list(d, mus, off, h):
    P = []
    for nu in mus:
        n = 0
        while off[nu] + n < 0:
            c = fq(h[nu][n])
            if c != 0: P.append((nu, -(mpf(off[nu].numerator)/off[nu].denominator + n), c))
            n += 1
    return P
def T1(d, mus, off, h, P, mu, n, eps):
    k = d['kappa']; E = mpf(off[mu].numerator)/off[mu].denominator + n
    tot = 0
    for nu in range(k):
        cc = cos(2*pi*mu*nu/k)
        for (nu2, D, c) in P:
            if nu2 == comp(k, nu):
                tot += cc*c*(D/E)**mpf(1.25)*besseli(2.5, 4*pi*sqrt(D*E))
    return 2*pi*eps*tot/sqrt(k)
def tailbound(d, P, mu, off, n):
    """2 pi sum_{c>=2} sum_{nu mod kappa} sum_m |c_nu(m)| (D/E)^{5/4} I_{5/2}(x/c), x = 4 pi sqrt(DE);
    I_{5/2}(y) <= (y/2)^{5/2} cosh(y)/Gamma(7/2) for c >= 3, exact for c = 2."""
    k = d['kappa']; E = mpf(off[mu].numerator)/off[mu].denominator + n
    tot = 0
    for (nu2, D, c) in P:
        mult = 1 if nu2 == 0 or 2*nu2 == k else 2
        x = 4*pi*sqrt(D*E)
        s = besseli(2.5, x/2) + (x/2)**2.5*cosh(x/3)/gamma(3.5)*(zeta(2.5) - 1 - mpf(2)**-2.5)
        tot += mult*abs(c)*(D/E)**1.25*s
    return 2*pi*tot
if __name__ == "__main__":
    import sys
    M = 60
    eps = {'X5':1,'X6':-1,'X8':-1,'X10':1,'X4,3':1,'X6,4':1,'X3,3':-1,'X4,4':-1,'X6,6':-1,'X6,2':1,'X4,2':-1}
    for name in (sys.argv[1:] or MODELS):
        d, mus, off, h, x, labels = solve(name, M)
        P = polar_list(d, mus, off, h)
        neg = [(mu, n, int(h[mu][n].p)) for mu in mus for n in range(M+1) if off[mu] + n > 0 and h[mu][n] < 0]
        print(name, "polar:", [(nu, mp.nstr(D, 5), mp.nstr(c, 8)) for nu, D, c in P])
        print("   negative non-polar coefficients (n<=%d):" % M, neg)
        for mu in mus:
            for n in [10, 30, 60]:
                ex = fq(h[mu][n]); t = T1(d, mus, off, h, P, mu, n, eps[name]); tb = tailbound(d, P, mu, off, n)
                print("   mu=%d n=%d exact/T1-1=%s  tail/T1=%s" % (mu, n, mp.nstr(ex/t - 1, 5), mp.nstr(tb/abs(t), 5)))
