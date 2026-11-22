"""Certified analysis (ball arithmetic) of the coefficients c_mu(n) of a weakly holomorphic
vector-valued modular form of weight w < 0 in the theta representation of level K, given its
exact coefficients.  nu = 1 - w is the Bessel index of the Rademacher expansion.

  * R(E): majorant of |c_mu(n)/L(E) - 1|, E = off_mu + n, where L(E) is the leading Rademacher term
          L(E) = 2 pi |c_00| K^{-1/2} (D0/E)^{nu/2} I_nu(4 pi sqrt(D0 E)).
          R(E) = sum_{(v,D) != (0,D0)} w_D I(s_D x)/I(x) + sqrt(K) sum_{(v,D)} w_D (I(s_D x/2) + K3 I(s_D x/3))/I(x).
  * n0: first n with R < 1 (so c_mu(n) > 0 for n >= n0);
  * nD: analytic threshold for c(n+1) > c(n):   x >= 4(nu+1), (1-b) x >= 2 nu + 3, 2R/(1-R) < a/(4 sqrt(E+1));
  * nE: analytic threshold for log-concavity:   x' >= x_nu (root of x^2-(6nu+7)x-4(nu+1)^2), (1-b) x' >= 2nu+5,
          32 R(E') (E'+2)^2/sqrt(E') < a (1 - R(E')), E' = E - 1.
Below the thresholds the exact rational coefficients decide."""
from fractions import Fraction as Fr
from flint import arb, ctx
import sys
ctx.prec = 256
PI = arb.pi()

def A(x): return arb(x.numerator)/x.denominator

class Engine:
    def __init__(self, K, w, off, h, eps):
        """eps = sign of the S-matrix, B_{mu nu} = eps K^{-1/2} e(-mu nu/K); the leading Rademacher term is
        eps*c_00 times a positive function, so we require eps*c_00 > 0."""
        self.K, self.off, self.h = K, off, h
        self.nu = A(1 - w)
        nu = self.nu
        self.K3 = arb(3)**nu*(nu.zeta() - 1 - arb(2)**(-nu))
        P = []
        for v in off:
            n = 0
            while off[v] + n < 0:
                if h[v][n] != 0: P.append((v, -(off[v] + n), h[v][n]))
                n += 1
        self.P = P
        D0 = max(D for v, D, c in P)
        top = [(v, c) for v, D, c in P if D == D0]
        assert len(top) == 1 and top[0][0] == 0, top
        self.D0, self.c00 = D0, top[0][1]
        assert eps*self.c00 > 0, 'leading term must be positive'
        T = []
        for v, D, c in P:
            mult = 1 if v == 0 or 2*v == K else 2
            wt = A(abs(c))/A(abs(self.c00))*(A(D)/A(D0))**(nu/2)*mult
            s = (A(D)/A(D0)).sqrt()
            if D != D0: T.append((wt, s))
            T.append((wt*arb(K).sqrt(), s/2))
            T.append((wt*arb(K).sqrt()*self.K3, s/3))
        self.T = T
        self.bmax = max(b for wt, b in T)
        self.a = 4*PI*A(D0).sqrt()
        nuf = float(nu.mid())
        self.xlc = (6*nuf + 7 + ((6*nuf + 7)**2 + 16*(nuf + 1)**2)**0.5)/2

    def R(self, E):
        x = self.a*A(E).sqrt()
        Ix = x.bessel_i(self.nu)
        return sum((wt*(b*x).bessel_i(self.nu)/Ix for wt, b in self.T), arb(0)), x

    def L(self, E):
        x = self.a*A(E).sqrt()
        return 2*PI*A(abs(self.c00))/arb(self.K).sqrt()*(A(self.D0)/A(E))**(self.nu/2)*x.bessel_i(self.nu)

    def analyse(self, M):
        nu = self.nu; nuf = float(nu.mid()); res = {}
        for mu in self.off:
            o = self.off[mu]; c = self.h[mu]
            first = next(n for n in range(M) if o + n >= 0)
            n = first
            while True:
                E = o + n
                if E > 0:
                    r, x = self.R(E)
                    if r < 1: break
                n += 1
            n0 = n
            nonpos = [m for m in range(first, n0) if not c[m] > 0]
            # growth
            n = first
            while True:
                E = o + n
                if E > 0:
                    r, x = self.R(E)
                    if x >= 4*(nuf + 1) and (1 - self.bmax)*x >= 2*nuf + 3 and r < 1 and 2*r/(1 - r) < self.a/(4*A(E + 1).sqrt()): break
                n += 1
            nD = n
            # log-concavity
            n = first + 1
            while True:
                Ep = o + n - 1
                if Ep > 0:
                    r, x = self.R(Ep)
                    if x >= self.xlc and (1 - self.bmax)*x >= 2*nuf + 5 and r < 1 and 32*r*A(Ep + 2)**2/A(Ep).sqrt() < self.a*(1 - r): break
                n += 1
            nE = n
            assert max(n0, nD + 1, nE + 1) <= M, "increase M"
            nongrowth = [m for m in range(first, nD) if not c[m + 1] > c[m]]
            nonlc = [m for m in range(first + 1, nE) if not c[m]**2 > c[m - 1]*c[m + 1]]
            res[mu] = dict(first=first, n0=n0, nD=nD, nE=nE, nonpos=nonpos, nongrowth=nongrowth, nonlc=nonlc)
        return res
