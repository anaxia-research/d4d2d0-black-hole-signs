"""Numerical check of the S-transformation
   h_mu(-1/tau) = (-i tau)^w sum_{nu in Z/K} B_{mu nu} h_nu(tau),  B_{mu nu} = eps K^{-1/2} e(-mu nu/K),
for the constructed forms (h_{-nu} = h_nu); determines the sign eps for each model, and checks T."""
from mpmath import mp, mpf, mpc, exp, pi, sqrt, cos
from fractions import Fraction as Fr
from models3 import MODELS, derived
from forms3 import solve, offsets
mp.dps = 30
def hval(h, off, mu, tau, M):
    q = exp(2j*pi*tau)
    o = off[mu]
    return exp(2j*pi*tau*mpf(o.numerator)/o.denominator)*sum(mpf(h[mu][n].numerator)/h[mu][n].denominator*q**n for n in range(M+1))
def comp(K, nu):
    nu %= K; return min(nu, K - nu)
if __name__ == "__main__":
    M = 60
    for name, d in MODELS.items():
        if name.startswith("V"): continue
        pol = None
        if name == 'B3': pol = {0: [-18], 1: [486]}
        if name in ('G7', 'Rey'): continue
        h, off = solve(name, M, pol)
        K = d['K']; chi, chiD, N, w = derived(d)
        res = []
        for tau in [mpc(0, 1.1), mpc(0.3, 0.9), mpc(-0.2, 1.3)]:
            for mu in h:
                lhs = hval(h, off, mu, -1/tau, M)
                rhs = (-1j*tau)**(mpf(w.numerator)/w.denominator)*sum(cos(2*pi*mu*nu/K)*hval(h, off, comp(K, nu), tau, M) for nu in range(K))/sqrt(K)
                res.append(lhs/rhs)
        eps = res[0]
        print(name, d['label'], "w =", w, " eps =", mp.nstr(eps, 15), " max deviation:", mp.nstr(max(abs(r - eps) for r in res), 3),
              " (-1)^(chi_D+1) =", (-1)**(int(chiD) + 1), " sign of eps*c_00:", "+" if (eps.real*h[0][0]) > 0 else "-")
