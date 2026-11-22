"""Numerical check of the S-transformation h_mu(-1/tau) = (-i tau)^{-3/2} sum_{nu mod kappa} B_{mu nu} h_nu(tau),
with the conjectured B_{mu nu} = eps * kappa^{-1/2} * cos(2 pi mu nu / kappa); determine eps for each model."""
from mpmath import mp, mpf, mpc, exp, pi, sqrt, cos
from vvmf import solve
from models import MODELS
mp.dps = 30
def hval(h, off, mu, tau, M):
    q = exp(2j*pi*tau)
    return exp(2j*pi*tau*mpf(off[mu].numerator)/off[mu].denominator)*sum(mpf(int(h[mu][n].p))/int(h[mu][n].q)*q**n for n in range(M+1))
def comp(kappa, nu):
    nu %= kappa; return min(nu, kappa - nu)
if __name__ == "__main__":
    M = 40
    for name in MODELS:
        d, mus, off, h, x, labels = solve(name, M)
        k = d['kappa']; res = []
        for tau in [mpc(0, 1.1), mpc(0.3, 0.9), mpc(-0.2, 1.3)]:
            for mu in mus:
                lhs = hval(h, off, mu, -1/tau, M)
                rhs = (-1j*tau)**(-1.5)*sum(cos(2*pi*mu*nu/k)*hval(h, off, comp(k, nu), tau, M) for nu in range(k))/sqrt(k)
                res.append(lhs/rhs)
        eps = res[0]
        print(name, "eps =", mp.nstr(eps, 12), " max deviation:", mp.nstr(max(abs(r - eps) for r in res), 3), " predicted (-1)^(chi_D+1):", (-1)**(int(k/6 + d['c2']/12) + 1))
