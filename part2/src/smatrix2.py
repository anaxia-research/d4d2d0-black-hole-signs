"""Numerical check of the S-matrix of the modular part hhat of the rank-2 series:
hhat_mu(-1/tau) = (-i tau)^{-3/2} sum_nu S_{mu nu} hhat_nu(tau), predicted
S_{mu nu} = eps2 (2 kappa)^{-1/2} e(-mu nu/(2 kappa)) summed over +-nu, eps2 = (-1)^{chi_{2D}+1}."""
import mpmath as mp
from r2common import get
mp.mp.dps = 40
def ev(S, tau):
    o, p = S
    q = mp.exp(2j*mp.pi*tau)
    return mp.exp(2j*mp.pi*tau*mp.mpf(o.numerator)/o.denominator)*sum(mp.mpf(int(p[i].p))/int(p[i].q)*q**i for i in range(p.degree() + 1))
def check(name, M=80):
    d, out, h1, chi2D, chi_2D = get(name, M)
    k = d['kappa']; mus = d['mus']; L = 2*k
    eps = (-1)**(int(chi_2D) + 1)
    worst = 0
    for tau in (mp.mpc(0.13, 1.1), mp.mpc(-0.21, 0.95), mp.mpc(0.05, 1.3)):
        lhs = [ev(out[mu][0], -1/tau) for mu in mus]
        rhs = []
        for mu in mus:
            s = 0
            for nu in mus:
                nus = {nu % L, (-nu) % L}
                s += sum(mp.exp(-2j*mp.pi*mu*n/L) for n in nus)*ev(out[nu][0], tau)
            rhs.append(eps/mp.sqrt(L)*(-1j*tau)**mp.mpf(-1.5)*s)
        for a, b in zip(lhs, rhs):
            worst = max(worst, abs(a - b)/abs(a))
    print(f"{name}: chi_2D={chi_2D} eps2={eps:+d}  max relative deviation {mp.nstr(worst, 5)}")
    # also check the mock part is NOT modular by itself (sanity) and the full h2 is not either
for n in ('X10', 'X8'): check(n)

def check1(name, M=60):
    """rank one: h_mu(-1/tau) = (-i tau)^{-3/2} eps1 kappa^{-1/2} sum_nu e(-mu nu/kappa) h_nu(tau)"""
    from vvmf import solve
    from fractions import Fraction as Fr
    d, mus, off, h, x, lab = solve(name, M)
    k = d['kappa']; chi_D = Fr(k, 6) + Fr(d['c2'], 12); eps = (-1)**(int(chi_D) + 1)
    worst = 0
    for tau in (mp.mpc(0.13, 1.1), mp.mpc(-0.21, 0.95), mp.mpc(0.05, 1.3)):
        for mu in mus:
            lhs = ev((off[mu], h[mu]), -1/tau)
            s = 0
            for nu in range(k):
                nn = min(nu, (-nu) % k)
                s += mp.exp(-2j*mp.pi*mu*nu/k)*ev((off[nn], h[nn]), tau)
            rhs = eps/mp.sqrt(k)*(-1j*tau)**mp.mpf(-1.5)*s
            worst = max(worst, abs(lhs - rhs)/abs(lhs))
    print(f"{name} rank one: chi_D={chi_D} eps1={eps:+d}  max relative deviation {mp.nstr(worst, 5)}")
for n in ('X10', 'X8'): check1(n)
