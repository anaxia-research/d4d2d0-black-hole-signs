"""Rank-2 sign theorem: certified error terms.

c_{2,mu}(E) = hat c_mu(E) + m_mu(E)  (modular part + mock part), E = exponent of q.
 * hat c: exact Rademacher expansion of the weakly holomorphic vector-valued form hhat (lattice 2 kappa,
   S = eps2 (2kappa)^{-1/2} e(-mu nu/2kappa)); |hat c/Lhat - 1| <= Rhat(E), Rhat as in Part I with kappa -> 2 kappa.
 * mock part: |m_mu(E)| <= Mstar(E) = K_G * sum K_a K_b * sqrt(2 D1/E) * exp(4 pi sqrt(2 D1 E)), from
   |m(E)| <= e^{2 pi E y} G^+(y) h_a^+(y) h_b^+(y) at y = sqrt(2 D1/E) <= y0, with
   G^+(y) <= K_G/y^2 (H(N) <= N) and h^+_nu(y) <= K_nu y^{3/2} e^{2 pi D1/y} (S-transform of the rank-1 series).
 * integer invariants: extra |1/4 c_1(E/2)| <= Istar(E) = K/4 (2D1/E)^{3/4} exp(2 pi sqrt(2 D1 E)).
For E >= E* all three ratios are decreasing, so Rhat + (Mstar + Istar)/Lhat < 1 at E* gives positivity for all E >= E*.
"""
from fractions import Fraction as Fr
from flint import arb, ctx
from rig1 import polar, terms as terms1, R, q2a, PI, NU, K3
from vvmf import solve
from r2common import MODELS2, get
import sys
ctx.prec = 256

def terms2(L, P):
    D0 = max(D for nu, D, c in P); c00 = [c for nu, D, c in P if D == D0][0]
    T = []
    for nu, D, c in P:
        mult = 1 if nu == 0 or 2*nu == L else 2
        w = q2a(abs(c))/q2a(abs(c00))*(q2a(D)/q2a(D0))**(arb(5)/4)*mult
        s = (q2a(D)/q2a(D0)).sqrt()
        if D != D0: T.append((w, s))
        T.append((w*arb(L).sqrt(), s/2))
        T.append((w*arb(L).sqrt()*K3, s/3))
    return D0, c00, T

def Lhat(L, chi, D0, E):
    x = 4*PI*(q2a(D0)*E).sqrt()
    return 2*PI*chi/arb(L).sqrt()*(q2a(D0)/E)**(arb(5)/4)*x.bessel_i(NU)

def rank1_K(name, y0, M1=80):
    """K_nu with h^+_nu(y) <= K_nu y^{3/2} e^{2 pi D1/y} for 0 < y <= y0."""
    d, mus, off, h, xs, lab = solve(name, M1)
    k = d['kappa']; D1 = Fr(k + d['c2'], 24); chi_D = Fr(k, 6) + Fr(d['c2'], 12)
    P = polar(d, mus, off, h); D0r, c00r, Tr = terms1(d, P)
    assert D0r == D1
    Y0 = 1/q2a(y0); y0a = q2a(y0)
    A = {}
    for nu in mus:
        s = arb(0)
        for n in range(M1 + 1):
            c = h[nu][n]
            if c != 0: s += q2a(abs(c))*(-2*PI*q2a(off[nu] + n)*Y0).exp()
        Ef = q2a(off[nu] + M1 + 1)
        assert Ef >= 4*q2a(D1)/Y0**2 and Ef >= 1
        r, _ = R(Tr, D1, Ef); assert r <= 1
        B = 4*PI*q2a(chi_D)/arb(k).sqrt()*q2a(D1)**(arb(5)/4)
        tail = B*(4*PI*(q2a(D1)*Ef).sqrt() - 2*PI*Y0*Ef).exp()/(1 - (-PI*Y0).exp())
        A[nu] = (s + tail)*(-2*PI*q2a(D1)*Y0).exp()
    # all non-positive coefficients (exact, n <= M1; Part I: none beyond)
    neg = {nu: [(off[nu] + n, h[nu][n]) for n in range(M1 + 1) if h[nu][n] < 0] for nu in mus}
    for nu in mus:
        n = next(n for n in range(M1) if off[nu] + n > 0)
        while not R(Tr, D1, q2a(off[nu] + n))[0] < 1: n += 1
        assert n <= M1
    base = sum((A[nu]*(1 if nu == 0 or 2*nu == k else 2) for nu in mus), arb(0))/arb(k).sqrt()
    fac = y0a**(-arb(3)/2)*(-2*PI*q2a(D1)/y0a).exp()
    assert q2a(y0) < 4*PI*q2a(D1)/3
    K = {}
    for nu in mus:
        N = sum((q2a(abs(c))*max(arb(1), (-2*PI*q2a(E)*y0a).exp()) for E, c in neg[nu]), arb(0))
        K[nu] = base + 2*N*fac
    return dict(D1=D1, K=K, neg=neg, off=off, h=h, d=d)

def mock_terms(name, mu):
    """list of (lambda, a, b) with m_mu = +- sum G_lambda h_{1,a} h_{1,b}"""
    if name == 'X10': return [((mu + 1) % 2, 0, 0)]
    return [((mu - 2*mu1 + 2) % 4, mu1 % 2, (mu - mu1) % 2) for mu1 in range(2)]

def KG(name, y0):
    y0 = q2a(y0)
    if name == 'X10': g0, g2 = arb(1)/12, 4/PI**2
    else: g0, g2 = arb(1)/8, 9/PI**2
    return g0*y0**2 + g2

def setup(name, M=120, y0=Fr(1, 2)):
    d, out, h1, chi2D, chi_2D = get(name, M)
    k = d['kappa']; L = 2*k; mus = d['mus']
    off = {mu: out[mu][0][0] for mu in mus}; hh = {mu: out[mu][0][1] for mu in mus}
    P = polar(d, mus, off, hh)
    D0, c00, T = terms2(L, P)
    assert D0 == Fr(chi2D, 24)
    eps = (-1)**(int(chi_2D) + 1); assert Fr(int((eps*c00).p), int((eps*c00).q)) == chi_2D, (c00, chi_2D)
    r1 = rank1_K(name, y0)
    D1 = r1['D1']; Kg = KG(name, y0)
    delta = 4*PI*(q2a(D0).sqrt() - (2*q2a(D1)).sqrt())
    h1n = h1 if isinstance(h1, dict) else {0: h1}
    return dict(name=name, d=d, out=out, h1=h1n, L=L, mus=mus, P=P, D0=D0, chi=chi_2D, T=T, r1=r1, Kg=Kg, y0=y0,
                delta=delta, bmax=max(b for w, b in T), M=M)

def Mstar(S, mu, E):
    r1 = S['r1']; K = r1['K']; D1 = q2a(r1['D1'])
    s = sum((K[a]*K[b] for lam, a, b in mock_terms(S['name'], mu)), arb(0))
    return S['Kg']*s*(2*D1/E).sqrt()*(4*PI*(2*D1*E).sqrt()).exp()

def Istar(S, mu, E):
    ic = MODELS2[S['name']]['intcorr']
    if mu not in ic: return arb(0)
    r1 = S['r1']; D1 = q2a(r1['D1'])
    return r1['K'][ic[mu]]/4*(2*D1/E)**(arb(3)/4)*(2*PI*(2*D1*E).sqrt()).exp()

def Rtot(S, mu, E, integer=False):
    rh, x = R(S['T'], S['D0'], E)
    lh = Lhat(S['L'], q2a(S['chi']), S['D0'], E)
    m = Mstar(S, mu, E) + (Istar(S, mu, E) if integer else 0)
    return rh, m/lh, rh + m/lh, x

def Emin(S):
    """smallest E from which Rhat, Mstar/Lhat, Istar/Lhat (times sqrt(E+1)) are all decreasing and y <= y0"""
    D1 = q2a(S['r1']['D1'])
    return max(2*D1/q2a(S['y0'])**2, (7/S['delta'])**2, arb(1))

if __name__ == "__main__":
    for name in (sys.argv[1:] or ['X10', 'X8']):
        S = setup(name)
        print(name, "D0", S['D0'], "chi_2D", S['chi'], "D1", S['r1']['D1'], "delta", S['delta'].mid().str(5),
              "b", S['bmax'].mid().str(5), "K_G", S['Kg'].mid().str(5), "K", {k: v.mid().str(6) for k, v in S['r1']['K'].items()},
              "Emin", Emin(S).mid().str(5))
        for mu in S['mus']:
            o = S['out'][mu][0][0]
            for n in range(0, 400, 20):
                E = q2a(o + n)
                if E < Emin(S): continue
                rh, mr, tot, x = Rtot(S, mu, E, integer=True)
                print(f"  mu={mu} E={float(E.mid()):.2f} Rhat={float(rh.upper()):.3e} M/L={float(mr.upper()):.3e} tot={float(tot.upper()):.3e}")
                if tot < 1e-3: break
