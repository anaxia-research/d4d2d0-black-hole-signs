"""Closed forms and the Wronskian criterion.
 (1) h^{(1,21)} = 25 E4/eta^11, h^{(2,29)} = h^{(2,20)} = 9 E4/eta^12;
 (2) kappa = 2: W = -(1/2)[E4, theta]_1 = -2 E4 theta' + (1/4) E4' theta  (' = q d/dq), h^{(1,8)} = 49 W/eta^16,
     h^{(1,6)} = 25 W/eta^16, h^{(2,12)} = 25 W/eta^17, h^{(3,18)} = h^{(3,15)} = 9 W/eta^18;
 (3) bicubic: F2 = 2 theta/eta^6, F1 = (E6^2 theta/24 - E6 [E4,theta]_1/4)/eta^30 - 54 theta/eta^6;
 (4) Wronskian det(D^j theta_mu)_{0<=j<d, 0<=mu<=K/2} = c * eta^{24 Lambda} for K = 1..7."""
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq, fmpq_mat
from vvmf import E4, E6, theta, serre, mul, eta_inv_pow, trunc
from forms3 import solve
M = 60
def qd(o, p):  # q d/dq of q^o p  -> (o, p')
    return fmpq_poly([p[i]*(fmpq(o.numerator, o.denominator) + i) for i in range(p.degree() + 1)])
def coeffs(p, n): return [Fr(int(p[i].p), int(p[i].q)) for i in range(n)]
log = open('../logs/closed3.log', 'w')
def say(*a): print(*a); print(*a, file=log)
e4, e6 = E4(M + 40), E6(M + 40)
e4p = qd(Fr(0), e4)
def bracket(o, th, wth):   # [E4, theta]_1 = 4 E4 theta' - wth E4' theta
    return trunc(4*mul(e4, qd(o, th), M + 40) - fmpq(wth.numerator, wth.denominator)*mul(e4p, th, M + 40), M + 40)
# (1)
for name, c, Nn in (('Q5', 25, 11), ('M229', 9, 12), ('M220', 9, 12)):
    h, off = solve(name, M)
    g = mul(e4, eta_inv_pow(Nn, M + 40), M + 40)
    say(name, "= %d E4/eta^%d :" % (c, Nn), coeffs(g*c, M + 1) == h[0][:M + 1])
# (2)
for name, c, Nn in (('R7', 49, 16), ('H5', 25, 16), ('M212', 25, 17), ('M318', 9, 18), ('M315', 9, 18)):
    h, off = solve(name, M)
    ok = True
    for mu in (0, 1):
        o, th, wth = theta(2, mu, M + 40)
        W = trunc(-bracket(o, th, wth)/2, M + 40)
        g = mul(W, eta_inv_pow(Nn, M + 40), M + 40)*c
        # align: h_mu = q^{off} sum c(n) q^n ; g = q^{o - Nn/24} * series
        s = off[mu] - (o - Fr(Nn, 24)); assert s.denominator == 1 and s >= 0
        s = int(s)
        gc = coeffs(g, M + 1 + s)
        ok &= all(x == 0 for x in gc[:s]) and gc[s:s + M + 1] == h[mu][:M + 1]
    say(name, "= %d W/eta^%d with W = -(1/2)[E4,theta]_1:" % (c, Nn), ok)
# W has integer coefficients
okW = True
for mu in (0, 1):
    o, th, wth = theta(2, mu, M + 40)
    W = trunc(-bracket(o, th, wth)/2, M + 40)
    okW &= all(x.denominator == 1 for x in coeffs(W, M + 40))
say("W integral (n < %d):" % (M + 40), okW)
# (3)
F1, off = solve('B3', M, {0: [1], 1: [0]})
F2, _ = solve('B3', M, {0: [0], 1: [1]})
ok2 = ok1 = True
def aligned(G, Nn, o, mu, n):
    """coefficient of q^{off_mu + n} in G = q^{o - Nn/24} * series"""
    s = off[mu] + n - (o - Fr(Nn, 24)); assert s.denominator == 1
    s = int(s)
    return Fr(int(G[s].p), int(G[s].q)) if s >= 0 else Fr(0)
for mu in (0, 1):
    o, th, wth = theta(3, mu, M + 40)
    g2 = mul(th, eta_inv_pow(6, M + 40), M + 40)*2
    br = bracket(o, th, wth)
    num = trunc(mul(mul(e6, e6, M + 40), th, M + 40)/24 - mul(e6, br, M + 40)/4, M + 40)
    g1 = mul(num, eta_inv_pow(30, M + 40), M + 40)
    ok2 &= all(aligned(g2, 6, o, mu, n) == F2[mu][n] for n in range(M + 1))
    ok1 &= all(aligned(g1, 30, o, mu, n) - 27*aligned(g2, 6, o, mu, n) == F1[mu][n] for n in range(M + 1))
    # g1 must have no terms below the unitarity bound
    s0 = int(off[mu] - (o - Fr(30, 24)))
    ok1 &= all(g1[i] == 0 for i in range(s0))
say("F2 = 2 theta/eta^6:", ok2, "  F1 = (E6^2 theta/24 - E6[E4,theta]_1/4)/eta^30 - 54 theta/eta^6:", ok1)
# (4) Wronskian
for K in range(1, 8):
    d = K//2 + 1
    rows = []
    th = [[theta(K, mu, 40)] for mu in range(d)]
    for mu in range(d):
        for j in range(d - 1):
            o, p, w = th[mu][-1]; th[mu].append(serre(o, p, w, 40))
    # leading exponent of each row and the determinant of the matrix of q-series: compute det as a series
    # det of d x d matrix of power series (after removing q^{o_mu} from row mu)
    o_mu = [th[mu][0][0] for mu in range(d)]
    Lam = sum(o_mu)
    from itertools import permutations
    def sign(p):
        s = 1; p = list(p)
        for i in range(len(p)):
            while p[i] != i:
                j = p[i]; p[i], p[j] = p[j], p[i]; s = -s
        return s
    det = fmpq_poly([])
    for perm in permutations(range(d)):
        term = fmpq_poly([sign(perm)])
        for mu in range(d):
            term = mul(term, th[mu][perm[mu]][1], 40)
        det += term
    lead = next(i for i in range(41) if det[i] != 0)
    K_w = sum(th[0][j][2] for j in range(d))
    # compare with c eta^{24 Lambda}: eta^{24L} = q^L prod (1-q^n)^{24 L}
    c = det[lead]
    eta24 = eta_inv_pow(-int(24*Lam), 40) if (24*Lam).denominator == 1 else None
    same = eta24 is not None and lead == 0 and all(det[i] == c*eta24[i] for i in range(40))
    say("K=%d d=%d Lambda=%s weights sum=%s 12*Lambda=%s leading coefficient %s at q^%s; det = c eta^{24 Lambda}: %s" % (K, d, Lam, K_w, 12*Lam, c, Lam + lead, same))
