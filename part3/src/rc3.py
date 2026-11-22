"""Rankin-Cohen expressions of the numerators of the one-polar forms:
U_{kappa,c2} * eta^N = alpha [E4,theta]_2 + beta [E6,theta]_1 + gamma E4^2 theta   (kappa = 4, 5)."""
from fractions import Fraction as Fr
from math import comb
from flint import fmpq_poly, fmpq, fmpq_mat
from vvmf import E4, E6, theta, mul, eta_inv_pow, trunc
from forms3 import solve
M = 30
def qd(o, p):
    return fmpq_poly([p[i]*(fmpq(o.numerator, o.denominator) + i) for i in range(p.degree() + 1)])
def binom(x, k):  # x rational
    r = Fr(1)
    for i in range(k): r *= (x - i)
    for i in range(1, k + 1): r /= i
    return r
def RC(f, kf, g, og, kg, n):
    """[f, g]_n with f = scalar series (offset 0), g = q^og * series"""
    out = fmpq_poly([])
    fd = [f]; gd = [g]
    for _ in range(n): fd.append(qd(Fr(0), fd[-1])); gd.append(qd(og, gd[-1]))
    for r in range(n + 1):
        s = n - r
        c = (-1)**r*binom(Fr(n + kf - 1), s)*binom(n + kg - 1, r)
        out += mul(fd[r], gd[s], M)*fmpq(c.numerator, c.denominator)
    return trunc(out, M)
e4, e6 = E4(M), E6(M)
for (k, c2) in [(2, 8), (3, 6), (4, 4), (5, 2)]:
    N = 4*k + c2
    h, off = solve('V%d_%d' % (k, c2), 12, {mu: ([1] if mu == 0 else []) for mu in range(k//2 + 1)})
    cand = {}
    for mu in range(k//2 + 1):
        o, th, w = theta(k, mu, M)
        if k <= 3:
            cand.setdefault('[E4,th]_1', {})[mu] = RC(e4, 4, th, o, w, 1)
            cand.setdefault('E6 th', {})[mu] = mul(e6, th, M)
        else:
            cand.setdefault('[E4,th]_2', {})[mu] = RC(e4, 4, th, o, w, 2)
            cand.setdefault('[E6,th]_1', {})[mu] = RC(e6, 6, th, o, w, 1)
            cand.setdefault('E4^2 th', {})[mu] = mul(mul(e4, e4, M), th, M)
        cand.setdefault('_o', {})[mu] = o
    names = [n for n in cand if n != '_o']
    # target: h * eta^N  -> compare coefficient lists of numerator: numerator = q^{o} series
    rows, rhs = [], []
    for mu in range(k//2 + 1):
        o = cand['_o'][mu]
        num = mul(eta_inv_pow(-N, M), fmpq_poly([0]), M)
        # build h_mu * eta^N as series in q^{o}: h_mu = q^{off}*sum c(n) q^n, eta^N = q^{N/24} prod
        s = off[mu] + Fr(N, 24) - o; assert s.denominator == 1 and s >= 0; s = int(s)
        hs = fmpq_poly([0]*s + [fmpq(x.numerator, x.denominator) for x in h[mu]])
        tgt = mul(hs, eta_inv_pow(-N, M), 12)
        for i in range(12):
            rows.append([Fr(int(cand[n][mu][i].p), int(cand[n][mu][i].q)) for n in names]); rhs.append(Fr(int(tgt[i].p), int(tgt[i].q)))
    A = fmpq_mat(len(rows), len(names) + 1, [fmpq(x.numerator, x.denominator) for r, v in zip(rows, rhs) for x in r + [v]])
    R, rank = A.rref()
    sol = {}
    for i in range(rank):
        j = next(j for j in range(len(names) + 1) if R[i, j] != 0)
        assert j < len(names), "not in span"
        sol[names[j]] = Fr(int(R[i, len(names)].p), int(R[i, len(names)].q))
    print("(%d,%d): eta^%d U =" % (k, c2, N), " + ".join("(%s)*%s" % (v, n) for n, v in sol.items()))
