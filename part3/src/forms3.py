"""Exact construction of the generating series of the quotient models as vector-valued modular
forms h_mu = eta^{-N} sum_l P_l(E4,E6) D^l theta^{(K)}_mu, and the linear algebra of polar parts.

Conventions.  h_mu = q^{off_mu} sum_{n>=0} c_mu(n) q^n, where off_mu = -chi(D)/24 + {(chi(D)-N)/24 + (mu+K/2)^2/(2K)}
is the exponent allowed by the unitarity bound q0hat <= chi(D)/24.  A coefficient is polar if off_mu + n < 0."""
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq_mat, fmpq
from vvmf import E4, E6, theta, serre, mul, eta_inv_pow, trunc
from models3 import MODELS, derived

def F2q(x): return fmpq(x.numerator, x.denominator)

def offsets(d):
    chi, chiD, N, w = derived(d); K = d['K']
    return {mu: Fr(-chi, 24) + (Fr(chi - N, 24) + Fr((2*mu + K)**2, 8*K)) % 1 for mu in range(K//2 + 1)}

def basis(d, M):
    """All functions eta^{-N} E4^a E6^b D^l theta_mu of the right weight, as q-series in the
    normalisation of `offsets` (h_mu = q^{off_mu} * series); entries before the allowed range are kept
    and returned separately as 'beyond' coefficients (they must vanish).
    Returns (mus, off, cols, labels, pre) with cols[j][mu] a list of rationals, pre[mu] = number of
    leading coefficients beyond the unitarity bound."""
    chi, chiD, N, w = derived(d); K = d['K']
    mus = list(range(K//2 + 1))
    off = offsets(d)
    wnum = w + Fr(N, 2)                       # weight of the numerator
    Mx = M + 2*N                               # work to a larger order to accommodate shifts
    e4, e6 = E4(Mx), E6(Mx)
    th = {mu: [theta(K, mu, Mx)] for mu in mus}
    w0 = th[0][0][2]
    lmax = int((wnum - w0) // 2)
    for mu in mus:
        for l in range(lmax):
            o, p, ww = th[mu][-1]; th[mu].append(serre(o, p, ww, Mx))
    ei = eta_inv_pow(N, Mx)
    pw4 = [fmpq_poly([1])]; pw6 = [fmpq_poly([1])]
    for _ in range(int(wnum)//4 + 2): pw4.append(mul(pw4[-1], e4, Mx))
    for _ in range(int(wnum)//6 + 2): pw6.append(mul(pw6[-1], e6, Mx))
    cols, labels = [], []
    pre = {}
    for mu in mus:
        o = th[mu][0][0] - Fr(N, 24)          # exponent of the first stored coefficient
        s = off[mu] - o; assert s.denominator == 1
        pre[mu] = int(s)                       # >= 0 : number of coefficients beyond the bound
    for l in range(lmax + 1):
        rest = wnum - th[0][l][2]
        if rest < 0 or rest.denominator != 1: continue
        rest = int(rest)
        for b in range(rest//6 + 1):
            if (rest - 6*b) % 4: continue
            a = (rest - 6*b)//4
            mf = mul(pw4[a], pw6[b], Mx)
            col = {}
            for mu in mus:
                ser = mul(mul(mf, th[mu][l][1], Mx), ei, Mx)
                col[mu] = [Fr(int(ser[i].p), int(ser[i].q)) for i in range(M + pre[mu] + 1)]
            cols.append(col); labels.append((l, a, b))
    return mus, off, cols, labels, pre

def polar_positions(d, off):
    """list of (mu, n) of polar coefficients, most polar first within each mu."""
    P = []
    for mu, o in off.items():
        n = 0
        while o + n < 0: P.append((mu, n)); n += 1
    return P

def space(d, M):
    """Basis of the space V (forms with no term beyond the unitarity bound), as rational coefficient
    vectors h[mu][n], n = 0..M, together with the polar map."""
    mus, off, cols, labels, pre = basis(d, M)
    J = len(cols)
    rows = [[cols[j][mu][i] for j in range(J)] for mu in mus for i in range(pre[mu])]
    if rows:
        A = fmpq_mat(len(rows), J, [F2q(x) for r in rows for x in r])
        # kernel of A
        R, rank = A.rref()
        piv = []
        for i in range(rank):
            for j in range(J):
                if R[i, j] != 0: piv.append(j); break
        free = [j for j in range(J) if j not in piv]
        ker = []
        for f in free:
            v = [Fr(0)]*J; v[f] = Fr(1)
            for i, j in enumerate(piv):
                x = R[i, f]; v[j] = -Fr(int(x.p), int(x.q))
            ker.append(v)
    else:
        ker = [[Fr(int(i == j)) for i in range(J)] for j in range(J)]
    V = []
    for v in ker:
        h = {mu: [sum(v[j]*cols[j][mu][pre[mu] + n] for j in range(J)) for n in range(M + 1)] for mu in mus}
        V.append(h)
    return mus, off, V

def rref_rows(vecs):
    """row-reduced basis of the span of rational vectors"""
    if not vecs: return []
    A = fmpq_mat(len(vecs), len(vecs[0]), [F2q(x) for r in vecs for x in r])
    R, rank = A.rref()
    return [[Fr(int(R[i, j].p), int(R[i, j].q)) for j in range(len(vecs[0]))] for i in range(rank)]

def solve(name, M, polar=None):
    """the unique form in V with the given polar vector (dict mu -> list), or None if none exists"""
    d = MODELS[name]
    mus, off, V = space(d, M)
    P = polar_positions(d, off)
    pol = polar or d['polar']
    target = [Fr(pol[mu][n]) for (mu, n) in P]
    # solve sum x_i phi(V_i) = target
    cols = [[h[mu][n] for (mu, n) in P] for h in V]
    k = len(V)
    A = fmpq_mat(len(P), k + 1, [F2q(x) for i in range(len(P)) for x in [cols[j][i] for j in range(k)] + [target[i]]])
    R, rank = A.rref()
    piv = []
    for i in range(rank):
        for j in range(k + 1):
            if R[i, j] != 0: piv.append(j); break
    if k in piv: return None
    x = [Fr(0)]*k
    for i, j in enumerate(piv): x[j] = Fr(int(R[i, k].p), int(R[i, k].q))
    return {mu: [sum(x[j]*V[j][mu][n] for j in range(k)) for n in range(M + 1)] for mu in mus}, off

if __name__ == "__main__":
    import sys
    for name, d in MODELS.items():
        mus, off, V = space(d, 12)
        P = polar_positions(d, off)
        img = rref_rows([[h[mu][n] for (mu, n) in P] for h in V])
        print(name, d['label'], "offsets", {m: str(o) for m, o in off.items()}, "#polar", len(P), "dim V", len(V), "rank of polar map", len(img))
