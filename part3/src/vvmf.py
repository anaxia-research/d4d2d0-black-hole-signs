"""Exact construction of the Abelian D4-D2-D0 generating series h_mu of a one-parameter
hypergeometric model as the unique weight -3/2 vector-valued modular form with given polar part.
Each h_mu = q^{-Delta_mu} * sum_n c_mu(n) q^n, with Delta_mu in Q, n >= 0."""
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq_mat, fmpq, fmpz_poly
MODELS = {}

def sigma(n, k): return sum(d**k for d in range(1, n+1) if n % d == 0)
def eis(k, c, M):
    return fmpq_poly([1] + [c*sigma(n, k-1) for n in range(1, M+1)])
def E2(M): return eis(2, -24, M)
def E4(M): return eis(4, 240, M)
def E6(M): return eis(6, -504, M)
def trunc(p, M): return fmpq_poly([p[i] for i in range(min(p.degree()+1, M+1))]) if p.degree() >= 0 else p
def mul(a, b, M): return trunc(a*b, M)
def eta_inv_pow(N, M):
    """prod (1-q^n)^{-N} to order M"""
    # log-derivative recursion: a_n = (N/n) sum_{k=1}^n sigma(k) a_{n-k}
    a = [Fr(1)] + [Fr(0)]*M
    s = [0] + [sigma(k, 1) for k in range(1, M+1)]
    for n in range(1, M+1):
        a[n] = Fr(N, n)*sum(s[k]*a[n-k] for k in range(1, n+1))
    return fmpq_poly([int(x) for x in a])

def theta(kappa, mu, M):
    """theta_mu^{(kappa)} = q^{o} * poly ; returns (o, poly, weight)"""
    terms = {}
    jr = range(-2*M-5, 2*M+6)
    for j in jr:
        a2 = Fr(2*(kappa*j + mu) + kappa, 2)        # kappa*k with k = j + mu/kappa + 1/2
        e = a2*a2/(2*kappa)
        if kappa % 2 == 0:
            c = Fr((-1)**((mu + kappa//2) % 2))
        else:
            # -i kappa (-1)^{kappa k} k = -(-1)^{(kappa+1)/2} (-1)^{kappa j + mu} kappa k
            c = Fr(-(-1)**(((kappa+1)//2) % 2) * (-1)**((kappa*j + mu) % 2)) * a2
        terms[e] = terms.get(e, 0) + c
    o = min(terms)
    poly = [Fr(0)]*(M+1)
    for e, c in terms.items():
        d = e - o; assert d.denominator == 1
        if d <= M: poly[int(d)] += c
    w = Fr(1, 2) if kappa % 2 == 0 else Fr(3, 2)
    return o, fmpq_poly([fmpq(x.numerator, x.denominator) for x in poly]), w

def serre(o, p, w, M):
    """D_w f = q d/dq f - w/12 E2 f, f = q^o p"""
    dp = fmpq_poly([p[i]*(fmpq(o.numerator, o.denominator) + i) for i in range(p.degree()+1)])
    return o, trunc(dp - mul(p, E2(M), M)*fmpq(w.numerator, w.denominator)/12, M), w + 2

def build(name, M):
    d = MODELS[name]; kappa, c2 = d['kappa'], d['c2']; N = 4*kappa + c2
    mus = list(range(0, kappa//2 + 1))
    wF = Fr(N, 2) - Fr(3, 2)
    e4, e6 = E4(M), E6(M)
    # D^l theta for each mu
    Dth = {mu: [theta(kappa, mu, M)] for mu in mus}
    lmax = int((wF - Dth[0][0][2]) // 2)
    for mu in mus:
        for l in range(lmax):
            o, p, w = Dth[mu][-1]; Dth[mu].append(serre(o, p, w, M))
    cols, labels = [], []
    pw4 = [fmpq_poly([1])]; pw6 = [fmpq_poly([1])]
    for _ in range(int(wF)//4 + 2): pw4.append(mul(pw4[-1], e4, M))
    for _ in range(int(wF)//6 + 2): pw6.append(mul(pw6[-1], e6, M))
    for l in range(lmax + 1):
        rest = wF - Dth[0][l][2]
        if rest < 0 or rest.denominator != 1: continue
        rest = int(rest)
        for b in range(rest//6 + 1):
            if (rest - 6*b) % 4: continue
            a = (rest - 6*b)//4
            mf = mul(pw4[a], pw6[b], M)
            cols.append({mu: (Dth[mu][l][0], mul(mf, Dth[mu][l][1], M)) for mu in mus}); labels.append((l, a, b))
    ei = eta_inv_pow(N, M)
    off = {mu: Dth[mu][0][0] - Fr(N, 24) for mu in mus}
    # h_mu columns: q^{off} * (poly * ei); re-index so that h_mu = q^{-chi/24 + frac(mu^2/2kappa + mu/2)} * sum_n c(n) q^n
    shift = {}
    for mu in mus:
        target = Fr(-(kappa + c2), 24) + (Fr(mu*mu, 2*kappa) + Fr(mu, 2)) % 1
        s = target - off[mu]; assert s.denominator == 1 and s >= 0
        shift[mu] = int(s); off[mu] = target
    H = [{mu: mul(c[mu][1], ei, M) for mu in mus} for c in cols]
    return d, mus, off, H, labels, shift

def shiftdown(p, s, M):
    for i in range(s): assert p[i] == 0, "term beyond the unitarity bound"
    return fmpq_poly([p[i] for i in range(s, M+1)])

def solve(name, M=12):
    d, mus, off, H, labels, shift = build(name, M)
    rows, rhs = [], []
    for mu in mus:
        for n, v in enumerate([0]*shift[mu] + d['polar'][mu]):
            rows.append([H[j][mu][n] for j in range(len(H))]); rhs.append(v)
    A = fmpq_mat(len(rows), len(H) + 1, [x for r, v in zip(rows, rhs) for x in r + [fmpq(v)]])
    R, rank = A.rref()
    piv = []
    for i in range(rank):
        for j in range(len(H) + 1):
            if R[i, j] != 0: piv.append(j); break
    assert len(H) not in piv, "inconsistent polar data"
    x = [fmpq(0)]*len(H)
    for i, j in enumerate(piv): x[j] = R[i, len(H)]
    h = {mu: shiftdown(sum((H[j][mu]*x[j] for j in range(len(H)) if x[j] != 0), fmpq_poly([])), shift[mu], M) for mu in mus}
    return d, mus, off, h, x, labels

if __name__ == "__main__":
    import sys
    for name in (sys.argv[1:] or MODELS):
        d, mus, off, h, x, labels = solve(name, 14)
        ok = True
        for mu in mus:
            got = [h[mu][n] for n in range(len(d['polar'][mu]) + len(d['check'][mu]))]
            want = d['polar'][mu] + d['check'][mu]
            if [int(g.p) if g.q == 1 else g for g in got] != want:
                ok = False
                print(name, mu, "offset", off[mu], "\n  got ", [str(g) for g in got], "\n  want", want)
        print(name, "offsets", {mu: str(off[mu]) for mu in mus}, "matches AFKPS data:", ok)
