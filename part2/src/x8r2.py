from rank2 import *
from fractions import Fraction as Fr
def pw(p, k, M):
    r = fmpq_poly([1])
    for _ in range(k): r = trunc(r*p, M)
    return r
def Gseries(M):
    """G_lambda = (1/2) V_2[H]_lambda, lambda in Z/4:
    V_2[H]_lambda = sum_{N = -lambda^2 mod 8} H(N) q^{N/8} + 2 [lambda even] H_{lambda/2}(2 tau)."""
    Hn = hurwitz(8*M + 16)
    G = {}
    for lam in range(4):
        r = (-lam*lam) % 8
        o = Fr(r, 8)
        p = [Fr(0)]*(M + 1)
        for n in range(M + 1):
            p[n] += Hn[8*n + r]
        if lam % 2 == 0:
            # H_{lam/2}(2 tau) = sum_{N = -(lam/2) mod 4} H(N) q^{N/2}
            s = (lam//2) % 2
            for N in range(0, 8*M + 16):
                if (N + s) % 4 == 0:
                    e = Fr(N, 2) - o
                    if e.denominator == 1 and 0 <= e <= M:
                        p[int(e)] += 2*Hn[N]
                    else:
                        assert e.denominator == 1 or True
        G[lam] = ser(o, [x/2 for x in p])
    return G
def build_x8(M):
    MM = M + 10
    e4, e6 = E4(MM), E6(MM)
    def mono(a, b): return trunc(pw(e4, a, MM)*pw(e6, b, MM), MM)
    A = -(mono(11,1)*27489243085 + mono(8,3)*189114116937 + mono(5,5)*177632766471 + mono(2,7)*21771164387)/69657034752
    B = (mono(12,0)*4132634222 + mono(9,2)*90332924791 + mono(6,4)*173271553815 + mono(3,6)*50628642229 + mono(0,8)*895654223)/5804752896
    C = -(mono(10,1)*4768083767 + mono(7,3)*25127775147 + mono(4,5)*17062071045 + mono(1,7)*1415010841)/241864704
    ei = eta_inv_pow(104, MM)
    d, mus, off, h, x, lab = solve('X8', MM)
    h1 = {0: (off[0], h[0]), 1: (off[1], h[1])}
    G = Gseries(MM)
    out = {}
    for mu in range(3):
        th = theta2(2, mu, MM); d1 = Dser(th, Fr(1, 2), MM); d2 = Dser(d1, Fr(5, 2), MM)
        terms = [(t[0] - Fr(104, 24), trunc(trunc(P*ei, MM)*t[1], MM)) for P, t in ((A, th), (B, d1), (C, d2))]
        hhat = sadd(sadd(terms[0], terms[1], MM), terms[2], MM)
        mock = None
        for mu1 in range(2):
            lam = (mu - 2*mu1 + 2) % 4
            a, b = mu1 % 2, (mu - mu1) % 2
            t = smul(G[lam], smul(h1[a], h1[b], MM), MM)
            mock = t if mock is None else sadd(mock, t, MM)
        mock = sscale(mock, (-1)**(mu + 2))
        out[mu] = (hhat, mock, sadd(hhat, mock, MM))
    return out, h1
if __name__ == "__main__":
    out, h1 = build_x8(12)
    for mu in range(3):
        hh, mk, tot = out[mu]
        o, p = tot
        k = 0
        while p[k] == 0 and k < 5: k += 1
        print(mu, "offset", o, [str(p[i]) for i in range(12)])
