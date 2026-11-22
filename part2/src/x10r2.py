from rank2 import *
from fractions import Fraction as Fr
def pw(p, k, M):
    r = fmpq_poly([1])
    for _ in range(k): r = trunc(r*p, M)
    return r
def build_x10(M):
    e4, e6 = E4(M + 10), E6(M + 10)
    def mono(a, b): return trunc(pw(e4, a, M + 10)*pw(e6, b, M + 10), M + 10)
    A = (mono(12,0)*5397523 + mono(9,2)*70149738 - mono(6,4)*12112656 - mono(3,6)*61127530 - mono(0,8)*2307075)/46438023168
    B = (-mono(10,1)*10826123 - mono(7,3)*14574207 + mono(4,5)*20196255 + mono(1,7)*5204075)/1934917632
    ei = eta_inv_pow(100, M + 10)
    d, mus, off, h, x, lab = solve('X10', M + 10)
    h1 = (off[0], h[0])
    h1sq = smul(h1, h1, M + 10)
    H = Hseries(M + 10)
    out = {}
    for mu in (0, 1):
        th = theta2(1, mu, M + 10)
        dth = Dser(th, Fr(1, 2), M + 10)
        t1 = (th[0] - Fr(100, 24), trunc(trunc(A*ei, M + 10)*th[1], M + 10))
        t2 = (dth[0] - Fr(100, 24), trunc(trunc(B*ei, M + 10)*dth[1], M + 10))
        hhat = sadd(t1, t2, M + 10)
        mock = sscale(smul(H[(mu + 1) % 2], h1sq, M + 10), (-1)**(mu + 1))
        out[mu] = (hhat, mock, sadd(hhat, mock, M + 10))
    return out, h1
if __name__ == "__main__":
    out, h1 = build_x10(12)
    for mu in (0, 1):
        hh, mk, tot = out[mu]
        o, p = tot
        # strip leading zeros to the published offset
        k = 0
        while p[k] == 0: k += 1
        print(mu, "offset", o + k, [str(p[k + i]) for i in range(9)])
        print("   hhat offset", hh[0], "mock offset", mk[0])
