"""Rank-2 (two units of D4 charge) generating series of X10 and X8 (Alexandrov-Feyzbakhsh-Klemm-Pioline,
CNTP 20 (2026), eqs. (4.3), (4.11)): h2 = hhat + mock part, hhat in M_2(Y), mock part built from Hurwitz
class numbers and the rank-1 series h1.  Series are stored as (offset, fmpq_poly) with h = q^offset * poly(q)."""
from fractions import Fraction as Fr
from flint import fmpq_poly, fmpq
from vvmf import E2, E4, E6, eta_inv_pow, mul, trunc, solve, serre
import math

def hurwitz(Nmax):
    """H(N) for 0 <= N <= Nmax (H(0) = -1/12, H(N) = 0 unless N = 0, 3 mod 4)."""
    H = [Fr(0)]*(Nmax + 1); H[0] = Fr(-1, 12)
    a = 1
    while 3*a*a <= Nmax:
        for b in range(-a + 1, a + 1):
            # c >= a, N = 4ac - b^2
            c = a
            while True:
                N = 4*a*c - b*b
                if N > Nmax: break
                if c >= a and not (c == a and b < 0):
                    w = Fr(1)
                    if a == c and b == 0: w = Fr(1, 2)
                    elif a == b == c: w = Fr(1, 3)
                    H[N] += w
                c += 1
        a += 1
    return H

def ser(off, coeffs):
    return (Fr(off), fmpq_poly([fmpq(x.numerator, x.denominator) if isinstance(x, Fr) else x for x in coeffs]))

def smul(A, B, M):
    return (A[0] + B[0], trunc(A[1]*B[1], M))

def sadd(A, B, M):
    """add two series whose offsets differ by an integer"""
    d = A[0] - B[0]; assert d.denominator == 1
    if d < 0: A, B = B, A; d = -d
    d = int(d)
    sh = fmpq_poly([0]*d + [A[1][i] for i in range(A[1].degree() + 1)]) if A[1].degree() >= 0 else fmpq_poly([])
    return (B[0], trunc(sh + B[1], M))

def sscale(A, c): return (A[0], A[1]*c)

def Hseries(M):
    """H_mu(tau) = sum_{N = -mu mod 4} H(N) q^{N/4}, mu in Z/2; returned with integer-spaced exponents."""
    H = hurwitz(4*M + 8)
    H0 = ser(0, [H[4*n] for n in range(M + 1)])
    H1 = ser(Fr(3, 4), [H[4*n + 3] for n in range(M + 1)])
    return {0: H0, 1: H1}

def theta2(kappa, mu, M):
    """tilde theta_mu^{(kappa)} = sum_{k in Z + mu/(2kappa)} q^{kappa k^2}, as (offset, poly), weight 1/2"""
    terms = {}
    for j in range(-3*M - 5, 3*M + 6):
        e = Fr((2*kappa*j + mu)**2, 4*kappa)
        terms[e] = terms.get(e, 0) + 1
    o = min(terms)
    p = [0]*(M + 1)
    for e, c in terms.items():
        d = e - o; assert d.denominator == 1
        if d <= M: p[int(d)] += c
    return (o, fmpq_poly(p))

def Dser(A, w, M):
    o, p = serre(A[0], A[1], w, M)[:2]
    return (o, p)
