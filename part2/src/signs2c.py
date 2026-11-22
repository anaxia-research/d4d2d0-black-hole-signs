"""Sharpest pointwise majorant: at each grid point E take y = y0 = (rational >= sqrt(2 D1/E)), with the
constants K(y0), K_G(y0) recomputed for that y0.  Used to certify R_tot(E) < 1 pointwise below the monotone range."""
from fractions import Fraction as Fr
from math import ceil, sqrt
from flint import arb
from rigorous2 import setup, rank1_K, KG, mock_terms, Lhat, R, q2a, PI
from r2common import MODELS2
from signs2 import intseries
import sys
_cache = {}
def Kfor(name, y0):
    if (name, y0) not in _cache: _cache[(name, y0)] = (rank1_K(name, y0)['K'], KG(name, y0))
    return _cache[(name, y0)]
def Rpoint(S, mu, E, integer=True):
    name = S['name']; D1 = S['r1']['D1']
    y0 = Fr(ceil(sqrt(2*float(D1)/float(E))*1000), 1000)
    if not (y0 <= 2 and y0 < 4*Fr(355, 113)*D1/3*Fr(99, 100)): return None
    K, Kg = Kfor(name, y0); y = q2a(y0); Ea = q2a(E); D1a = q2a(D1)
    m = Kg*sum((K[a]*K[b] for lam, a, b in mock_terms(name, mu)), arb(0))*y*(2*PI*Ea*y + 4*PI*D1a/y).exp()
    ic = MODELS2[name]['intcorr']
    if integer and mu in ic: m += K[ic[mu]]/4*y**(arb(3)/2)*(2*PI*(Ea/2)*y + 2*PI*D1a/y).exp()
    rh, x = R(S['T'], S['D0'], Ea)
    return rh + m/Lhat(S['L'], q2a(S['chi']), S['D0'], Ea)
def run(name, out=sys.stdout):
    S = setup(name, M=80)
    Emono = max(2*q2a(S['r1']['D1'])/q2a(S['y0'])**2, (6/S['delta'])**2)
    res = {}
    for mu in S['mus']:
        ser = intseries(S, mu, 80)
        first = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
        ok = []; i = first
        while True:
            E = ser[i][0]
            r = Rpoint(S, mu, E) if E > 0 else None
            ok.append(r is not None and r < 1)
            if q2a(E) >= Emono and ok[-1]: break
            i += 1
        n0 = max([k for k, v in enumerate(ok) if not v], default=-1) + 1
        r0 = Rpoint(S, mu, ser[first + n0][0])
        res[mu] = (n0, ser[first + n0][0], r0, i - first)
        print(f"{name} mu={mu}: R_tot<1 at every n >= n*+{n0} (E >= {ser[first+n0][0]}, value {float(r0.upper()):.3f}); pointwise to n*+{i-first}, monotone beyond", file=out)
    return res
if __name__ == "__main__":
    with open('../logs/signs2c.log', 'w') as f:
        for name in ('X10', 'X8'): run(name, f)
    print(open('../logs/signs2c.log').read())
