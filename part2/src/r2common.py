"""Common interface for the two rank-2 models."""
from fractions import Fraction as Fr
from x10r2 import build_x10
from x8r2 import build_x8
MODELS2 = {'X10': dict(kappa=1, c2=34, build=build_x10, mus=[0, 1], intcorr={1: 0}),
           'X8':  dict(kappa=2, c2=44, build=build_x8,  mus=[0, 1, 2], intcorr={0: 1, 2: 0})}
def normalize(S):
    """(offset, poly) -> (offset of first coefficient kept, list) with the offset shifted so that offset in (-1, 0]... keep raw"""
    return S
def get(name, M):
    d = MODELS2[name]
    out, h1 = d['build'](M)
    k = d['kappa']; chi2D = 8*k + 2*d['c2']; chi_2D = Fr(8*k, 6) + Fr(2*d['c2'], 12)
    return d, out, h1, chi2D, chi_2D
