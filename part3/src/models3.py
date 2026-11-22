"""Data of the quotient threefolds of McGovern, "New examples of Abelian D4D2D0 indices",
SciPost Phys. 20 (2026) 015 (arXiv:2412.01149), Table 1, Sections 3.7, 4 and 5.

For a model: b2 = second Betti number, K = |Lambda*/Lambda| (order of the discriminant group),
chi = chi(D) = kappa_111 + c_{2,1} (Euler number of the divisor), N = power of 1/eta in the basis
eta^{-N} P(E4,E6) D^l theta^{(K)}, N = 4K + c_{2,1} + 1 - b2, tors = order of H^2(Y,Z)_tors.
polar[mu] = the polar coefficients of h_mu, most polar first (from formula (3.44) of McGovern,
or from (3.54) for the two obstructed models); check[mu] = further coefficients printed by McGovern
(as coefficients of h_mu, i.e. including the overall prefactor).  Components mu = 0..floor(K/2)."""

MODELS = {
 # one-parameter models (b2 = 1)
 'Q5':  dict(label='(1,21)', geom=r'$\mathbb P^4[5]/\mathbb Z_5$', b2=1, K=1, kappa=1, c2=10, tors=5, chiY=-40,
             polar={0: [25]},
             check={0: [25*x for x in [251, 4877, 49378, 360005]]}),
 'B3':  dict(label='(1,25)', geom=r'$\mathbb P^5[3,3]/\mathbb Z_3$', b2=1, K=3, kappa=3, c2=18, tors=3, chiY=-48,
             polar={0: [-18], 1: [None]},      # second polar term unknown (= t)
             check={0: [9*x for x in [174, 119052, 5318746, 117995460]],
                    1: [9*x for x in [26838, 1669194, 44349552, 738224496]]}, t_mcg=486),
 'R7':  dict(label='(1,8)', geom=r'Pfaffian$/\mathbb Z_7$', b2=1, K=2, kappa=2, c2=8, tors=7, chiY=-14,
             polar={0: [49], 1: []},
             check={0: [49*x for x in [136, 2081, 18152, 117028]], 1: [49*x for x in [56, 1136, 10912, 75072, 414304]]}),
 'H5':  dict(label='(1,6)', geom=r'symmetroid$/\mathbb Z_5$', b2=1, K=2, kappa=2, c2=8, tors=5, chiY=-10,
             polar={0: [25], 1: []},
             check={0: [25*x for x in [136, 2081, 18152, 117028]], 1: [25*x for x in [56, 1136, 10912, 75072, 414304]]}),
 # multiparameter models
 'M229': dict(label='(2,29)', geom=r'$\mathbb P^2\times\mathbb P^2[3,3]/\mathbb Z_3$', b2=2, K=1, kappa=0, c2=12, tors=3,
             polar={0: [9]}, check={0: [9*x for x in [252, 5130, 54760, 419895, 2587788, 13630694]]}),
 'M220': dict(label='(2,20)', geom=r'CICY$/\mathbb Z_3$', b2=2, K=1, kappa=0, c2=12, tors=3,
             polar={0: [9]}, check={0: [9*x for x in [252, 5130, 54760, 419895, 2587788, 13630694]]}),
 'M212': dict(label='(2,12)', geom=r'$\mathbb P^4\times\mathbb P^4[1^5,1^5]/\mathbb Z_5$', b2=2, K=2, kappa=1, c2=10, tors=5,
             polar={0: [25], 1: []},
             check={0: [25*x for x in [137, 2219, 20508, 139755, 779254]], 1: [25*x for x in [56, 1192, 12160, 88424, 514888, 2564184]]}),
 'M318': dict(label='(3,18)', geom=r'$(\mathbb P^2)^3[3,3,3]/\mathbb Z_3$', b2=3, K=2, kappa=0, c2=12, tors=3,
             polar={0: [9], 1: []},
             check={0: [9*x for x in [138, 2358, 23004, 165117, 967374]], 1: [9*x for x in [56, 1248, 13464, 103136, 631488, 3298752]]}),
 'M315': dict(label='(3,15)', geom=r'CICY$/\mathbb Z_3$', b2=3, K=2, kappa=0, c2=12, tors=3,
             polar={0: [9], 1: []},
             check={0: [9*x for x in [138, 2358, 23004, 165117, 967374]], 1: [9*x for x in [56, 1248, 13464, 103136, 631488, 3298752]]}),
 # the two models whose predicted polar terms (McGovern (3.55)) we show cannot be completed
 'G7':  dict(label='Gr(2,7)/Z7', geom=r'$\mathrm{Gr}(2,7)\cap\mathbb P^{13}/\mathbb Z_7$', b2=1, K=6, kappa=6, c2=12, tors=7,
             polar={0: [-98], 1: [882], 2: [343], 3: [0]}, check={0: [], 1: [], 2: [], 3: []}),
 'Rey': dict(label='Reye/Z5', geom=r'Reye congruence$/\mathbb Z_5$', b2=1, K=7, kappa=7, c2=10, tors=10,
             polar={0: [-200], 1: [1600], 2: [500], 3: [0]}, check={0: [], 1: [], 2: [], 3: []}),
}

def derived(d):
    """chi(D), chi_D, N, weight."""
    from fractions import Fraction as Fr
    K, b2, c2 = d['K'], d['b2'], d['c2']
    # chi(D) = kappa_111 + c_{2,1} (kappa = kappa_111 for p = e_1)
    chi = d['kappa'] + c2
    chiD = Fr(d['kappa'], 6) + Fr(c2, 12)
    N = 4*K + c2 + 1 - b2
    w = -1 - Fr(b2, 2)
    return chi, chiD, N, w

# All b2 = 1 invariants (kappa, c2) with chi_D = kappa/6 + c2/12 a positive integer and at most two polar
# terms (enumerated in classify3.py).  'U' entries: one polar term; 'P' entries: two polar terms, the
# second polar coefficient being a free parameter.  Polar coefficients normalised: pure D4 term = 1.
def virtual(kappa, c2):
    return dict(label='(%d,%d)' % (kappa, c2), geom='', b2=1, K=kappa, kappa=kappa, c2=c2, tors=1,
                polar={mu: [] for mu in range(kappa//2 + 1)}, check={mu: [] for mu in range(kappa//2 + 1)})
for k, c in [(1, 10), (1, 22), (2, 8), (3, 6), (4, 4), (5, 2), (1, 34), (1, 46), (2, 20), (3, 18), (11, 2)]:
    MODELS['V%d_%d' % (k, c)] = virtual(k, c)
