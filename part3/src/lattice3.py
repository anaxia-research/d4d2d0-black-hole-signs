"""The lattice (H^2(Y,Z)_free, k_1jk) of the multiparameter models is isomorphic to <K> + <-1>^{b2-1}."""
import itertools, numpy as np
G = {'(2,29)': [[0,1],[1,1]], '(2,20)': [[0,1],[1,3]], '(2,12)': [[1,2],[2,2]],
     '(3,18)': [[0,1,1],[1,1,2],[1,2,1]], '(3,15)': [[0,1,1],[1,1,3],[1,3,3]]}
K = {'(2,29)': 1, '(2,20)': 1, '(2,12)': 2, '(3,18)': 2, '(3,15)': 2}
out = open('../logs/lattice3.log', 'w')
for name, g in G.items():
    g = np.array(g); r = len(g)
    target = np.diag([K[name]] + [-1]*(r - 1))
    found = None
    for ent in itertools.product(range(-2, 3), repeat=r*r):
        P = np.array(ent).reshape(r, r)
        if round(abs(np.linalg.det(P))) != 1: continue
        if (P.T @ g @ P == target).all(): found = P; break
    line = "%s: Gram %s ~ diag%s via P = %s" % (name, g.tolist(), target.diagonal().tolist(), None if found is None else found.T.tolist())
    print(line); print(line, file=out)
