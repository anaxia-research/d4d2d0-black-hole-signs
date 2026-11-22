"""Tables (LaTeX) and figure data for the paper."""
from fractions import Fraction as Fr
from flint import arb
from models3 import MODELS, derived
from forms3 import solve, offsets, polar_positions
from rig3 import Engine, A
P = '../logs/'
M = 80
def fr(x): return str(x) if x.denominator == 1 else r'%s\frac{%d}{%d}' % ('-' if x < 0 else '', abs(x.numerator), x.denominator)
def frs(x): return '$' + fr(x) + '$'
mcg = [n for n in MODELS if not n.startswith('V')]
# ---- models table
with open(P + 'tab_models.txt', 'w') as f:
    for name in mcg:
        d = MODELS[name]; chi, chiD, N, w = derived(d); off = offsets(d)
        f.write(r'%s & %d & %d & %d & %d & %d & %d & %d & %d & $%s$ & %s\\' % (d['label'].replace('Gr(2,7)/Z7', r'$\mathrm{Gr}/\Z_7$').replace('Reye/Z5', r'Reye$/\Z_5$'),
                d['b2'], d['K'], d['kappa'], d['c2'], d['tors'], chi, chiD, N, fr(w), ', '.join(frs(off[m]) for m in sorted(off))) + '\n')
# ---- solve all forms used later
forms = {}
for name in mcg:
    if name in ('G7', 'Rey'): continue
    d = MODELS[name]; pol = {0: [-18], 1: [486]} if name == 'B3' else None
    h, off = solve(name, M, pol); forms[d['label']] = (name, h, off, 1)
for (k, c2) in [(1, 10), (1, 22), (2, 8), (3, 6), (4, 4), (5, 2)]:
    name = 'V%d_%d' % (k, c2); d = MODELS[name]; chi, chiD, N, w = derived(d); eps = (-1)**(int(chiD) + 1)
    h, off = solve(name, M, {mu: ([1] if mu == 0 else []) for mu in range(k//2 + 1)})
    forms['U(%d,%d)' % (k, c2)] = (name, {mu: [eps*x for x in h[mu]] for mu in h}, off, eps)
# ---- thresholds table
with open(P + 'tab_thr.txt', 'w') as f:
    for lab, (name, h, off, eps) in forms.items():
        d = MODELS[name]; chi, chiD, N, w = derived(d)
        eng = Engine(d['K'], w, off, h, (-1)**(int(chiD) + 1))
        res = eng.analyse(M)
        cells = '; '.join('%d: %d/%d/%d/%d' % (mu, r['first'], r['n0'], r['nD'], r['nE']) for mu, r in res.items())
        lab2 = lab if lab.startswith('U') else lab
        f.write(r'%s & %s & %s & %s\\' % (('$U_{%s}$' % lab[1:]) if lab.startswith('U') else lab, frs(1 - w), eng.bmax.mid().str(3), cells) + '\n')
# ---- examples table
with open(P + 'tab_examples.txt', 'w') as f:
    for lab in ('(1,21)', '(1,8)', '(3,18)'):
        name, h, off, eps = forms[lab]; d = MODELS[name]; chi, chiD, N, w = derived(d)
        eng = Engine(d['K'], w, off, h, (-1)**(int(chiD) + 1))
        for n in range(1, 7):
            E = off[0] + n; L = eng.L(E); r, x = eng.R(E)
            rel = arb(int(h[0][n]))/L - 1
            f.write(r'%s & %d & %s & %d & %s & %s & %s\\' % (lab if n == 1 else '', n, frs(E), h[0][n], L.mid().str(8, radius=False),
                    rel.mid().str(3, radius=False), r.mid().str(3, radius=False)) + '\n')
        if lab != '(3,18)': f.write(r'\midrule' + '\n')
# ---- error figure data
import math
for lab, fn in (('(1,21)', 'err_q5'), ('(1,8)', 'err_r7'), ('(3,18)', 'err_m318')):
    name, h, off, eps = forms[lab]; d = MODELS[name]; chi, chiD, N, w = derived(d)
    eng = Engine(d['K'], w, off, h, 1)
    with open(P + fn + '.dat', 'w') as f:
        f.write('E err maj\n')
        for n in range(1, 41):
            E = off[0] + n; L = eng.L(E); r, x = eng.R(E)
            rel = abs(arb(int(h[0][n]))/L - 1)
            f.write('%.5f %.4f %.4f\n' % (float(E), math.log10(float(rel.mid())) if float(rel.mid()) > 0 else -99, math.log10(float(r.mid()))))
# ---- growth figure data: log10 c_0(n) vs sqrt(E) for the six U forms
with open(P + 'grow_u.dat', 'w') as f:
    labs = ['U(1,10)', 'U(1,22)', 'U(2,8)', 'U(3,6)', 'U(4,4)', 'U(5,2)']
    f.write('n ' + ' '.join('s%d y%d' % (i, i) for i in range(6)) + '\n')
    for n in range(1, 61):
        row = [str(n)]
        for lab in labs:
            name, h, off, eps = forms[lab]
            E = off[0] + n
            row += ['%.5f' % math.sqrt(float(E)), '%.5f' % math.log10(float(h[0][n]))]
        f.write(' '.join(row) + '\n')
# ---- entropy constants
with open(P + 'tab_entropy.txt', 'w') as f:
    for lab, (name, h, off, eps) in forms.items():
        d = MODELS[name]; chi, chiD, N, w = derived(d); nu = 1 - w; K = d['K']
        D0 = Fr(chi, 24); c00 = abs(h[0][0])
        C = arb(int(c00))*A(D0)**((2*A(nu) - 1)/4)/(2*arb(K)).sqrt()
        f.write(r'%s & %d & %s & %s & %s & %s & %s\\' % (('$U_{%s}$' % lab[1:]) if lab.startswith('U') else lab, chi, frs(nu), frs(-(nu/2 + Fr(1, 4))),
                frs((4*nu*nu - 1)/8), C.mid().str(6, radius=False), (C.log()).mid().str(6, radius=False)) + '\n')
# ---- first coefficients (appendix)
with open(P + 'tab_first.txt', 'w') as f:
    for lab, (name, h, off, eps) in forms.items():
        if not lab.startswith('U'): continue
        for mu in sorted(h):
            f.write(r'$U_{%s}$ & %d & %s & %s\\' % (lab[1:], mu, frs(off[mu]), ', '.join(str(eps*h[mu][n]) for n in range(7))) + '\n')
with open(P + 'tab_first2.txt', 'w') as f:
    for (k, c2) in [(1, 34), (1, 46), (2, 20), (3, 18), (11, 2)]:
        name = 'V%d_%d' % (k, c2); d = MODELS[name]; off = offsets(d); Pp = polar_positions(d, off)
        for i, gl in ((0, 'G'), (1, 'G\'')):
            pol = {mu: [0]*sum(1 for (m, n) in Pp if m == mu) for mu in off}
            m, n = Pp[i]; pol[m][n] = 1
            h, off = solve(name, 8, pol)
            for mu in sorted(h):
                if k == 11 and mu not in (0, 2, 3, 5): continue
                f.write(r'$%s_{%d,%d}$ & %d & %s & %s\\' % (gl, k, c2, mu, frs(off[mu]), ', '.join(str(h[mu][n]) for n in range(6))) + '\n')
print(open(P + 'tab_thr.txt').read()); print(open(P + 'tab_examples.txt').read()); print(open(P + 'tab_entropy.txt').read())
