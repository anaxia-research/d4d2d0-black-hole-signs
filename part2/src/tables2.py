"""LaTeX tables and figure data for Part II (written to ../logs)."""
from fractions import Fraction as Fr
from flint import arb
from rigorous2 import setup, Rtot, Lhat, Mstar, Istar, R, q2a, PI
from signs2 import intseries, F
from signs2c import Rpoint
from vvmf import solve
from r2common import MODELS2
import math, os
P = '../logs/'
def fr(x):
    x = Fr(x)
    if x.denominator == 1: return f"${x.numerator}$"
    s = '-' if x < 0 else ''
    return f"${s}\\frac{{{abs(x.numerator)}}}{{{x.denominator}}}$"
def num(x):
    x = Fr(x)
    if x.denominator == 1: return f"${x.numerator}$"
    return f"${x.numerator}/{x.denominator}$"
SS = {n: setup(n, M=420) for n in ('X10', 'X8')}
CHIY = {'X10': -288, 'X8': -296}
# models table
rows = []
for n, S in SS.items():
    d = S['d']; k = d['kappa']; c2 = d['c2']
    rows.append(f"$X_{{{n[1:]}}}$ & ${k}$ & ${c2}$ & ${CHIY[n]}$ & ${k+c2}$ & ${Fr(k,6)+Fr(c2,12)}$ & ${8*k+2*c2}$ & ${S['chi']}$ & ${'+1' if (-1)**(int(S['chi'])+1) > 0 else '-1'}$ & {fr(S['D0'])} & {fr(2*S['r1']['D1'])} & ${float(S['delta'].mid()):.4f}$\\\\")
open(P + 'tab_models.txt', 'w').write('\n'.join(rows) + '\n')
# polar table
rows = []
for n, S in SS.items():
    for mu in S['mus']:
        tot, hh, mk = S['out'][mu][2], S['out'][mu][0], S['out'][mu][1]
        o = tot[0]; first = True
        for i in range(tot[1].degree() + 1):
            E = o + i
            if E >= 0: break
            ct = F(tot[1][i]); ch = F(hh[1][int(E - hh[0])]) if E >= hh[0] else Fr(0); cm = F(mk[1][int(E - mk[0])]) if E >= mk[0] else Fr(0)
            if ct == 0 and ch == 0 and cm == 0 and E < -S['D0'] + 0: continue
            lab = f"$X_{{{n[1:]}}}$ & ${mu}$" if first else " & "
            first = False
            rows.append(f"{lab} & {fr(-E)} & {num(ct)} & {num(ch)} & {num(cm)}\\\\")
        rows.append("\\addlinespace[2pt]")
open(P + 'tab_polar.txt', 'w').write('\n'.join(rows[:-1]) + '\n')
# exceptions + first coefficients
exc = []; firsts = []
for n, S in SS.items():
    for mu in S['mus']:
        ser = intseries(S, mu, 30)
        f0 = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
        for j in range(f0, f0 + 6):
            E, c, ci = ser[j]
            if not c > 0 or not ci > 0:
                exc.append(f"$X_{{{n[1:]}}}$ & ${mu}$ & $n^\\ast+{j-f0}$ & {fr(-E)} & {num(c)} & {num(ci)}\\\\")
        for j in range(f0, f0 + 4):
            E, c, ci = ser[j]
            lab = f"$X_{{{n[1:]}}}$ & ${mu}$" if j == f0 else " & "
            firsts.append(f"{lab} & {fr(E)} & {num(c)} & {num(ci) if ci != c else '$=$'}\\\\")
        firsts.append("\\addlinespace[2pt]")
open(P + 'tab_exc.txt', 'w').write('\n'.join(exc) + '\n')
open(P + 'tab_first.txt', 'w').write('\n'.join(firsts[:-1]) + '\n')
# thresholds (values from the logs of signs2c / growth2)
import signs2c, growth2, io
th = []
for n, S in SS.items():
    r0 = signs2c.run(n, io.StringIO()); S2, g = growth2.run(n, io.StringIO())
    for mu in S['mus']:
        lab = f"$X_{{{n[1:]}}}$" if mu == 0 else ""
        ser = intseries(S, mu, 12); f0 = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
        th.append(f"{lab} & ${mu}$ & {fr(ser[f0][0])} & ${r0[mu][0]}$ & ${r0[mu][3]}$ & ${g[mu]['n2']}$ & ${g[mu]['n3']}$\\\\")
open(P + 'tab_thr.txt', 'w').write('\n'.join(th) + '\n')
# constants
cons = []
for n, S in SS.items():
    K = S['r1']['K']
    cons.append(f"$X_{{{n[1:]}}}$ & ${float(S['bmax'].mid()):.4f}$ & ${float(S['Kg'].mid()):.4f}$ & ${', '.join(format(math.ceil(float(v.upper())*1000)/1000, '.3f') for v in K.values())}$ & ${float(((6/S['delta'])**2).mid()):.2f}$ & ${float((arb(S['chi'].numerator)*(8*S['d']['kappa']+2*S['d']['c2'])/(24*arb(4*S['d']['kappa']).sqrt())).mid()):.4f}$\\\\")
open(P + 'tab_const.txt', 'w').write('\n'.join(cons) + '\n')
# examples X10 mu=0 and X8 mu=2
ex = []
for n, mu in (('X10', 0), ('X8', 2)):
    S = SS[n]; ser = intseries(S, mu, 40); f0 = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
    for j in list(range(f0, f0 + 7)) + [f0 + 10, f0 + 20, f0 + 30]:
        E, c, ci = ser[j]; Ea = q2a(E)
        lh = Lhat(S['L'], q2a(S['chi']), S['D0'], Ea)
        m = F(S['out'][mu][1][1][int(E - S['out'][mu][1][0])])
        rp = Rpoint(S, mu, E, integer=False)
        lab = f"$X_{{{n[1:]}}}$, $\\mu={mu}$" if j == f0 else ""
        ex.append(f"{lab} & ${j-f0}$ & {fr(E)} & ${float(c):.4e}$ & ${float(lh.mid()):.4e}$ & ${float((q2a(c)/lh - 1).mid()):+.2e}$ & ${float((q2a(m)/lh).mid()):+.2e}$ & {('$' + format(float(rp.upper()), '.2e') + '$') if rp is not None else '--'}\\\\")
    ex.append("\\addlinespace[2pt]")
open(P + 'tab_examples.txt', 'w').write('\n'.join(ex[:-1]) + '\n')
# figure data
for n, mu, nmax, tag in (('X10', 0, 160, 'x10'), ('X8', 0, 110, 'x8')):
    S = SS[n]; ser = intseries(S, mu, nmax); f0 = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
    with open(P + f'err2_{tag}.dat', 'w') as fh:
        fh.write("E tot majt mock majm hat majh\n")
        for j in range(f0, len(ser)):
            E, c, ci = ser[j]; Ea = q2a(E)
            lh = Lhat(S['L'], q2a(S['chi']), S['D0'], Ea)
            ch = F(S['out'][mu][0][1][int(E - S['out'][mu][0][0])]); cm = F(S['out'][mu][1][1][int(E - S['out'][mu][1][0])])
            rh, _ = R(S['T'], S['D0'], Ea)
            rp = Rpoint(S, mu, E, integer=False)
            if rp is None: continue
            l10 = lambda a: float(a.mid().log())/math.log(10)
            fh.write(f"{float(E):.4f} {l10(abs(q2a(c)/lh - 1)):.4f} {l10(rp):.4f} {l10(abs(q2a(cm))/lh):.4f} {l10(rp - rh):.4f} {l10(abs(q2a(ch)/lh - 1)):.4f} {l10(rh):.4f}\n")
# anatomy: log10 |hat c| and |m| against sqrt(E), X10 mu=0
S = SS['X10']; mu = 0; ser = intseries(S, mu, 400); f0 = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
with open(P + 'anat_x10.dat', 'w') as fh:
    fh.write("s hat mock\n")
    for j in range(f0, len(ser), 4):
        E = ser[j][0]
        ch = F(S['out'][mu][0][1][int(E - S['out'][mu][0][0])]); cm = F(S['out'][mu][1][1][int(E - S['out'][mu][1][0])])
        l = lambda v: math.log10(abs(v.numerator)) - math.log10(v.denominator)
        fh.write(f"{math.sqrt(float(E)):.4f} {l(ch):.4f} {l(cm):.4f}\n")
# threshold plot: R_tot for all 5 components on E in [0,20]
with open(P + 'thr_all.dat', 'w') as fh:
    cols = []
    data = {}
    for n, S in SS.items():
        for mu in S['mus']:
            ser = intseries(S, mu, 30)
            pts = []
            for E, c, ci in ser:
                if 0 < E <= 20:
                    rp = Rpoint(S, mu, E, integer=True)
                    if rp is not None: pts.append((float(E), float(rp.mid().log())/math.log(10)))
            data[f"{n}m{mu}"] = pts
    for key, pts in data.items():
        with open(P + f'thr_{key}.dat', 'w') as g:
            g.write("E r\n" + ''.join(f"{a:.4f} {b:.4f}\n" for a, b in pts))
print("tables written")
