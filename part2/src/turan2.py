"""Higher-order Turan inequality (hyperbolicity of the cubic Jensen polynomial) for the rank-2 sequences."""
from rigorous2 import setup
from signs2 import intseries
import sys
def run(name, M=400, out=sys.stdout):
    S = setup(name, M=M)
    for mu in S['mus']:
        ser = intseries(S, mu, M)
        first = next(i for i, (E, c, ci) in enumerate(ser) if E >= 0)
        for kind, j in (('rational', 1), ('integer', 2)):
            a = [ser[i][j] for i in range(len(ser))]
            fails = []
            for n in range(first + 1, len(a) - 2):
                lhs = 4*(a[n]**2 - a[n-1]*a[n+1])*(a[n+1]**2 - a[n]*a[n+2]); rhs = (a[n]*a[n+1] - a[n-1]*a[n+2])**2
                if not lhs > rhs: fails.append(n - first)
            print(f"{name} mu={mu} {kind}: higher-order Turan fails for n-n* in {fails or 'none'} (checked n-n* <= {len(a)-3-first})", file=out)
if __name__ == "__main__":
    with open('../logs/turan2.log', 'w') as f:
        for name in ('X10', 'X8'): run(name, out=f)
    print(open('../logs/turan2.log').read())
