"""Higher-order Turan inequality 4(a_n^2-a_{n-1}a_{n+1})(a_{n+1}^2-a_n a_{n+2}) > (a_n a_{n+1}-a_{n-1}a_{n+2})^2: exact check, n <= 200."""
from vvmf import solve
from models import MODELS
M = 200
for name in MODELS:
    d, mus, off, h, x, l = solve(name, M + 12)
    out = []
    for mu in mus:
        first = next(n for n in range(M) if off[mu] + n >= 0)
        a = [int(h[mu][n].p) for n in range(M + 3)]
        fails = [n for n in range(first + 1, M) if not 4*(a[n]**2 - a[n-1]*a[n+1])*(a[n+1]**2 - a[n]*a[n+2]) > (a[n]*a[n+1] - a[n-1]*a[n+2])**2]
        out.append(f"mu={mu}: fails {fails[:12]}{'...' if len(fails) > 12 else ''} (last {fails[-1] if fails else '-'})")
    print(name, "; ".join(out), flush=True)
