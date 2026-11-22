"""Extended exact data (n <= M): integrality, positivity beyond the first non-polar term,
strict monotonicity and log-concavity c(n)^2 > c(n-1) c(n+1)."""
from vvmf import solve
from models import MODELS
import sys
M = int(sys.argv[1]) if len(sys.argv) > 1 else 200
for name in MODELS:
    d, mus, off, h, x, labels = solve(name, M + 12)
    integ = all(h[mu][n].q == 1 for mu in mus for n in range(M + 1))
    rep = []
    for mu in mus:
        first = next(n for n in range(M) if off[mu] + n >= 0)
        c = [int(h[mu][n].p) for n in range(M + 1)]
        pos = all(c[n] > 0 for n in range(first + 1, M + 1))
        mono = all(c[n + 1] > c[n] for n in range(first, M))
        lc_fail = [n for n in range(first + 1, M) if not c[n]**2 > c[n-1]*c[n+1]]
        lc_fail2 = [n for n in range(first + 2, M) if not c[n]**2 > c[n-1]*c[n+1]]
        rep.append((mu, pos, mono, lc_fail[:6], len(str(c[M]))))
    print(name, "integral:", integ, " ".join(f"[mu={mu} pos>first:{p} mono:{m} logconc-fail:{l} digits(c_M):{dg}]" for mu, p, m, l, dg in rep), flush=True)
