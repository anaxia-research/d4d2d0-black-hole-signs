# Signs of D4-D2-D0 black hole indices

Code and computation logs for the series *Signs of D4-D2-D0 Black Hole Indices*.

The Abelian D4-D2-D0 indices of a one-parameter Calabi-Yau threefold count the microstates of
Maldacena-Strominger-Witten black holes with one unit of D4-brane charge. Their generating series is a
vector-valued modular form of weight -3/2. This repository proves, for every charge, the sign, the growth and
the log-concavity of these indices, and the black hole entropy with a certified error.

## Part I: the hypergeometric threefolds (`part1/`)

Models: X5, X6, X8, X10, X4,3, X6,4, X3,3, X4,4, X6,6, X6,2, X4,2. Polar coefficients from
Alexandrov-Feyzbakhsh-Klemm-Pioline-Schimannek, *Quantum geometry, stability and modularity*,
Commun. Number Theory Phys. 18 (2024) 49-151.

| Result | What is computed | Log |
|---|---|---|
| Theorem A (exact construction) | each series as an exact rational combination of `E4^a E6^b D^l theta / eta^N`; agreement with all 161 published non-polar coefficients | `construction.log` |
| S-matrix | `h(-1/tau) = (-i tau)^{-3/2} eps kappa^{-1/2} cos(2 pi mu nu/kappa) h(tau)`, `eps = (-1)^{chi_D+1}`, checked to 1e-29 | `smatrix.log` |
| Theorem B (Rademacher with certified error) | `|c/L - 1| < R` certified at 2097 pairs | `certify_error.log`, `rademacher_c1.log` |
| Theorem C (signs) | `R < 1` certified from the threshold on; exact signs below it; 18 negative first non-polar coefficients | `rigorous.log`, `tab_exc.txt` |
| Theorem D (growth) | analytic threshold certified in ball arithmetic; exact checks below it | `rigorous.log` |
| Theorem E (log-concavity) | analytic threshold certified; exact checks below it | `logconcave.log` |
| Extended data | integrality, positivity, growth, log-concavity for n <= 200; higher-order Turan inequalities | `extended_200.log`, `turan.log` |
| Tables and figure data | thresholds, constants, worked examples, error curves | `tab_*.txt`, `*.dat` |

All exact computations use rational arithmetic (FLINT); all real inequalities are certified in
midpoint-radius ball arithmetic (Arb, 256 bits).

### Reproduce

```sh
pip install python-flint==0.9.0 mpmath
./part1/run_all.sh      # about 15 seconds
```

```mermaid
flowchart TD
  P[Polar data] --> A[Theorem A: exact modular forms]
  A --> S[Weil representation, pure D4 term]
  S --> B[Theorem B: Rademacher with certified error R]
  Bes[Bessel lemmas] --> B
  B --> C[Theorem C: signs]
  B --> D[Theorem D: growth]
  B --> E[Theorem E: log-concavity]
  C --> F[Theorem F: entropy]
  D --> F
  E --> F
```

## License

MIT (see `LICENSE`).
