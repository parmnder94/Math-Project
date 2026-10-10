# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds Python, NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, the exact fraction.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Numerical estimate written as a decimal | Only an exact fraction string is accepted; decimals raise. |
| Rationalising a numerical estimate | Monte Carlo of the SDE reaches about 1e-2 on the finite heat (large O(m) corrections, a 1/m divergence to subtract), far from what fixes a fraction with a 14-digit denominator. Exact theory is needed. |
| Reusing the answer of the original (author-supplied) problem | -141231719/146966400 belongs to different parameters; rejected (`local_runs/verifier_runs.txt`). |
| Plausible shortcuts | Kinetic current in the relative entropy (+0.441), x as a 1D Markov diffusion (+0.740), density proportional to T (-1.129), current part of the force dropped in the heat (-1.919), both reservoirs weighted by 1/T (+0.620), observed path taken as reversible (+2.542); see `ablations.json`. |
| Pattern-matching the structure from the statement | The force is printed only in expanded trigonometric form; the fast covariance, the current h, the density exponent and the divergence coefficient kappa are not printed anywhere. |
| Looking up the problem | The model and its numbers are new; the cited small-mass framework is general. |
| Strings with spaces, booleans, NaN | Whitespace is stripped; booleans, floats and non-fraction strings raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |
