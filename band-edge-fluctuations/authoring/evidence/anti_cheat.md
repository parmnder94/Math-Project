# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only (model, parameters, output format). The environment image holds nothing but Python, NumPy, SciPy, SymPy and mpmath.
What stays sealed in `tests/`: `truth.json`, the four constants to about 1e-9 relative.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`, so every test errors and the reward is 0. |
| Monte Carlo of finite networks | At N = 200 with 80,000 samples, m carries a statistical error of about 0.02 (about 20%), v about 0.5% and k more than 100%. There is also O(1/N) bias. The gate is 1e-6 relative. |
| Brute-force class enumeration without resummation | Near the band edge the class contributions decay like 0.94^l, so enumeration to length 10 is still 81% off in v and 6% off in m (`ablations.json`). |
| Leading order only (m = 0), Gaussian fluctuations (k = 0), or a wrong normalisation | Relative errors of order 1. |
| Guessing, or copying from the literature | Neither the model nor the parameters appear elsewhere, and none of the constants has a closed form. |
| Strings, NaN, Infinity, booleans, extra keys | Keys are normalised and values parsed with `float`. Non-finite values and booleans raise. Extra keys are ignored. |
| Looking for the truth in the environment | There is none. `environment/` contains only a Dockerfile. |
| Writing into `/logs/verifier` from the agent container | The verifier runs in a separate container and reads only the declared artifact. |
| Partial theory (no power-cycle structure, unit weights for powers, resummation only up to total power 3) | See `ablations.json`: every one misses the 1e-6 gate. |

The verifier never executes agent code. It parses one JSON file.
