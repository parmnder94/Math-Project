# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only (model, parameters, output format). The environment image holds nothing but Python, NumPy, SciPy, SymPy and mpmath.
What stays sealed in `tests/`: `truth.json`, the five constants to about 1e-8 relative or better.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`, so every test errors and the reward is 0. |
| Monte Carlo of finite graphs | The standard deviation of X_N is about 1.7. The 1/N effect of m1 is about 1e-3 of it at N = 10, so m1 is out of reach, and v1 is resolved only to tens of percent. The gate is 1e-6 relative on every key. |
| Exact averages at small N plus extrapolation | Brute force over S_N × S_N is possible only up to N ≈ 9. Extrapolating the 1/N series from there gives m1 and v1 to about 1e-3, three orders of magnitude short of the gate. |
| Leading order only (m1 = v1 = 0), or a wrong normalisation | Relative errors of order 1. |
| Partial 1/N theory (critical-subgroup count without the falling-factorial term, no disconnected cross term, no −E_aE_b on connected cycles, no 1/N term for proper powers) | See `ablations.json`: every route misses the 1e-6 gate, most by orders of magnitude. |
| Stopping the pair sums too early | Pairs to total length 10 already miss the gate (1.5e-6 on v1). |
| Guessing, or copying from the literature | Neither the model nor the parameters appear elsewhere, and none of the constants has a closed form. |
| Strings, NaN, Infinity, booleans, extra keys | Keys are normalised and values parsed with `float`. Non-finite values and booleans raise. Extra keys are ignored. |
| Looking for the truth in the environment | There is none. `environment/` contains only a Dockerfile. |
| Writing into `/logs/verifier` from the agent container | The verifier runs in a separate container and reads only the declared artifact. |

The verifier never executes agent code. It parses one JSON file.
