# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds nothing but Python, NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, the exact fraction.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Monte Carlo plus rounding to a fraction | Simulation reaches about 1% on the heat rate, and nothing on the path relative entropy, so it cannot single out the exact fraction. |
| Plausible shortcuts | See `ablations.json`: 19/210, 0, 19/400 and 724/16575 are all rejected by exact comparison. |
| Looking up the source paper | The paper gives general formulas, not this instance; the drag rotation, the field and the fixed-field comparison must still be worked out. The instruction does not cite it. |
| Strings, NaN, Infinity, booleans | Values are parsed as exact fractions or finite floats; booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |

The verifier never executes agent code. It parses one JSON file.
