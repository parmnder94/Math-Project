# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds nothing but Python, NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, the constant K to about 1e-7 relative.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Smooth-profile anomaly formula | The local formula (density ~ \|grad T\|^2/T) has no finite value for a step; there is no number to submit. |
| Monte Carlo plus extrapolation in m | Even with an exact Ornstein-Uhlenbeck step, the simulation reaches about 1% at m = 0.0025, and it still carries an O(sqrt m) bias. The gate is 1e-4 relative. A plain Euler scheme is about 1% high at every mass. |
| One-dimensional layer, frozen drag frame | See `ablations.json`: 35% and 3% off. |
| Looking up a constant | The layer constant for this drag tensor, field and temperature ratio is not tabulated anywhere; the instruction cites nothing. |
| Strings, NaN, Infinity, booleans | Values are parsed as finite floats (or p/q strings); booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |

The verifier never executes agent code. It parses one JSON file.
