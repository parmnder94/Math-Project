# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds only Python, NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, the exact fraction.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Direct simulation plus rounding | Corrections decay like N^(-1/8): at N = 1e16 the error is still 1e-3, so simulation cannot single out the fraction (gate: exact, or 1e-12 as a decimal). |
| Plausible shortcuts | See `ablations.json`: endpoint frame dropped (0.769), connection omitted (0) or only -pi K_r (0.906), no phase reversal in the reversed pulse (0.905), carrier phases ignored (0.905). |
| Lower-order Floquet theory | Orders 1-7 cancel identically by construction; any truncation below eighth order leaves no coupling at all. |
| Looking up the reference | The paper gives the general high-frequency framework, not this protocol or number. |
| Strings, NaN, Infinity, booleans | Values are parsed as exact fractions or finite floats; booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |
