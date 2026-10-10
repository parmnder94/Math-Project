# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds only Python, NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, the limit to 30 digits.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Direct simulation plus extrapolation | Corrections decay like N^(-1/8): at N = 1e16 the error is still 1e-3, and extrapolation in eps reaches only about 3e-5. The gate is 1e-9 absolute. |
| Reverse-engineering the eighth-order coefficients from the constants | The prompt uses a round amplitude scaling (Delta T_N eps_N^8 = 1/10) and no tuned gap shift, so no constant encodes B_8; a relative error of 1e-6 in B_8 already moves p by 8e-8. |
| Plausible shortcuts | See `ablations.json`: endpoint frame dropped (0.474), connection omitted (0.493) or only -pi K_r (0.749), no phase reversal in the reversed pulse (0.745), carrier phases ignored (0.724), (Q-P) part of B_8 omitted (0.981). |
| Lower-order Floquet theory | Orders 1-7 cancel identically; any truncation below eighth order leaves no coupling (p = 1). |
| Looking up the reference | The cited framework is general; neither this protocol nor its number appears in the literature. |
| Strings, NaN, Infinity, booleans | Values are parsed as finite floats; booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |
