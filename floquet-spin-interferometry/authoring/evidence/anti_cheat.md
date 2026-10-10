# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds only Python, NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, with p_inf and c1 to 30 digits and c2, c3 to 18.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Direct simulation plus extrapolation | Corrections decay like N^(-1/8): at N = 1e16 the error in p is still 1e-3, and extrapolation in eps reaches only about 3e-5. For c, the second-order coefficient is about -15, so the slope (p_N - p_inf)/eps_N is off by 1e-4 even at N = 1e40. The gate is 1e-9 absolute on every field. |
| Reverse-engineering the eighth-order coefficients from the constants | The prompt uses a round amplitude scaling (Delta T_N eps_N^8 = 1/10) and no tuned gap shift, so no constant encodes B_8; a relative error of 1e-6 in B_8 already moves p by 8e-8. |
| Plausible shortcuts | See `ablations.json`: endpoint frame dropped (0.474), connection omitted (0.493) or only -pi K_r (0.749), no phase reversal in the reversed pulse (0.745), carrier phases ignored (0.724), (Q-P) part of B_8 omitted (0.981). |
| Fitting simulated p_N for c2, c3 | The coefficients grow quickly (c4 is about 421, c5 about -58): the remainder (p_N - p_inf - c1 eps)/eps^2 at N = 1e40 still misses c2 by 4e-5, and c3 needs a further division by eps. The gate is 1e-9 absolute on each. |
| Effective-theory shortcuts for c2, c3 | Static first-order kicks around exp(-i G_r) give c2 = -9.160, c3 = 1.180; the exact Floquet Hamiltonian with a connection not dressed by the micromotion gives c1 = +0.143, c2 = -1.855, c3 = -1.455 (`ablations.json`). |
| Wrong treatment of the first correction | Stroboscopic gauge without the boundary kick gives c = 0; a kick at only one end gives +0.585 or -0.669; dropping the 2Q gap in the kick denominators gives -0.420 (`ablations.json`). |
| Lower-order Floquet theory | Orders 1-7 cancel identically; any truncation below eighth order leaves no coupling (p = 1). |
| Looking up the reference | The cited framework is general; neither this protocol nor its number appears in the literature. |
| Strings, NaN, Infinity, booleans | Values are parsed as finite floats; booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |
