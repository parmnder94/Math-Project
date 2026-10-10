# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds Python with NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, alpha_zz to 10 digits.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Compute the axion angle and report theta/2pi | That is only the Chern-Simons part: 0.2374, which is 0.036 (12 times the bound) from alpha_zz = 0.2014. |
| Surface (layer-resolved) Hall conductivity of a slab | Also reproduces theta/2pi = 0.2374; the cross-gap part is a bulk polarization that the layer decomposition misses. |
| Half-quantized value (0.5) or zero | Time reversal, inversion and their product are broken; neither applies. |
| Unconverged finite-field supercell | The forward difference at flux 1/32 gives 0.155; central differences at 1/32 to 1/64 miss by 7e-3 to 2e-3 until extrapolated. |
| Sign of the Peierls phase or of the orientation reversed | Gives -0.2014. |
| Looking up the reference | The general theory (Essin, Turner, Moore and Vanderbilt 2010; Malashevich, Souza, Coh and Vanderbilt 2010) is public; this model, its parameters and its value are not. |
| Strings, NaN, Infinity, booleans | Values are parsed as finite floats; booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |
