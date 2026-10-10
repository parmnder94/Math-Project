# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image holds Python with NumPy, SciPy, SymPy and mpmath. What stays sealed in `tests/`: `truth.json`, theta to 12 digits.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Quantized value (theta = pi or 0) | Time reversal, inversion and their product are all broken; theta = 1.885. |
| Closed-form d-vector formula (Neel mass only) | Valid only without the orbital-dependent exchange; gives 2.345. |
| Eigenvectors straight from a diagonalizer | Random phases make the Chern-Simons integrand meaningless; a smooth periodic gauge or a gauge-free second-Chern route is required. |
| Projection onto the orbital basis states | The overlap becomes singular (smallest eigenvalue 1e-5 on a 64^3 grid) because of the band inversion. |
| Gauge-free route along the first path one tries (exchange before Neel mass) | The gap closes (Weyl semimetal) on that path. |
| Coarse uniform grid | 64^3, 96^3 and 128^3 miss by 2e-3, 1e-4 and 1e-6 (`ablations.json`); the gate is 1e-8. |
| Sign or branch slip | -1.885 or values shifted by 2 pi are rejected; the interval (-pi, pi] is fixed by the instruction. |
| Looking up the reference | The cited framework (Qi-Hughes-Zhang 2008; Essin-Moore-Vanderbilt 2009) defines theta in general; this parameter set and its value do not appear in the literature. |
| Strings, NaN, Infinity, booleans | Values are parsed as finite floats; booleans and non-finite values raise. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |
