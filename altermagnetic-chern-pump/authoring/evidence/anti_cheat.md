# Anti-cheat notes (never executed by the pipeline)

What the agent sees: the instruction only. The environment image is a bare Python 3.12 image with no extra packages. What stays sealed in `tests/`: `truth.json`, the integer C2.

| Attempt | Why it fails |
|---|---|
| Nop (no file) | The `answer` fixture cannot open `/app/output/answer.json`; reward 0. |
| Degree on the l-torus only (covering x -> l ignored) | Gives -37. The map x -> l has det L = 2. |
| Signed north-minus-south mass sum used as the degree, then doubled for the covering | Gives -148: the signed sum is already twice the north-pole degree. |
| Mod-2 and mod-3 data treated as independent (Z_AB = Z_A Z_B / 81) | Gives -1106/9, not even an integer. |
| Pairing prefactor prod chi, or one of the Z terms, dropped | Gives 78, 210 or 138 (`ablations.json`). |
| Orientation of S^4 or of the torus reversed | Gives +74. |
| Numerical quadrature of the degree density | The integrand is sharply peaked and highly oscillatory; 24^4 to 40^4 grids on the x-torus give -109, -143, -119 (`quadrature_attempts.json`). |
| Looking up the reference | The cited framework (Qi, Hughes, Zhang 2008) defines second Chern numbers of four-dimensional models in general; this d-vector and its number do not appear in the literature. |
| Strings, floats, booleans | The value must be an integer: "-74" and -74.0 are accepted, -74.5 and booleans are rejected. |
| Writing into `/logs/verifier` | The verifier runs in a separate container and reads only the declared artifact. |

Note on one coincidence: the signed north-minus-south sum on the l-torus without halving and without the covering factor is also -74, because the two omitted factors of 2 cancel. Reaching the right integer that way needs two compensating mistakes.

The instruction asks for an exact argument without explicit enumeration of the residue sectors. The verifier checks only the integer and cannot see the method.
