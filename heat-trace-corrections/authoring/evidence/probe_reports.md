# Solver probes (isolated agent attempts during authoring)

Each probe ran as an autonomous agent in its own empty directory. It saw only the instruction, had general internet access, was told the 9000 s budget, and was forbidden to read authoring files. A "lighter" probe used a smaller model than the "frontier" probe. These records document how the task was calibrated. They are not part of the difficulty argument, which rests on the mathematics alone (README).

## Earlier versions (rejected as too easy)

| version | what was asked | outcome |
|---|---|---|
| 1 | gap-opening field of a spin-orbit triangle network (cusp of the limit spectrum) | frontier: solved in ~27 min |
| 2 | O(1) heat-trace statistics of the triangle network | frontier: solved in ~12 min |
| 3 | O(1) band-edge Lorentzian statistics (z, m, v, k), exact transfer resummation required | frontier: solved in ~36 min; the platform's lighter screening pass also solved it too often |
| 4 | the same on square-plaquette networks (c-bond to στσ⁻¹, 6×6 linearization required) | lighter: solved in ~31 min; frontier: solved in ~45 min |

Each of these needed known theorems: strong convergence, linearization, Nica and Linial–Puder, and transfer operators. That made them reachable by applying the literature. The shipped version asks for the 1/N coefficients, for which no off-the-shelf formula exists.

## Shipped version (heat-trace corrections)

**Frontier probe.** Solved in about 59 min, by an independent organisation of the theory:

* class sums of the heat kernel to word length 15, from a compiled depth-first pass;
* rank-1 and rank-2 "image graphs" in place of my quotient walk;
* exact finite-N moments for N ≤ 9 as a check.

Final answer against `tests/truth.json` (relative error):

| key | probe | truth | rel. error |
|---|---|---|---|
| z_inf | 1.3078808147368293 | 1.3078808147368284 | 7e-16 |
| m | −2.447016342390509 | −2.447016342390545 | 1.5e-14 |
| m1 | 0.011810119934744512 | 0.011810119934665837 | 6.7e-12 |
| v | 2.884123369340993 | 2.884123369340829 | 5.7e-14 |
| v1 | 0.2786938384300517 | 0.2786938384186302 | 4.1e-11 |

The probe's derivation is independent of the reference, so its answer is also a second derivation of the ground truth.

**Lighter probe.** Solved in about 68 min, also independently:

* an Ihara–Bass-type factorization, Tr e^{−βH_G} = 2N z_∞ + Σ_w φ(w) fix_w, checked on small graphs;
* rank-1 and rank-2 quotient counts for the mean, summed to word length 12 in C++;
* a core-count decomposition of the variance, validated on toy weightings against exact S_N character theory.

Its relative errors against the truth are 3.4e-11 (m1) and 6.0e-11 (v1); z_inf, m and v agree to 1e-13 or better.

**Small-model probe (Haiku class).** Not solved. It stopped after about 32 min with only z_inf correct:

* m was right to about 1e-5;
* m1 was rough (0.01192 against 0.011810);
* v came only from Monte Carlo (2.889 ± 1%);
* v1 was missing.

Its naive Poisson-cycle variance formula failed its own checks, and it did not find the quotient expansion for the 1/N terms.

**Calibration consequence.** Both model tiers solved this version within the 9000 s budget, in about an hour each. The agent timeout in `task.toml` is the remaining calibration lever. The shipped budget is the 3600 s minimum: the frontier probe needed about 58 min for its first correct answer file, the lighter probe about 68 min, and the small model gave up. Every probe was told 9000 s, so these times were not taken under time pressure.
