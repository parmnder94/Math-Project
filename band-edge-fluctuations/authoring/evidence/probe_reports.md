# Solver probes (isolated frontier-agent attempts during authoring)

Each probe ran as an autonomous agent in its own empty directory. It saw only the instruction, had general internet access, was told the 9000 s budget, and was forbidden to read authoring files.

## Probe 1: earlier draft (rejected as too easy)

* **Draft.** The same random triangle network and Hamiltonian, asking for the critical Zeeman field B_c at which the limiting new spectrum S(B) disconnects, and the energy E_c of gap opening (λ = 0.6, gate 1e-6 absolute).
* **Truth.** B_c = 3.4536676148355704, E_c = 0.6059655357419483, from a 40-digit triple-root solve of the cactus cavity equations.
* **Outcome.** Solved in about 1600 s, with B_c = 3.4536676148355707 and E_c = 0.6059655357419482.
* **Route.** Bordenave–Collins limit, then an independently derived triangle-cactus cavity system, then a Newton solve of the cusp conditions.
* **Consequence.** A spectral-edge question on this limit operator is within reach of current agents. The deliverable was moved to the finite-size statistics, which also need random-word-map theory.

## Probe 2: heat-trace version (rejected as too easy)

* **Draft.** The same network, asking for z_inf, m, v and k of X_N = Tr[exp(-0.3 H_N) Q_N] (gate 1e-5 relative).
* **Outcome.** The answer file written after about 740 s was already correct to 1e-9 relative or better: z_inf 6e-16, m 1.1e-9, v 1.5e-10, k 2e-16. The probe stopped on request at about 1110 s.
* **Route.** It wrote the exact trace expansion over words, used Nica's limits for word maps and the Poisson limits of Dumitriu–Johnson–Pal–Paquette, took cactus branch Green's functions on a contour, and summed classes up to cyclic length 6 and powers up to total length 8.
* **Why enumeration was enough.** For the heat kernel at beta = 0.3, class contributions fall about 8x per unit length, so direct enumeration reaches 1e-9.
* **Consequence.** The observable was moved to a Lorentzian-smoothed DOS just inside the upper band edge. There the class sums decay only like q^l with q > 0.9, and an exact resummation over conjugacy classes is required.

## Probe 3: final instruction

* **Draft.** The shipped instruction (Lorentzian DOS at x0 = 5.6, eta = 0.06; gate 1e-6 relative).
* **Outcome.** Solved. Values became reliable at about 2150 s; the probe was stopped on request at about 2590 s. Final answer versus `tests/truth.json` (relative error):

| key | probe 3 | truth | rel. error |
|---|---|---|---|
| z_inf | 0.1476936669067148 | 0.1476936669067 | 1e-13 |
| mean_correction | 0.1179457474038948 | 0.1179457472883 | 9.8e-10 |
| variance | 35.83466424496791 | 35.834664243632 | 3.7e-11 |
| third_cumulant | 1.5662943293476796 | 1.5662942790283 | 3.2e-8 |

* **Route (independent of the reference).** It used a triangle-cactus cavity with three 2x2 branch Green's functions instead of a linearization. For the finite-size terms it used a method-of-images formula, Delta_g(z) = -d/dz[log det(1-K_g) + log det(1-K_g^-1)], with 2x2 triangle-passage holonomies K_g, checked against brute-force unicyclic resolvents. It enumerated closed passage paths exactly to length 12 (about 22M paths). It summed the tails exactly with non-backtracking transfer operators on tensor powers, using truncation degree 4 for v, 5 for k and 7 for m.
* **What it found hardest.** The proper-power and 1/L, 1/(2L) weightings, and the long-cycle dominance of v near the edge.
* **Consequences.**
  * The task is solvable by a frontier agent within the budget. The target of 0-1 successes in 8 is **not** demonstrated: the one sampled attempt succeeded.
  * Because the probe's method is independent of the reference, its answer is also a second derivation of the ground truth. It agrees with `truth.json` to 1e-9 or better on z_inf, m and v. On k it agrees to 3.2e-8, which matches the remaining class-length truncation of the truth's k (the change from length 11 to 13 was 2.9e-7). This is 30x inside the gate.
