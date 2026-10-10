# floquet-spin-interferometry

A spin-3/2 qudit is used as a holonomic gate built from a calibrated multiphoton transition. The |±3/2⟩ and |±1/2⟩ manifolds are two drive photons apart. An AC-Stark compensation (a_N) and a second-harmonic tone (b_N) cancel every lower-order shift and coupling, so the first surviving coupling is of eighth order. The amplitude scaling ε_N⁸ ΔT_N = 1/10 fixes that coupling's pulse area as the pulses become adiabatic, while slow frames steer the quantization axis. A control qubit compares the three-pulse gate with its literal time reverse, an echo diagnostic of time-reversal asymmetry in the non-Abelian holonomy. The task asks for the N → ∞ probability of measuring |−⟩ on the control:

$$p_- = 0.686758897365029385815907\ldots$$

The answer is graded at an absolute error of 10⁻⁹.

## Difficulty

**What makes it hard.** The limit is controlled by a coupling that first appears at eighth order of the high-frequency expansion. The derivation has to:

1. **Find the cluster.** Recognise that the states |q, −2⟩ (q ∈ Q) and |p, 0⟩ (p ∈ P) form a degenerate quasienergy cluster in Floquet space.
2. **Show the low orders cancel.** Orders one through seven vanish: odd orders by photon parity, and orders 2, 4 and 6 exactly through the calibration.
3. **Compute the eighth-order block exactly.** There are 60 ordered operator words, and the result is B₈ = (23539/2160)(Q−P) + (8447/540)X. One must also justify that no energy-dependent term enters at this order. No constant in the instruction encodes this block, and a 10⁻⁶ relative error in it already moves p by 8×10⁻⁸.
4. **Track the carrier phase.** It enters the P–Q coupling as X → X_φ.
5. **Project the slow connection.** Its secular part is −(3π/2)K_r, including the D_r term.
6. **Get the endpoint frames right.** S_r = diag(iN_{3α}, I): the Q sector rotates with azimuth 3α.
7. **Handle literal time reversal.** It sends φ to −φ and flips the connection.

**Shortcuts give wrong values** (`authoring/evidence/ablations.json`):

| shortcut | p |
|---|---|
| endpoint frame D_r(1) set to the identity | 0.474 |
| slow connection omitted | 0.493 |
| connection −πK_r only (D_r term dropped) | 0.749 |
| reversed pulse keeps the carrier phase | 0.745 |
| carrier phases ignored | 0.724 |
| (Q−P) part of the eighth-order block omitted | 0.981 |
| eighth-order coupling omitted | 1 |

Simulation cannot substitute for the derivation. Corrections decay like ε_N ∝ N^(−1/8): at N = 10¹⁶ the error is still about 10⁻³.

**Who does this in practice.** This is the work of a researcher in quantum control or Floquet engineering who uses high-frequency expansions, quasi-degenerate (Feshbach or van Vleck) perturbation theory and geometric phases. For a specialist it is several days of careful derivation.

## Reference solution

`solution/solve.py` uses SymPy for the exact recurrence and mpmath for 40-digit exponentials, and runs in about 3 s:

1. **Recurrence.** The exact energy-zero Feshbach recurrence confirms B₁…B₇ = 0 and returns B₈.
2. **Pulses.** U_r = S_r e^{−iG_r} with G_r = (1/10)[(23539/2160)(Q−P) + (8447/540)X_{φ_r}] − (3π/2)K_r.
3. **Reversed pulses.** Ũ_r = e^{−iG̃_r} S_r†, where G̃_r has φ → −φ and +(3π/2)K_r.
4. **Endpoint frames.** S_r = (−(7i/3)J_r + (4i/3)J_r³)(Q − iK_r), exactly.
5. **Result.** w = Tr(B†AP)/2 and p₋ = (1 − Re w)/2.

## Verification

`tests/test_answer.py` reads `p` from `/app/output/answer.json` and requires an absolute error of at most 10⁻⁹ against `tests/truth.json`, which holds 30 digits. The verifier never runs agent code.

**Ground truth: how it is produced.** The reference solver computes the truth from the derivation the task author supplied (following Eckardt & Anisimovas, New J. Phys. 17, 093039 (2015)). The truth and the reference therefore share one derivation.

**Independent checks** (`authoring/evidence/exact_checks.json`, `authoring/evidence/simulation_convergence.json`):

* **B₈.** A separate exact Floquet recurrence (`authoring/provenance/floquet.py`) reproduces B₁…B₇ = 0 and B₈.
* **Evaluation.** An independent script that computes the endpoint frames by numerical expm agrees to 30 digits. The same pipeline, given the constants of the earlier rational version of this task, reproduces its exact fraction to 3×10⁻⁴¹.
* **Direct simulation of the literal protocol** (`sim4.py` on `sim2.py` and `sim3.py`; `sim_brute.py` steps every period at small N):
  * The one-period propagator is computed in 40-digit arithmetic, and the slow terms in its interaction picture.
  * This reaches N = 10¹⁶. The deviation from the truth shrinks steadily with ε_N.
  * Polynomial extrapolation to ε → 0 brackets the truth to about 3×10⁻⁵.

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/evidence/ablate.py solution <out.json>` recomputes the shortcut values.
* `python3 authoring/provenance/sim4.py 1e10 1e16` re-runs the simulation, at about 80 s per point.
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image. Oracle and nop runs are in `authoring/evidence/local_runs/`.
