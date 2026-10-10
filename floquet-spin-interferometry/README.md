# floquet-spin-interferometry

A spin-3/2 qudit is used as a holonomic gate built from a calibrated multiphoton transition. The |±3/2⟩ and |±1/2⟩ manifolds are two drive photons apart. An AC-Stark compensation (a_N), a second-harmonic tone (b_N) and a 1/T_N gap correction cancel every lower-order shift and coupling, so the first surviving coupling is of eighth order. The amplitude scales as ε_N ∝ T_N^(−1/8), which keeps that coupling's pulse area fixed as the pulses become adiabatic, while slow frames steer the quantization axis. A control qubit compares the three-pulse gate with its literal time reverse, an echo diagnostic of time-reversal asymmetry in the non-Abelian holonomy. The task asks for the exact N → ∞ probability of measuring |−⟩ on the control:

$$p_- = \frac{1467167323933}{5157470703125} \approx 0.284474194500808 .$$

## Difficulty

**What makes it hard.** The limit is controlled by a coupling that first appears at eighth order in the high-frequency expansion. Getting it takes several steps:

1. **Floquet space.** Recognise that the states |q, −2⟩ (q ∈ Q) and |p, 0⟩ (p ∈ P) form a degenerate quasienergy cluster.
2. **Cancellation through seventh order.** Verify that the counterterms cancel every order through seven. Odd orders vanish by photon parity; orders 2, 4 and 6 cancel exactly.
3. **Eighth order.** Compute the full eighth-order block, 60 ordered operator words: B₈ = (23539/2160)(Q−P) + (8447/540)X. Justify that no energy-dependent term enters at this order.
4. **Carrier phase.** Track how the carrier phase enters the P–Q coupling (X → X_φ).
5. **Slow-frame connection.** Project the connection secularly: −(3π/2)K_r, including the D_r term.
6. **Endpoint frames.** Get the endpoint frames S_r = diag(iN_{3α}, I) right, including the Q-sector azimuth 3α.
7. **Reversed pulses.** Show that time reversal sends φ to −φ and flips the connection.
8. **Exact algebra.** Do the final matrix algebra exactly. The effective generators have spectrum {±2, ±1/2}.

**Shortcuts give wrong values** (`authoring/evidence/ablations.json`):

| shortcut | value |
|---|---|
| endpoint frame D_r(1) set to the identity | 0.769298… |
| slow connection omitted | 0 |
| connection −πK_r only (D_r term dropped) | 0.905920… |
| reversed pulse keeps the carrier phase | 0.904794… |
| carrier phases ignored | 0.905000… |

Simulation cannot substitute for the derivation. Corrections decay like ε_N ∝ N^(−1/8). At N = 10¹⁶ the error is still about 10⁻³, while the answer is an exact fraction.

**Who does this in practice.** This is the work of a researcher in quantum control or Floquet engineering who uses high-frequency expansions, quasi-degenerate (Feshbach/van Vleck) perturbation theory and geometric phases. For a specialist it is several days of careful derivation.

## Reference solution

`solution/solve.py` uses SymPy and exact arithmetic and runs in about 4 s:

1. **Floquet recurrence.** The exact energy-zero Feshbach recurrence W_n = R Σ_o V_o W_{n−o}, B_n = Π Σ_o V_o W_{n−o} confirms B₁…B₇ = 0 and returns B₈.
2. **Effective generator.** With Δ T_N ε_N⁸ = 1080π/(8447√3), the (Q−P) part of B₈ plus the gap shift is a pure scalar, and the X_φ coefficient is 2/√3. So G_r = π[−(γ/2)I + (2/√3)X_{φ_r} − (3/2)K_r].
3. **Pulses.** U_r = S_r e^{−iπM_r} and Ũ_r = e^{−iπM̃_r} S_r†, where:
   * exp(−iπM) = (4M² − I)/15 + (8i/15)(M³ − 4M), with the polynomial identity checked;
   * R_r(1) = −(7i/3)J_r + (4i/3)J_r³;
   * D_r(1) = Q − iK_r.
4. **Result.** w = Tr(B†AP)/2 = (2223136055259 + 2895894695712i)/5157470703125, so p₋ = (1 − Re w)/2.

## Verification

`tests/test_answer.py` reads `p` from `/app/output/answer.json`. It accepts an exact fraction equal to the truth, or a decimal within 10⁻¹². The verifier never runs agent code.

**Ground truth: how it is produced.** The value comes from the derivation the task author supplied (golden solution, following Eckardt & Anisimovas, New J. Phys. 17, 093039 (2015)). `solution/solve.py` recomputes it in exact arithmetic, so the truth and the reference share one derivation.

**Independent checks:**

* **Eighth-order coefficient.** A separate implementation of the Floquet recurrence (`authoring/provenance/floquet.py`) reproduces B₁…B₇ = 0 and B₈ exactly.
* **Final algebra.** A separate exact computation reproduces w and p₋. The polynomial exponential matches numerical expm to 10⁻¹². The fraction is reduced: its denominator is 5¹⁵·13² and the gcd is 1.
* **Direct simulation of the literal protocol** (`authoring/evidence/simulation_convergence.json`, scripts `sim.py`, `sim2.py`, `sim3.py`):
  * The one-period propagator is computed in 40-digit arithmetic, and the slow terms in its interaction picture.
  * This reaches N = 10¹⁶ (ε = 0.0066). Brute-force stepping agrees at N = 10³ and 10⁴.
  * The deviation from the exact value is linear in ε at small ε. Polynomial extrapolation to ε → 0 gives 0.284474 ± 2×10⁻⁵; the best fit is within 2×10⁻⁸ of 0.2844741945.

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/evidence/ablate.py solution <out.json>` recomputes the shortcut values.
* `python3 authoring/provenance/sim3.py 1e9 1e12 1e16` re-runs the simulation, at about 80 s per point.
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image. Oracle and nop runs are in `authoring/evidence/local_runs/`.
