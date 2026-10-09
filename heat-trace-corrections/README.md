# heat-trace-corrections

A spin-1/2 particle hops on a random 4-regular graph built from two independent uniform random permutations σ, τ of N sites (bonds v→σ(v) and v→τ(v)). The hoppings are spin-orbit matrices T_a = exp(iλσ_x) and T_b = exp(iλσ_y), and there is a Zeeman field B. The observable is the heat trace on the non-trivial subspace:

X_N = Tr[exp(−βH_N) Q_N],  λ = 0.6, B = 1.5, β = 0.3.

The task asks for the first two orders of the large-N expansions of its mean and variance:

E X_N = 2N z_∞ + m + m₁/N + O(N⁻²),  Var X_N = v + v₁/N + O(N⁻²).

| constant | truth (`tests/truth.json`) |
|---|---|
| z_∞ | 1.3078808147368284 |
| m | −2.447016342390545 |
| m₁ | 0.011810119934665837 |
| v | 2.884123369340829 |
| v₁ | 0.2786938384186302 |

## Difficulty

**No dataset.** The model is fully specified, and the five targets are exact constants of it. The ground truth is computed by the reference machinery (`solution/solve.py`) run at larger cutoffs. Its correctness rests on independent checks, listed under Verification.

**Leading orders (known theory).**

* By strong convergence, H_N becomes P = Bσ_z + T_aλ(a) + T_bλ(b) + h.c. on F_2.
* For every N, X_N = (N−1) tr c_e + Σ_{w≠e} tr(c_w)(fix(w(σ,τ)) − 1), where c_w = ⟨δ_w, e^{−βP}δ_e⟩ is the heat kernel on the free group.
* Nica's limit E fix(u^d) → τ(d) and the Linial–Puder Poisson cycle limit give m and v.

**The 1/N coefficients (the core of the task).** No off-the-shelf formula gives m₁ and v₁. Both need the exact expansion of word-map moments, E fix_w(N) and E[fix_{w₁} fix_{w₂}](N), to order 1/N, summed against the class sums of the heat kernel over every class and every pair of classes.

* **The Linial–Puder quotient expansion.** E[Π_i fix_{w_i}](N) = Σ_Γ (N)_V/((N)_{E_a}(N)_{E_b}). The sum runs over the folded quotients Γ of the disjoint union of the cycle graphs C_{w_i}, with V vertices and E_a, E_b edges. Each term is N^{χ(Γ)}(1 + O(1/N)).
* **Single words.** For w = u^d, E fix_w = τ(d) + a₁(w)/N + O(N⁻²) with
  a₁(w) = #{quotients with χ = −1} − Σ_{χ = 0 quotients} E_a E_b.
  The second term comes from expanding the falling factorials of the cycle quotients themselves. It is easy to miss, and the critical-subgroup count alone does not cover proper powers.
* **Pairs of words.** Cov(fix_{w₁}, fix_{w₂}) = C₀ + C₁/N + O(N⁻²), where C₀ = #{connected χ = 0 quotients} (= σ(gcd) for related words, else 0). C₁ has three parts:
  * the connected χ = −1 quotients of C_{w₁} ⊔ C_{w₂};
  * −E_aE_b on each connected χ = 0 quotient;
  * a cross term −Σ(E_a(Γ₁)E_b(Γ₂) + E_b(Γ₁)E_a(Γ₂)) over pairs of χ = 0 quotients of the two words separately. It arises because (N)_{p+q} ≠ (N)_p(N)_q at order 1/N for disconnected quotients. Without it v₁ is off by 180%.
* **Enumeration and truncation.** The quotients must be enumerated exactly (every partition once) for all classes up to length 11 and all pairs of classes up to total length 12, proper powers included. The pair sums alternate in sign and fall off slowly, so pairs that stop at total length 10 already miss the 1e-6 gate.

None of this can be checked by simulation at the required precision. Monte Carlo of the finite graphs cannot resolve m₁: its effect m₁/N is below 1e-3 of the standard deviation of X_N at N ≥ 10. It resolves v₁ only to tens of percent. Exact averages over S_N × S_N exist only for N ≤ 7, where the 1/N series cannot be extrapolated to 1e-6.

**Who does this in practice.** This is the work of a research mathematician in random matrix theory and combinatorial group theory: someone who works on word measures on symmetric groups, Stallings graphs and finite-size corrections for random graph covers. For a specialist it is several days of derivation and validation.

## Reference solution

`solution/solve.py` uses NumPy only and runs in about 80 s on 2 cores (87 s in the environment container). Its largest deviation from the truth is 1.0e-8 relative (v₁).

1. **Class sums.** On the 4-regular Cayley tree the resolvent factorizes, G_z(w,e) = Φ_{s₁}⋯Φ_{s_k}G_z(e,e), with branch resolvents found by Newton's method. The conjugator sums are resummed exactly, F_z(y) = tr(Ω Φ_y). The heat kernel is the trapezoid contour integral (1/2πi)∮ e^{−βz} F_z dz on |z| = 12 with 48 nodes, which is exponentially convergent.
2. **Quotients.** A depth-first walk enumerates the folded quotients of one cycle, or of the union of two cycles, with χ ≥ −1. Classes are labelled in order of first appearance, so each partition appears exactly once. This gives a₁, C₀ and C₁. The counts are cached on orbits of the eight signed letter permutations.
3. **Assembly.** Primitive classes up to length 11 with their powers, and pairs of classes up to total length 12:
   * z_∞ = tr c_e/2 and m = Σ F (τ(d) − 1) − tr c_e;
   * m₁ = Σ F a₁;
   * v = Σ F F C₀ and v₁ = Σ F F C₁.

## Verification

`tests/test_answer.py` parses `/app/output/answer.json`.

* Keys are normalised, and booleans, NaN and Infinity are rejected.
* Each of the five constants must be within a **relative error of 1e-6** of `tests/truth.json`.
* The verifier never runs agent code.

**Ground truth: how it is produced.** `tests/truth.json` is **not** independent of the reference code. `authoring/provenance/truth.py` runs the reference solver `solution/solve.py` itself, with larger cutoffs: classes up to length 13 and pairs up to total length 14, against 11 and 12 in the shipped solution. The convergence record is in `authoring/evidence/truth_report.json`. Relative to the truth, the settings below give:

| setting | m₁ | v₁ |
|---|---|---|
| classes 12, pairs 13 | 3.4e-11 | 4.3e-9 |
| classes 11, pairs 12 (reference solution) | 2.4e-9 | 1.0e-8 |

**Why the truth can be trusted: independent checks.** These use separate code paths (`authoring/provenance/checks.py`, results in `authoring/evidence/checks.json`), plus two independent derivations:

* **Heat kernel.** The contour trace equals the exact moment series Σ (−β)ⁿ τ(Pⁿ)/n! (walks on reduced words) to 9e-16.
* **Word maps.** The quotient expansion of E fix and E[fix fix] equals exact averages over all of S₅ × S₅ for 16 single words and pairs, as exact rationals.
* **Finite N.** E X_N and Var X_N from the class sums plus the finite-N expansion match brute force over S_N × S_N for N = 4…7. The brute force takes σ over conjugacy-class representatives, weighted. The mean agrees to ≤ 6e-12 and the variance to ≤ 4e-8 (pair truncation).
* **Coefficient extraction.** m₁ and v₁ from the pruned χ ≥ −1 enumeration match Richardson extrapolation of the full (unpruned) expansion. That extrapolation is evaluated in exact rational arithmetic at N = 2·10⁴ … 1.6·10⁵, with the same class truncation, and agrees to 3e-10 (m₁) and 3e-9 (v₁).
* **Independent derivations.** Two solvers derived the theory separately (`authoring/evidence/probe_reports.md`). One used a compiled enumeration with rank-1 and rank-2 image graphs, the other an Ihara–Bass factorization with core counts. Both reproduce all five constants, agreeing with `tests/truth.json` to about 1e-11 or better.
* **Consistency with exact small-N data.** N(Var X_N − v) = 0.345, 0.350, 0.333, 0.324 for N = 4…7. This approaches v₁ = 0.279 with a 1/N² term, and the exact means behave the same way for m₁.

**Gate calibration** (`authoring/evidence/ablations.json`; worst relative error over the five keys, against the truth):

| route | worst relative error | key |
|---|---|---|
| 1/N of E fix from χ = −1 quotients only (no −E_aE_b) | 40 | m₁ |
| 1/N of Cov without the disconnected cross term | 1.8 | v₁ |
| 1/N of Cov without −E_aE_b on connected χ = 0 quotients | 0.64 | v₁ |
| leading order only (m₁ = v₁ = 0) | 1.0 | m₁, v₁ |
| proper powers given no 1/N term | 5.3e-6 | m₁ |
| classes and pairs truncated at length 10 | 1.5e-6 | v₁ |

Every route misses the 1e-6 gate.

**Rebuild.**

* `python3 authoring/provenance/truth.py <dir>` regenerates the truth and its convergence record.
* `python3 authoring/provenance/checks.py <file>` re-runs the independent checks.
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image.
* Oracle and nop runs are in `authoring/evidence/local_runs/`.
