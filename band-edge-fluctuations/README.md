# band-edge-fluctuations

A spin-1/2 electron hops on a random "triangle network". The network is built from two independent uniform random permutations σ, τ of N sites, with bonds v→σ(v), v→τ(v) and v→σ(τ(v)). The hoppings are bond-dependent spin-orbit matrices T_γ = exp(iλ s_γ), and there is a Zeeman field B. The observable is the Lorentzian-broadened density of states just inside the top of the band, on the non-trivial subspace:

X_N = Tr[h(H_N) Q_N],  h(x) = η/((x − x₀)² + η²),  x₀ = 5.6, η = 0.06 (λ = 0.6, B = 1.5; the band top is at ≈ 5.93).

The task asks for four constants:

* z_∞, the large-N density: E X_N = 2N z_∞ + m + o(1);
* m, the O(1) mean correction;
* v = lim Var X_N;
* k = lim E(X_N − E X_N)³.

| constant | truth (z_inf, m, v to about 1e-9 relative; k to about 3e-8) |
|---|---|
| z_∞ | 0.1476936669067 |
| m | 0.1179457472883 |
| v | 35.834664243632 |
| k | 1.5662942790283 |

## Difficulty

**No dataset.** The model is fully specified, and the four targets are exact constants of it. The ground truth comes from the mathematics (see Verification), not from any solver output.

The solution needs four steps, each research-level:

1. **The limit object.** By Bordenave–Collins strong convergence, H_N becomes P = Bσ_z + T_aλ(a) + T_bλ(b) + T_cλ(a)λ(b) + h.c. in the regular representation of F_2. That is a tree of triangles with SU(2) hoppings. z_∞ needs its exact resolvent, from an operator-valued tree recursion, a linearization or a triangle cavity.
2. **An exact trace expansion.** For every N, X_N = (N−1) tr c_e + Σ_{w≠e} tr(c_w)(fix(w(σ,τ)) − 1). Here c_w = ⟨δ_w, h(P)δ_e⟩, a Green's-function coefficient of a group word.
3. **Random-word-map statistics.**
   * Write w = g u^d g⁻¹ with u primitive. Then E fix(w) → τ(d), the number of divisors of d (Nica).
   * The fluctuations form a compound-Poisson limit: independent Poisson(1/j) cycle counts for each primitive conjugacy class, with each class identified with its inverse class (Linial–Puder).
   * Hence the limits are the class sums v = Σ σ(gcd(n,n′)) Σ_κ T_κ(n)T_κ(n′) and k = Σ J₂(gcd) Σ_κ T_κ T_κ T_κ, where T_κ(n) sums tr c_w over all conjugates of u^{±n}.
4. **Exact resummation (the new obstacle).**
   * Just inside the band edge, the per-letter transfer operators have spectral radius close to 1. The radius is 0.94 for the mixed z/z̄ transfer that controls v. Class contributions therefore decay only like 0.94^ℓ.
   * Direct enumeration is hopeless. It is the route a frontier agent used to solve the earlier, faster-converging heat-trace version of this task in under 15 minutes (`authoring/evidence/probe_reports.md`). Even enumeration to class length 10 is still 81% relative off in v.
   * Every term of total power ≤ 4 has to be summed over all cyclically reduced words by transfer operators on tensor powers of the branch propagators:
     * Ω is inserted at each rotation's cut, using tr(ΩΦ_SΦ_P) = tr(Φ_P Ω Φ_S);
     * inverse rotations are handled through transposed inverse-letter slots;
     * cyclic closures are needed for powers;
     * the result is then restricted to primitive classes by Möbius inversion.
   * Total power 4 is not optional. Resumming only up to total power 3 leaves k off by 9.4e-5 relative.

None of steps 3–4 can be checked by simulation at the required precision. Monte Carlo of the finite networks (N = 200, 80,000 samples) pins v only to about 0.5%, m to about 20% and k not at all (standard error larger than k). Conceptual slips are invisible to it but break the 1e-6 gate (table below):

* dropping the j | d power-cycle structure;
* using E fix(u^d) = 2 for every proper power;
* using a Gaussian CLT.

**Probe history.** Each probe was an isolated frontier agent working from the instruction alone, with general internet access.

| version | what was asked | outcome |
|---|---|---|
| 1 | gap-opening field of the same model (cusp problem) | solved to 1e-15 in ~27 min |
| 2 | heat-trace statistics (fast convergence) | solved to ≤1e-9 in ~12 min |
| 3 (shipped) | this version | **solved**: reliable at ~36 min, final errors ≤ 3.2e-8 relative, by an independent method-of-images route |

**Honest difficulty assessment.** Only one attempt was sampled on the shipped version, and it succeeded. The 0–1/8 target is therefore **not** demonstrated. The task is still a research-grade derivation (four theory steps plus an exact resummation), and every shortcut route fails the gate (table below).

## Reference solution

`solution/solve.py` uses NumPy only and runs in about 18 minutes on 2 cores (1094 s in the environment container). Its largest deviation from the truth is 2.5e-7 relative (k), 4× inside the gate.

1. **4×4 self-adjoint linearization.** P = L₀ + Y*Y, with Y = 1·λ(a⁻¹) + T_c·λ(b) and Q = −1. Operator-valued branch resolvents H_s on the 4-regular tree are found by Newton's method at z and z̄. Then G(w,e) = Φ_{s₁}⋯Φ_{s_k}G(e,e), with Φ_s = H_sA_s.
2. **Conjugator sums.** F(y) = Σ_g tr c_{gyg⁻¹} = tr(Ω_(first,last) Φ_y) comes from one adjoint superoperator solve.
3. **Mean.** Σ_y F(y^d) for d ≤ 4 uses Kronecker-power transfers with cyclic closure.
4. **Variance and third cumulant.** Every class-sum product of total power ≤ 4 uses a marked transfer with inserted Ω's, combining z/z̄ slots and forward/inverse rotations, followed by Möbius inversion to primitive classes.
5. **Tails.** All terms of total power ≥ 5 converge like (0.17)^ℓ and are summed over primitive classes up to length 12.

## Verification

`tests/test_answer.py` parses `/app/output/answer.json`.

* Keys are normalised, and booleans, NaN and Infinity are rejected.
* Each of the four constants must be within a **relative error of 1e-6** of `tests/truth.json`.
* The verifier never runs agent code.

**Ground truth** comes from `authoring/provenance/truth.py`, with its report in `truth_report.json`.

* `tests/truth.json` is the exact-resummation assembly with primitive-class tails to length 13 (`authoring/evidence/truth_runs/lmax13.txt`, 1707 s).
* Going from length 11 to 13 changes m by 7.7e-9, v by 1.4e-10 and k by 2.9e-7, relative (`truth_runs/lmax11.txt`).
* **Independent derivation.** Probe 3 used a different route (cactus cavity plus method-of-images holonomy determinants, with no linearization). It agrees with the truth to 1e-13 (z_inf), 9.8e-10 (m), 3.7e-11 (v) and 3.2e-8 (k). The k difference is consistent with the truth's remaining tail, and it is 30× inside the gate.

**Further checks.**

* Every transfer construction (Kronecker squares and cubes, two- and three-slot marked transfers with powers, inverse rotations, z/z̄ mixes) reproduces brute-force enumeration *word length by word length* to ≤ 3e-13 (`authoring/provenance/validation/`, log in `authoring/evidence/validation_logs.txt`).
* The conjugator resummation was checked against explicit conjugator enumeration on independent cactus-geodesic Green's functions (no linearization), agreeing to 2e-13 (`authoring/evidence/superseded_heat_trace/truth_report.json`, part B).
* The word-map statistics (E fix → τ(d), Cov → σ(gcd)) were checked by simulation.
* Monte Carlo of the finite networks agrees with the predicted E X_N and Var X_N (`authoring/evidence/mc_summary.json`).

**Gate calibration** (`authoring/evidence/ablations.json`; worst relative error over the four keys, against the reference engine):

| route | worst relative error | key |
|---|---|---|
| Gaussian fluctuations (k = 0) | 1.0 | k |
| power cycles ignored (j = 1 only) | 3.7e-2 | k |
| E fix(u^d) = 2 for every proper power | 1.3e-2 | m |
| E X_N = 2(N−1)z + m convention | 2.5 | m |
| exact resummation only to total power 3 | 9.4e-5 | k |
| no resummation, classes ≤ 6 | 0.91 | v |
| no resummation, classes ≤ 8 | 0.87 | v |
| no resummation, classes ≤ 10 | 0.81 | v |

Every route misses the 1e-6 gate. The smallest miss is 9.4e-5, which is 94× the gate.

**Rebuild.**

* `python3 authoring/provenance/truth.py <out>` regenerates the truth (assemblies at lengths 11, 12 and 13, plus the 6×6 cross-check), in a few hours on 2 cores.
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image.
* Oracle and nop runs are in `authoring/evidence/local_runs/`.
