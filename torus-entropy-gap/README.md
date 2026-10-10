# torus-entropy-gap

A Brownian particle on the torus (x, y) is coupled to two heat reservoirs acting along a frame that rotates with the angle θ = x + y + ½cos x. The reservoirs have frictions (1, 3) and temperatures (T, 2T), with T = 1 + (4/5)cos x. A magnetic term ½JV and an explicit trigonometric force F also act. The task asks for the exact value of

Δ = Ṡ₀ − σₓ,

where Ṡ₀ is the finite part of the reservoir entropy flow as the mass m → 0 (Ṡ_m = κ/m + Ṡ₀ + o(1)), and σₓ is the time-reversal relative-entropy rate of the path of the single coordinate x of the overdamped limit.

$$\Delta = -\frac{814956740459}{15500795001120} \approx -0.0525751576$$

Only the exact fraction is accepted.

## Origin

This is a hardened version of a problem supplied by the task author, whose original statement and step-by-step derivation give −141231719/146966400. The reference pipeline here reproduces that value exactly with the original parameters; see Verification. Compared with the original:

* the force is given only as explicit trigonometric expressions, so the fast covariance and the structure F = div Σ + Σ∇(log T²) + Bh are no longer printed;
* the 1/(9m) subtraction is replaced by "the constant κ for which the limit exists";
* phrases that pointed to the method ("complete continuous record", "holding J fixed") and the literature reference are removed;
* the parameters are new: the density is ∝ T² rather than ∝ T, the reservoir temperatures are not proportional to the frictions (so the entropy weights are anisotropic), and the magnetic strength, frictions, amplitude and deformation all change.

## Difficulty

**Why it is hard.**

1. **Fast covariance.** The scaled velocity w = √m V is locally Gaussian with the non-Maxwellian covariance Σ = T R C₀ Rᵀ, from the Lyapunov equation B₀C₀ + C₀B₀ᵀ = 2Q₀ with B₀ = diag(1,3) − J/2 and Q₀ = diag(1,6). This gives C₀ = [[55,−6],[−6,103]]/52.
2. **Hidden structure of the force.** F = div Σ + Σ∇(log T²) + Bh, with h = (1/(3T²), T′/(2T)). Nothing in the statement names this.
3. **Overdamped limit.** The limit L_eff = −⟨A L⁻¹ A⟩ has drift MF + (∂ₗM)Σₗ (with a noise-induced part) and diffusion D = T R D₀ Rᵀ, where M = B⁻¹ and D₀ = [[168,−24],[−24,100]]/169. Its stationary density is T²/⟨T²⟩, with a nonzero current j = ρ v and v = (1/(3T²), T′/(2T) + (21/26)T′).
4. **Divergent and finite heat.** In Itô form, Ṡ_m = (1/m)E[wᵀWw/T − 4], with W = R diag(1, 3/2) Rᵀ.
   * A quadratic Poisson observable Φ = wᵀZw gives κ = 2 tr(Q₀Z₀) − 4 = 3/104, and the exact identity Ṡ_m − κ/m = m^(−1/2) E[AΦ].
   * The finite part needs the first kinetic correction to the stationary density: Ṡ₀ = −∫ f₀ A L⁻¹ A Φ = 3283241117173/1291732916760. The current terms do not vanish.
5. **Partial observation.** The x path alone is not Markov.
   * Its quadratic variation gives 2T a(θ), and its covariation with a(θ_t) gives a′(θ). Together these fix θ mod π.
   * The symmetry y → y + π makes the remaining twofold lift carry no relative entropy.
   * So σₓ equals the full-position rate ∫ ρ vᵀD⁻¹v = 1012560749/390300768.

**Shortcuts give wrong values** (`authoring/evidence/ablations.json`):

| shortcut | Δ |
|---|---|
| kinetic current ρh used in the relative entropy | +0.4415 |
| x treated as a 1D Markov diffusion with y-averaged coefficients | +0.7396 |
| density taken ∝ T instead of ∝ T² | −1.1285 |
| current part Bh of the force dropped in the finite heat | −1.9194 |
| both reservoirs weighted by 1/T (T₂ = 2T ignored) | +0.6202 |
| observed path taken as reversible (σₓ = 0) | +2.5417 |
| **true value** | **−0.0526** |

**Who does this in practice.** A researcher in stochastic thermodynamics or multiscale stochastic analysis who studies entropy production of coarse-grained and partially observed Langevin systems. For a specialist it is several days of derivation and exact computation.

## Reference solution

`solution/solve.py` uses SymPy only, runs in about 1 s and works entirely in exact rational arithmetic:

1. **Representation.** Fields are polynomials in cos x, sin x, cos θ, sin θ and 1/T, differentiated by the chain rule.
2. **Fast operators.** L⁻¹ is solved in the reservoir frame, where L has constant coefficients. Gaussian averages use Isserlis' theorem.
3. **Limiting generator.** −⟨A L⁻¹ A⟩ is applied to the coordinates. It is checked against the closed-form drift, and the stationarity of ρ = T²/⟨T²⟩ is checked exactly (⟨(div j)²⟩ = 0).
4. **Exact averages.** The θ-average is taken monomial by monomial (θ is uniform in y at fixed x). Odd powers of sin x drop out, and the rest becomes a Laurent polynomial in T whose moments ⟨Tⁿ⟩ are exact (Legendre polynomials, √(1−a²) = 3/5).
5. **Force check.** The explicit force of the instruction is checked against the structured form to 30 digits.

## Verification

`tests/test_answer.py` reads `Delta` from `/app/output/answer.json` as a fraction string and requires exact equality with `tests/truth.json`. Decimals are rejected. The verifier never runs agent code.

**Ground truth: how it is produced.** The reference solver computes it from first principles: the homogenized generator, the Hilbert expansion and the Girsanov rate. Steps 3 and 5 of the Difficulty list are justified analytically, in the solver docstring and in the solution explanation.

**Independent checks:**

* **Original problem.** The task author derived −141231719/146966400 for the original parameters by a different route: explicit cubic-moment tensors in the reservoir frame, rotation averages and density-weighted integrals. The same solver reproduces that value exactly, and also both of its parts, 41249/28800 and 8657741/18370800 (`ablations.json`, last entry).
* **Exact internal checks.** Stationarity of ρ, the closed-form drift, and the equality of the explicit and structured force.
* **Monte Carlo of the literal SDE.** `authoring/provenance/mc_heat2.py` simulates the underdamped dynamics at finite mass (exact frozen-coefficient Gaussian steps) and estimates Ṡ_m − κ/m. The results are in `authoring/evidence/monte_carlo.json`.

**Rebuild.**

* `python3 solution/solve.py`
* `python3 authoring/evidence/ablate.py solution <out.json>`
* `python3 authoring/provenance/mc_heat2.py <m> <dt/m> <particles> <tau> <seed>`

Reference for the small-mass limit: Hottovy, McDaniel, Volpe and Wehr, Commun. Math. Phys. 336, 1259–1283 (2015).
