# magnetic-entropy-anomaly

A charged Brownian particle moves on a torus with a rotating drag tensor $\Gamma(y)$ and a uniform magnetic field. The medium has a **temperature step**: $T = 2$ on the stripe $0 < x < \pi$ and $T = 1$ on $\pi < x < 2\pi$. The heat entropy flow $\dot S_m$ diverges like $m^{-1/2}$ as the mass goes to zero. The task asks for

$$K = \lim_{m\to0}\sqrt m\,\dot S_m = 0.08339296 .$$

The answer is graded at a relative error of $10^{-4}$.

## Difficulty

**What makes it hard.** For a smooth temperature profile the small-mass limit of the entropy flow is a standard calculation (Hilbert expansion, local Maxwellian, cubic Hermite sector). A temperature step breaks it:

* The local anomaly formula (density proportional to $|\nabla T|^2/T$) diverges.
* The heat is exchanged in kinetic boundary layers of width $\sqrt m$ along the two lines where $T$ jumps.

Getting $K$ takes four steps:

1. **Bulk.** Show that the density is uniform in each stripe. Show that the exact momentum balance of the layer makes the pressure $\rho T$ continuous across the step, so $\rho T = P = 1/(3\pi^2)$.
2. **Reduce the heat to the layers.** With no work done on the particle, the heat released at $x$ is $(\Gamma{:}\Pi - T\,\mathrm{tr}\Gamma\,\rho)/m$. An energy balance in the layer turns each line's contribution into the energy flux $q = \int u_1|u|^2 p\,du$ at the step. Each line contributes $\int dy\,(q/2)(1/T_{\rm right} - 1/T_{\rm left})$.
3. **Solve the layer.** Solve the half-space (transmission) kinetic problem with a 2D velocity, the matrix drag in its rotated frame and the Lorentz force. Do this at every angle of the drag frame, then average along the lines.
4. **Precision.** The layer solution is singular at zero normal velocity, so velocity expansions converge only algebraically. Reaching $10^{-4}$ needs a controlled extrapolation.

**Shortcuts give wrong values or none** (`authoring/evidence/ablations.json`):

| shortcut | value |
|---|---|
| smooth-profile anomaly formula, or overdamped heat count | no finite value |
| one-dimensional layer (transverse velocity ignored) | 0.05414 (35% off) |
| drag frame frozen at $y = 0$ | 0.08592 (3% off) |
| Euler–Maruyama simulation at small mass | about 1% high |

**Who does this in practice.** This is the work of a researcher in kinetic theory or stochastic thermodynamics who handles boundary layers of Fokker–Planck equations (Milne-type problems). It is several days of analysis and careful numerics.

## Reference solution

`solution/solve.py` uses NumPy and SciPy and runs in about a minute on 2 cores. It does the following:

* **Layer equation.** It writes the layer equation $u_1\partial_\xi p = \nabla_u\cdot[(\Gamma - J)u\,p + T(\xi)\Gamma\nabla_u p]$ in a tensor Hermite basis at reference temperature $3/2$.
* **Bounded solutions.** It takes the bounded solutions on each side from ordered real Schur forms (stable invariant subspaces) plus the Maxwellian.
* **Matching.** It imposes continuity at the step.
* **Energy flux.** It reads off $q$.

The layer flux has symmetries: $q$ at $x = \pi$ equals $-q$ at $x = 0$. So $K = -\langle q_a\rangle_\theta/(3\pi)$, averaged over the frame angle with 8 points.

## Verification

`tests/test_answer.py` reads `K` from `/app/output/answer.json`. It requires a relative error of at most $10^{-4}$ and never runs agent code.

**Ground truth: how it is produced.** `tests/truth.json` comes from the reference solver itself, run at larger normal Hermite orders. `authoring/provenance/truth.py` runs K1 = 120 and 160 with Richardson extrapolation at the measured exponent 2.87. The convergence record is in `authoring/evidence/layer_convergence.json`.

* At one angle, K1 = 40 to 400 converges like $K_1^{-2.9}$.
* The 120/160 extrapolation agrees with the 240/320/400 sequence to 6e-8.
* The transverse order is converged to 1e-14, and 8 angles are exact to 1e-10.

**Independent checks** (`authoring/evidence/mc_summary.json`):

* **Consistency of the layer solution.** The matching residual with $\rho T = 1$ imposed on both sides is about 1e-15. The two symmetries hold to 1e-14.
* **Finite-mass solution (1D model).** The same step without field or drag rotation was solved exactly at finite mass on the circle by Hermite modes. It gives $\sqrt m\,\dot S_m$ = 0.054438, 0.054290 and 0.054216 at m = 0.01, 0.0025 and 0.000625. That converges like $\sqrt m$ to the layer constant 0.054142.
* **Simulation.** An exact Ornstein–Uhlenbeck-step Monte Carlo gives the following:
  * In 1D at m = 0.0025: 0.05438 ± 0.00036, against 0.05429.
  * In 2D at m = 0.0025: 0.08389 ± 0.00083, against 0.08339 plus a positive $O(\sqrt m)$ correction.

  An Euler–Maruyama scheme is about 1% high at every mass because it inflates the velocity variance.

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/provenance/truth.py solution` regenerates the truth.
* The scripts in `authoring/provenance/` re-run the checks: `circle1d.py`, `mc1d_ou.py`, `mc2d_ou.py` and `layer1d.py`.
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image. Oracle and nop runs are in `authoring/evidence/local_runs/`.
