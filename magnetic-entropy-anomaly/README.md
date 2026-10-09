# magnetic-entropy-anomaly

A charged Brownian particle moves on a torus with a temperature profile $T = 1 + \tfrac35\cos x$, a drag tensor $\Gamma(y)$ that rotates as $y$ changes, and a fixed perpendicular magnetic field. The task asks for the exact gap

$$\Delta = \lim_{m\to0}\dot S_m - \mathcal I_0 = \frac{17}{420},$$

between the heat entropy flow into the environment in the small-mass limit and the time-reversal relative entropy rate of the limiting position process, with the field held fixed.

## Difficulty

**What makes it hard.** As the mass goes to zero, the position settles into a simpler overdamped diffusion. One would expect the heat released to equal how irreversible that limiting motion is. With a temperature gradient this is false: velocity fluctuations that the overdamped motion cannot see keep dissipating heat (the "entropy anomaly"). Computing the gap exactly takes three steps, each easy to get wrong:

1. **The overdamped limit.** The rotating drag tensor and the magnetic field produce a noise-induced drift. The stationary density turns out to be uniform, but the stationary state still carries a circulating current, $\mathbf j = \rho_0 (J/3)\nabla T$.
2. **The irreversibility of that limit.** Because of the hidden current, the fixed-field relative entropy rate is not zero: $\mathcal I_0 = 1/20$.
3. **The heat flow as $m\to0$.** Its natural expression blows up like $1/\sqrt m$. Removing the singularity needs an exact rewriting and a Poisson equation for the fast velocity process in the third Hermite sector, where the magnetic field couples four cubic modes. The result is $\lim\dot S_m = 19/210$.

**Plausible shortcuts give wrong exact values** (`authoring/evidence/ablations.json`):

| shortcut | value obtained |
|---|---|
| uniform density taken as reversible (or the field flipped in the reversal) | 19/210 |
| no anomaly: heat rate equals the overdamped irreversibility | 0 |
| magnetic field dropped from the fast velocity process | 19/400 |
| drag replaced by its isotropic average $\tfrac32 I$ | 724/16575 |

Simulation cannot replace the derivation. The answer is an exact fraction, and Monte Carlo of the full dynamics reaches only about 1% even at small mass.

**Who does this in practice.** This is the work of a researcher in stochastic thermodynamics or applied probability who studies small-mass limits of Langevin equations. For such a specialist it is about a day of careful derivation.

## Reference solution

`solution/solve.py` repeats the derivation in exact SymPy arithmetic and runs in about a second:

1. **Fast velocity process.** With $\mathbf u = \sqrt m\,\mathbf V$, the velocity at fixed position is an Ornstein–Uhlenbeck process with drift matrix $\Gamma - J$ and Maxwellian of variance $T$.
2. **Limiting diffusion.** The mobility is $M = (\Gamma - J)^{-1} = [R\,\mathrm{diag}(2,1)R^{\mathsf T} + J]/3$. The density is uniform and the current is $\rho_0 (J/3)\nabla T$. Then $\mathcal I_0 = \int \mathbf j^{\mathsf T}D^{-1}\mathbf j/\rho_0 = I_T/4$, with $I_T = \langle T'^2/T\rangle = 1/5$.
3. **Heat rate.** The exact identity $\dot S_m = -\langle S\rangle/(2\sqrt m)$ holds with $S = (\nabla T\cdot\mathbf u)(|\mathbf u|^2 - 4T)/T^2$. Solving $\mathcal L_0\chi = S$ in the cubic Hermite sector of the local drag frame (a 4×4 system coupled by the field) gives $\lim\dot S_m = -\langle\chi S\rangle/4 = (19/42) I_T = 19/210$.

## Verification

`tests/test_answer.py` reads `delta` from `/app/output/answer.json`. It accepts an exact fraction equal to $17/420$, or a decimal within $10^{-12}$ of it. The verifier never runs agent code.

**Ground truth: how it is produced.** The value was derived analytically. The author's step-by-step solution follows the small-mass framework of Birrell, *Entropy Anomaly in Langevin–Kramers Dynamics with a Temperature Gradient, Matrix Drag, and Magnetic Field*, J. Stat. Phys. (2018), doi:10.1007/s10955-018-2162-2. `solution/solve.py` recomputes it in exact arithmetic, so the truth and the reference share one derivation.

**Independent checks:**

* **Hermite operator.** The 4×4 matrix of $-\mathcal L_0$ on the cubic sector was re-extracted numerically from the generator at random frame angles and temperatures (residual about 1e-14). Every inverse, contraction and integral was rechecked in SymPy.
* **Simulation.** A Monte Carlo simulation of the full underdamped dynamics agrees with the predicted heat-rate limit $19/210 \approx 0.0905$. At $m = 0.04$, refining the time step and extrapolating it to zero gives 0.0905. At a fixed step the values for $m = 0.01$–$0.04$ coincide, so the remaining mass dependence is below the statistical error of about 0.001 (`authoring/evidence/mc_summary.json`, script `authoring/provenance/mc.py`).

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/evidence/ablate.py solution <out.json>` recomputes the shortcut values.
* `python3 authoring/provenance/mc.py <m> <particles> <time> <seed>` re-runs the simulation.
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image; oracle and nop runs are in `authoring/evidence/local_runs/`.
