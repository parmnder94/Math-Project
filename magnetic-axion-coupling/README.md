# magnetic-axion-coupling

The axion angle $\theta$ sets the topological magnetoelectric effect of an insulator. An electric field induces a magnetization, and a magnetic field a polarization, both equal to $\theta e^2/(2\pi h)$ times the field. Time reversal or inversion pins $\theta$ to $0$ or $\pi$. In magnetic topological insulators where both are broken, together with their product, $\theta$ is no longer quantized, and its value controls the magnetoelectric and dynamical-axion responses that experiments look for.

The model is the standard four-band lattice model of the Bi₂Se₃ family. It carries a staggered (Néel) mass $m\,\tau_y$, which breaks time reversal and inversion, as in Li, Wang, Qi and Zhang, Nat. Phys. 6, 284 (2010). It also carries an orbital-dependent exchange from a uniform magnetization along $z$, which breaks their product, as in Yu et al., Science 329, 61 (2010). The task asks for the Chern–Simons axion angle of the two occupied bands, taken in $(-\pi,\pi]$:

$$\theta = 1.885246645979\ldots$$

It is graded at an absolute error of $10^{-3}$.

## Difficulty

**What makes it hard.** The difficulty is conceptual, not numerical precision.

1. **No formula to plug into.** The orbital-dependent exchange takes the Hamiltonian out of the Clifford ($\sum_a d_a\Gamma_a$) form, so the closed-form axion formula of Li et al. does not apply. $\theta$ has to come from the non-Abelian Berry connection of two occupied bands.
2. **A gauge problem.** The Chern–Simons integral is defined only in a gauge that is smooth and periodic over the whole Brillouin zone, while numerical eigenvectors carry random phases. There are two routes:
   * Build such a gauge: check that the occupied bundle is trivial, find trial states whose projection never becomes singular, and orthonormalize.
   * Integrate the gauge-invariant second Chern form along a path of Hamiltonians that stays gapped and starts at an insulator with known $\theta$, which requires its own $\mathbb Z_2$ analysis.
3. **Natural choices fail.**
   * Projecting onto the orbital basis states is singular, because the bands are inverted near the zone center: the smallest overlap eigenvalue is $10^{-5}$ on a $64^3$ grid.
   * A path that switches on the exchange before the Néel mass closes the gap and passes through a Weyl semimetal (`authoring/evidence/gauge_and_path_checks.json`).
4. **Sign and branch.** $\theta$ is defined modulo $2\pi$, and its sign depends on the orientation and on the Berry-connection convention.
5. **A small-gap region.** The parameters sit near the field-driven transition into a Weyl semimetal, where the axion response is largest. The direct gap is only 0.050 eV, on a ring around $(0,0,\pi)$. The calculation must be resolved and checked for convergence there; an unconverged $64^3$ grid is off by $2\times10^{-3}$.

**Physically wrong shortcuts** (`authoring/evidence/ablations.json`) all miss by at least 0.46:

| shortcut | θ |
|---|---|
| quantized value assumed | 3.14159 |
| Néel mass only (the closed Clifford formula) | 2.34523 |
| sign convention reversed | −1.88525 |
| orbital dependence of the exchange dropped | gap closes, θ undefined |

**Who does this in practice.** This is the work of a researcher in topological band theory or first-principles magnetoelectrics who predicts the magnetoelectric response of candidate magnetic topological insulators. For a specialist it is about one to two days of derivation, coding and checking.

## Reference solution

`solution/solve.py` (NumPy only) runs in about 35 s with 1.3 GB of memory:

1. **Gauge.** It projects the trial states $|s\rangle\otimes|\tau_y=-1\rangle$ onto the occupied subspace and Löwdin-orthonormalizes them, $U=P\Phi(\Phi^\dagger P\Phi)^{-1/2}$. The overlap stays above 0.43 everywhere, so $U$ is smooth and periodic.
2. **Grid.** Each coordinate is mapped by $k=k_0+(s-k_0)-c\sin(s-k_0)$, with $c=0.8$ and $k_0=(0,0,\pi)$. This smooth periodic map preserves orientation and concentrates points at the small-gap ring. The Chern–Simons density is a 3-form, so it is evaluated directly in $s$ with FFT derivatives on a uniform $128^3$ grid. A plain uniform grid works as well, only more slowly.
3. **Result.** $\theta=-\frac1{4\pi}\int\epsilon^{ijk}\operatorname{tr}[\mathcal A_i\partial_j\mathcal A_k-\tfrac{2i}3\mathcal A_i\mathcal A_j\mathcal A_k]$, wrapped into $(-\pi,\pi]$.

## Verification

`tests/test_answer.py` reads `theta` from `/app/output/answer.json` and requires an absolute error of at most $10^{-3}$ against `tests/truth.json`. The verifier never runs agent code.

**Why the bound is $10^{-3}$.**

* **Physically sufficient.** It is 0.05 % of $\theta$, finer than any magnetoelectric measurement resolves.
* **Every converged correct method passes.** Both gauge routes, with uniform or stretched grids from $72^3$ upward and the gauge-free path from 10 quadrature nodes upward, land within $1.3\times10^{-4}$. Only an unconverged $64^3$ uniform grid misses ($1.9\times10^{-3}$), and doubling the grid exposes it. See the table in `authoring/evidence/convergence.json`.
* **Every physically wrong answer fails.** The nearest one, the Néel mass alone, is 0.46 away, about 460 times the bound.

**Ground truth: how it is produced.** The truth comes from the author's Chern–Simons code (`authoring/provenance/cs_stretched.py`, the same method as the reference solver), run at $160^3$ points with two stretch strengths. The two agree to $5\times10^{-12}$. Framework: Qi, Hughes and Zhang, Phys. Rev. B 78, 195424 (2008); Essin, Moore and Vanderbilt, Phys. Rev. Lett. 102, 146805 (2009).

**Independent checks**, all well inside the bound:

* **Gauge-free route.** `authoring/provenance/path_second_chern.py` shares no code or gauge with the truth. It integrates the second Chern form, with $F_{ab}=iP[\partial_aP,\partial_bP]P$, along a gapped path (minimum gap 0.050 eV) from the time-reversal-symmetric strong topological insulator ($\theta=\pi$), switching on first the Néel mass and then the exchange. It lands $1.8\times10^{-5}$ from the truth with 10 Gauss nodes per segment on a $48^3$ grid, $3\times10^{-7}$ with 14 nodes on $96^3$, and $1.6\times10^{-9}$ with 28 nodes on $96^3$.
* **Literature formula.** At the Néel-only point, where the model has the Clifford form, the closed formula of Li et al. reproduces the code's $|\theta|=2.34523$ to $9\times10^{-7}$ (`authoring/evidence/literature_check.json`). The formula's overall sign follows its own Γ-matrix convention; the task's sign follows the definition written in the instruction.
* **Uniform grids.** Results agree with the truth to $10^{-4}$ from $72^3$ upward.
* **Sensitivity check, not independent.** A second trial gauge, an entangled spin–orbital trial subspace, changes the result by $8\times10^{-11}$.

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/provenance/cs_stretched.py 96,128,160 0.8` reruns the convergence series; use `0` as the stretch for uniform grids.
* `python3 authoring/provenance/path_second_chern.py 96 28 0.8` reruns the gauge-free check (about 40 min).
* `python3 authoring/provenance/clifford_formula_check.py` reruns the literature check.
* `python3 authoring/provenance/gap_scan.py` recomputes the gap data.
* `python3 authoring/evidence/ablate.py <out.json>` recomputes the shortcut values.
* `bash authoring/evidence/verifier_local_run.sh <dir>` reruns the verifier image. Oracle and nop runs are in `authoring/evidence/local_runs/`.
