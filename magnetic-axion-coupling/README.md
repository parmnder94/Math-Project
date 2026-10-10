# magnetic-axion-coupling

The linear magnetoelectric coefficient $\alpha_{ij}=\partial P_i/\partial B_j=\partial M_j/\partial E_i$ is the quantity that magnetoelectric experiments measure. In magnetic topological insulators it is the observable behind the "axion insulator" literature.

The model is the standard four-band lattice model of the Bi₂Se₃ family. It carries a staggered (Néel) mass, as in Li, Wang, Qi and Zhang, Nat. Phys. 6, 284 (2010), and an orbital-dependent exchange, as in Yu et al., Science 329, 61 (2010). The magnetic field enters through Peierls phases. The task asks for

$$\alpha_{zz} = 0.2013552\ldots\ e^2/h,$$

taken modulo $e^2/h$ in $(-1/2,1/2]$, and graded at an absolute error of $3\times10^{-3}$.

## Difficulty

**What makes it hard.**

1. **θ is not the answer.** The axion angle gives only the topological, isotropic part, $\theta e^2/2\pi h$, with $\theta/2\pi=0.2374$. The orbital magnetoelectric response of a band insulator also has a non-topological cross-gap (Kubo) part. The orbital dependence of the exchange makes it large here: $-0.036\,e^2/h$. The standard θ calculation (second Chern form along a gapped path, or a smooth-gauge Chern–Simons integral) is therefore off by 0.036, twelve times the bound.
2. **The full coefficient needs one of three routes:**
   * The complete multiband Kubo expression for the cross-gap term, plus θ. It is long and sign-sensitive.
   * The literal definition $\partial P_z/\partial B_z$ on magnetic (Hofstadter) supercells. This involves the Berry-phase polarization of $2q$ bands, the $e/q$ polarization quantum and its branch, magnetic-translation sampling, Wilson-loop discretization errors of order $10^{-3}$, and an extrapolation in the flux.
   * The converse route $\partial M_z/\partial E_z$, from the orbital magnetization of multiband slabs in a field, including the itinerant-circulation term.
3. **A tempting shortcut fails.** The top-surface (layer-resolved) anomalous Hall conductivity of a slab reproduces $\theta/2\pi$, not $\alpha_{zz}$ (`authoring/provenance/surface_hall.py`).
4. **Branch and orientation.** The result is defined modulo $e^2/h$, and the sign of the Peierls phase sets its sign.

**Shortcuts give wrong values** (`authoring/evidence/ablations.json`):

| shortcut | α_zz (e²/h) |
|---|---|
| axion angle only, θ/2π | 0.23737 |
| top-surface Hall conductivity of a slab | 0.23737 |
| half-quantized value | 0.5 |
| sign or field orientation reversed | −0.20136 |
| Hofstadter supercell, flux 1/32, forward difference | 0.1549 |

**Who does this in practice.** This is the work of a researcher in first-principles magnetoelectrics or topological band theory, the community that computes the Chern–Simons and Kubo contributions for candidate materials. For a specialist it is two to three days of derivation, coding and cross-checking.

## Reference solution

`solution/solve.py` (NumPy only) runs in about 15 s. It uses the Maxwell relation $\partial P_z/\partial B_z=\partial M_z/\partial E_z$ and works on the converse side:

1. **Slab in a field.** A slab of $L$ layers is placed in a uniform $E_z$, which adds $+eE_zz_l$ to the electron's on-site energy in layer $l$.
2. **Orbital magnetization.** $M_z=(e/\hbar)\,\mathrm{Im}\sum_{n\,\rm occ,\,m\,\rm unocc}(E_n+E_m)\langle n|\partial_xH|m\rangle\langle m|\partial_yH|n\rangle/(E_m-E_n)^2$, integrated over $k_x,k_y$. The slab Chern number is zero, so the chemical potential drops out.
3. **Slope.** $dM_z/dE_z$ per area is linear in $L$. The slope per layer, multiplied by $2\pi$, is $\alpha_{zz}$ in units of $e^2/h$. The slopes for $L=10\to14$, $14\to18$ and $18\to22$ are 0.2013557, 0.20135517 and 0.20135517.

## Verification

`tests/test_answer.py` reads `alpha_zz` from `/app/output/answer.json` and requires an absolute error of at most $3\times10^{-3}$ against `tests/truth.json`. The verifier never runs agent code.

**Why the bound is $3\times10^{-3}$** (1.5 % of α):

* **Wrong answers fail.** The most likely wrong answer, θ/2π (which the surface Hall conductivity also gives), is 0.036 away, twelve times the bound.
* **Correct answers pass.** Every converged correct method lands inside. The slab route agrees to $10^{-8}$. The supercell route agrees once its flux extrapolation is done; raw values at small $q$ are visibly unconverged.

**Ground truth: how it is produced.** The truth comes from the slab $\partial M_z/\partial E_z$ route (`authoring/provenance/slab_dMdE.py`, the same method as the reference solver), run at $80^2$ $k$-points and $L=18\to22$ layers.

**Independent checks** (`authoring/evidence/convergence.json`):

* **Literal definition.** `authoring/provenance/supercell_dPdB.py` computes $\partial P_z/\partial B_z$ on magnetic supercells and shares no code path with the truth. It uses the Landau gauge, the Berry-phase polarization of the $2q$ lowest bands and a central difference in the flux $\pm1/q$. With the Wilson-loop error removed and the series extrapolated in $1/q^2$, it agrees with the truth within its extrapolation uncertainty.
* **Normalization and sign.** In the Clifford limit $J_1=J_2=0$, where the cross-gap part vanishes, the slab result equals θ/2π from an independent Chern–Simons code to $4\times10^{-6}$ (`authoring/evidence/local_runs/clifford_limit_log.txt`).
* **Decomposition.** θ/2π = 0.2373567 from the Chern–Simons integral, and the same value from the surface Hall conductivity. So the cross-gap part is −0.0360.

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/provenance/slab_dMdE.py 80 14,18,22` reruns the truth series.
* `python3 authoring/provenance/supercell_dPdB.py 32,48,64 128` reruns the supercell check (about 1 h).
* `python3 authoring/provenance/clifford_limit.py`, `cs_stretched.py 64 0` and `surface_hall.py 24 64` rerun the decomposition checks.
* `python3 authoring/evidence/ablate.py <out.json>` recomputes the shortcut values.
* `bash authoring/evidence/verifier_local_run.sh <dir>` reruns the verifier image. Oracle and nop runs are in `authoring/evidence/local_runs/`.
