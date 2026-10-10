# altermagnetic-chern-pump

A gapped four-band BdG Hamiltonian $H=\sum_{a=1}^5 d_a\Gamma_a$ on the four-torus $(k_x,k_y,k_z,\phi)$ models a quantum-confined altermagnetic superconductor used as a four-dimensional pump. The first four components are $d_i=\sin(3\ell_i)$ in sheared coordinates $\ell=Lx$. The mass term $d_5$ is built from a pairing prefactor $\prod\chi(p_i,p_{i+1})$ and from two altermagnetic form factors $Z_1$ and $Z_2$, which couple period-2 and period-3 structure. The occupied-band invariant is the second Chern number, which here equals the degree of $\mathbf d/|\mathbf d|:T^4\to S^4$. The task asks for

$$C_2=-74.$$

The answer is graded by exact integer match.

## Difficulty

**What makes it hard.** The integer is fixed by a signed count of where the d-vector points to one pole of $S^4$. To get it exactly, the derivation has to:

1. **Prove the gap.** At the zeros of $d_1,\dots,d_4$, i.e. $\ell_i=\pi n_i/3$, rewrite $d_5$ in Chinese-remainder variables $b_i=n_i \bmod 2$ and $a_i=n_i \bmod 3$. There $\Xi(\omega^a,\omega^b)=\omega^{ab}$, $S=u^b$, $Z_1=\omega^A$, $Z_2=\omega^B$, and $\prod\chi=(-1)^{Q_2(b)}$. Then show that the bracket takes only the values $9/5$, $9/10$, $3/10$ and $-3/5$.
2. **Normalise the count.** Use the north pole as the regular value, with local chirality $(-1)^{\sum b_i}$. The signed north-minus-south sum is twice the degree.
3. **Evaluate coupled character sums.** Express the degree through $W_0=-8$, $W_A=-219$, $W_B=-255$ and the joint count $W_{AB}=-113$. The joint count does not factorise; it needs nine four-variable quadratic Gauss sums over $\mathbb F_3$ twisted by binary characters.
4. **Restore the covering.** $x\mapsto\ell$ has $\det L=2$, so it is an orientation-preserving double cover: $C_2=2D_\ell$.

**Shortcuts give wrong values** (`authoring/evidence/ablations.json`):

| shortcut | C2 |
|---|---|
| covering factor $\det L=2$ omitted | −37 |
| signed north-minus-south sum not halved | −148 |
| mod-2 and mod-3 data treated as independent | −1106/9 (not an integer) |
| pairing prefactor $\prod\chi$ dropped | 78 |
| $Z_2$ term dropped | 210 |
| $Z_1$ term dropped | 138 |
| orientation reversed | +74 |

Numerical quadrature of the degree density does not settle the integer: on $24^4$ to $40^4$ grids it returns −109, −143 and −119 (`authoring/evidence/quadrature_attempts.json`).

**Who does this in practice.** This is the work of a researcher in topological condensed matter or mathematical physics: someone who classifies topological superconductors and higher-dimensional pumps and computes second Chern numbers through mapping degrees and finite-field (Gauss-sum) counting. For a specialist it is most of a working day of careful derivation.

## Reference solution

`solution/solve.py` uses exact rational arithmetic in $\mathbb Q(\omega)$ (standard library only) and runs in about 1 s:

1. **Character-sum route.** It computes $W_0,W_A,W_B$ from the completed-square level-set counts $Z_A,Z_B\in\{24,33\}$, and $W_{AB}$ from the characters $\omega^{rA+sB}$. Then $D_\ell=(-81W_0+2W_A+2W_B-2W_{AB})/2=-37$.
2. **Internal cross-check.** It evaluates $d_5$ exactly at all zeros of $d_1,\dots,d_4$, straight from the instruction's definitions. This confirms the gap ($\min|d_5|=3/10$) and the north-pole count $-37$.
3. **Covering.** $\det L=2$ is computed exactly, giving $C_2=-74$.

## Verification

`tests/test_answer.py` reads `C2` from `/app/output/answer.json` and requires exact equality with `tests/truth.json`. Integer-valued strings and floats are accepted; non-integers and booleans are rejected. The verifier never runs agent code.

**Ground truth: how it is produced.** The value −74 is the result of the task author's derivation (golden solution; framework of Qi, Hughes and Zhang, Phys. Rev. B 78, 195424 (2008), https://doi.org/10.1103/PhysRevB.78.195424). The reference solver reproduces it exactly.

**Independent checks** (`authoring/evidence/exact_checks.json`):

* **Intermediates.** A separate floating-point script (`authoring/provenance/residue_sums.py`) reproduces the author's $W_0,W_A,W_B,W_{AB}$, all nine Gauss sums $U_{rs}$, $D_\ell=-37$ and $\det L=2$.
* **Other regular values.** `authoring/provenance/regular_values.py` counts the preimages of $\pm e_1,\dots,\pm e_4$ (a different reduction: three $d_i$ vanish, then one-dimensional root finding of $d_5$, about 2,300 roots each, with minimum Jacobian $0.07$). All eight give degree −37 on the $\ell$-torus, independently of the north-pole residue argument.

**Rebuild.**

* `python3 solution/solve.py` recomputes the answer.
* `python3 authoring/evidence/ablate.py <out.json>` recomputes the shortcut values.
* `python3 authoring/provenance/regular_values.py` re-runs the regular-value check (about 4 min; needs NumPy and SciPy).
* `bash authoring/evidence/verifier_local_run.sh <dir>` re-runs the verifier image. Oracle and nop runs are in `authoring/evidence/local_runs/`.
