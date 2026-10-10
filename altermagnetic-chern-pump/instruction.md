A quantum-confined altermagnetic superconductor defines a gapped four-dimensional BdG pump on $ T^4 $, with $ x=(k_x,k_y,k_z,\phi) $ and

$$
(\ell_1,\ell_2,\ell_3,\ell_4)=(k_x+k_y,\ k_y+k_z,\ k_z+\phi,\ -k_x-k_y-k_z+\phi).
$$

Let $ d_i=\sin(3\ell_i) $ for $ i=1,\dots,4 $, $ p_i=\cos(3\ell_i) $, $ u_i=e^{2i\ell_i} $, and $ \omega=e^{2\pi i/3} $. Define

$$
\chi(p,q)=\frac{1+p+q-pq}{2},\qquad S(p,u)=\frac{1+p+(1-p)u}{2},
$$

$$
\Xi(u,v)=\frac{(1+v+v^2)+u(1+\omega^2v+\omega v^2)+u^2(1+\omega v+\omega^2v^2)}{3}.
$$

With cyclic $ p_5=p_1 $, set

$$
Z_1=\Xi(u_1,u_2)\,\Xi(u_3,u_4)\prod_{i=1}^4S(p_i,u_i),\qquad
Z_2=\Xi(u_1,u_3)\,\Xi(u_2,u_4)\prod_{i=1}^4S(p_{i+1},u_i),
$$

and

$$
d_5=\prod_{i=1}^4\chi(p_i,p_{i+1})\left[\frac{1}{5}+\operatorname{Re}Z_1+\frac{3}{5}\operatorname{Re}Z_2\right].
$$

For the spectrally flattened four-band BdG Hamiltonian $ H=\sum_{a=1}^5d_a\Gamma_a $, where the $ \Gamma_a $ mutually anticommute, define the occupied-band invariant by

$$
C_2=\deg\!\left(\frac{\mathbf d}{\lvert\mathbf d\rvert}:T^4\to S^4\right)
$$

with the standard orientations. Determine the single integer $ C_2 $.

No numerical sampling, CAS, symbolic root solver, or explicit enumeration of the $ 6^4 $, $ 2^4 $, or $ 3^4 $ residue sectors is allowed. Prove that the spectrum is everywhere gapped and evaluate the oriented covering number exactly. Do not assume that the mod-2 and mod-3 residue contributions factorize independently.

Write `/app/output/answer.json` containing a JSON object with one integer field, `C2`, holding $ C_2 $. It is compared with the exact value.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
