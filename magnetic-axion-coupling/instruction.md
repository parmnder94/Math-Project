A magnetic topological insulator is described by the four-band tight-binding Bloch Hamiltonian

$$
H(\mathbf k)=M(\mathbf k)\,\tau_z+A_1\sin k_z\,\sigma_z\tau_x+A_2\left(\sin k_x\,\sigma_x+\sin k_y\,\sigma_y\right)\tau_x+m\,\tau_y+\sigma_z\left(J_1\frac{1+\tau_z}{2}+J_2\frac{1-\tau_z}{2}\right),
$$

$$
M(\mathbf k)=M_0-2B_1(1-\cos k_z)-2B_2(2-\cos k_x-\cos k_y),
$$

with $ \mathbf k=(k_x,k_y,k_z)\in[-\pi,\pi)^3 $ (all lattice constants equal to 1). The Pauli matrices $ \sigma_i $ act on spin and $ \tau_i $ on the orbital index, $ \sigma_i\tau_j $ stands for $ \sigma_i\otimes\tau_j $, and a single Pauli matrix is tensored with the identity on the other factor. The parameters, in eV, are

$$
M_0=0.28,\quad A_1=0.22,\quad A_2=0.40,\quad B_1=0.08,\quad B_2=0.50,\quad m=0.03,\quad J_1=0.20,\quad J_2=-0.05 .
$$

At half filling the two lowest bands are occupied at every $ \mathbf k $. Let $ |u_1(\mathbf k)\rangle,|u_2(\mathbf k)\rangle $ be orthonormal states spanning this occupied subspace, chosen smooth and periodic over the Brillouin zone, and let $ \mathcal A_j^{mn}(\mathbf k)=i\langle u_m(\mathbf k)|\partial_{k_j}u_n(\mathbf k)\rangle $. The axion angle is

$$
\theta=-\frac{1}{4\pi}\int_{[-\pi,\pi)^3}\mathrm d^3k\;\epsilon^{ijk}\operatorname{tr}\!\left[\mathcal A_i\,\partial_{k_j}\mathcal A_k-\frac{2i}{3}\mathcal A_i\mathcal A_j\mathcal A_k\right]\pmod{2\pi},
$$

with indices running over $ x,y,z $ and $ \epsilon^{xyz}=+1 $. What is $ \theta $, taken in the interval $ (-\pi,\pi] $?

Write `/app/output/answer.json` containing a JSON object with one numeric field, `theta`, holding $ \theta $ in radians. It is compared with the exact value and must agree to an absolute error of at most $ 10^{-3} $.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
