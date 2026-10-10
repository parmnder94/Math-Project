A magnetic topological insulator on a simple cubic lattice is described by the four-band tight-binding Bloch Hamiltonian

$$
H(\mathbf k)=M(\mathbf k)\,\tau_z+A_1\sin k_z\,\sigma_z\tau_x+A_2\left(\sin k_x\,\sigma_x+\sin k_y\,\sigma_y\right)\tau_x+m\,\tau_y+\sigma_z\left(J_1\frac{1+\tau_z}{2}+J_2\frac{1-\tau_z}{2}\right),
$$

$$
M(\mathbf k)=M_0-2B_1(1-\cos k_z)-2B_2(2-\cos k_x-\cos k_y),
$$

with $ \mathbf k $ in units of the inverse lattice constant. The Pauli matrices $ \sigma_i $ act on spin and $ \tau_i $ on the orbital index, $ \sigma_i\tau_j $ stands for $ \sigma_i\otimes\tau_j $, and a single Pauli matrix is tensored with the identity on the other factor. The parameters, in eV, are

$$
M_0=0.5,\quad A_1=0.6,\quad A_2=0.9,\quad B_1=0.3,\quad B_2=0.6,\quad m=0.25,\quad J_1=0.35,\quad J_2=-0.25 .
$$

The real-space model is fixed by $ H(\mathbf k)=\sum_{\mathbf R}t_{\mathbf R}\,e^{i\mathbf k\cdot\mathbf R} $, where $ t_{\mathbf R} $ is the $ 4\times4 $ block of matrix elements $ \langle\mathbf 0|\hat H|\mathbf R\rangle $ between the orbitals of the cell at the origin and those of the cell at lattice vector $ \mathbf R $, and all four orbitals sit at the lattice sites. The electrons, of charge $ -e $, fill the two lowest bands. A static uniform magnetic field $ \mathbf B=\nabla\times\mathbf A $ enters only through Peierls phases: every block $ \langle\mathbf R|\hat H|\mathbf R'\rangle $ is multiplied by $ \exp\!\left[-\frac{ie}{\hbar}\int_{\mathbf R'}^{\mathbf R}\mathbf A\cdot\mathrm d\mathbf l\right] $ along the straight segment, and there is no Zeeman coupling.

The magnetoelectric coefficient of this insulator is $ \alpha_{zz}=\partial P_z/\partial B_z $ at zero field, where $ \mathbf P $ is its electric polarization. It is defined modulo $ e^2/h $. What is $ \alpha_{zz} $ in units of $ e^2/h $, taken in the interval $ (-1/2,1/2] $?

Write `/app/output/answer.json` containing a JSON object with one numeric field, `alpha_zz`, holding $ \alpha_{zz} $ in units of $ e^2/h $. It is compared with the exact value and must agree to an absolute error of at most $ 3\times10^{-3} $.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
