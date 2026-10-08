Fix $N$ and draw two independent, uniformly random permutations $\sigma$ and $\tau$ of $\{0, \dots, N-1\}$. Every vertex $v$ of the network has three outgoing bonds: an $a$-bond to $\sigma(v)$, a $b$-bond to $\tau(v)$ and a $c$-bond to $\sigma(\tau(v))$. So $v$, $\tau(v)$ and $\sigma(\tau(v))$ always span a triangle.

A spin-1/2 electron moves on this network. Its Hamiltonian on $\mathbb{C}^2 \otimes \mathbb{C}^N$ (spin $\otimes$ site) is

$$H_N = B\,(\sigma_z \otimes I_N) + \sum_{\gamma \in \{a, b, c\}} \left( T_\gamma \otimes P_\gamma + T_\gamma^\dagger \otimes P_\gamma^{T} \right),$$

where $P_a e_v = e_{\sigma(v)}$, $P_b e_v = e_{\tau(v)}$, $P_c = P_a P_b$, and $T_\gamma = \cos(\lambda)\, I_2 + i \sin(\lambda)\, s_\gamma$ with $(s_a, s_b, s_c) = (\sigma_x, \sigma_y, \sigma_z)$, the Pauli matrices. Use $\lambda = 0.6$ and $B = 1.5$.

The subspace $\mathbb{C}^2 \otimes \operatorname{span}\{(1, \dots, 1)\}$ is invariant under $H_N$. Let $Q_N$ be the orthogonal projection onto its orthogonal complement. Near the top of the band we look at the Lorentzian-broadened density of states

$$X_N = \operatorname{Tr}\left[ h(H_N)\, Q_N \right], \qquad h(x) = \frac{\eta}{(x - x_0)^2 + \eta^2}, \qquad x_0 = 5.6, \quad \eta = 0.06.$$

$X_N$ is random through $\sigma$ and $\tau$. As $N \to \infty$:

- $\mathbb{E}[X_N] = 2N z_\infty + m + o(1)$
- $\operatorname{Var}(X_N) \to v$
- $\mathbb{E}\left[(X_N - \mathbb{E}[X_N])^3\right] \to k$

These limits exist. We need all four constants $z_\infty$, $m$, $v$ and $k$.

Write `/app/output/answer.json` containing a JSON object with four numeric fields:

- `z_inf`: $z_\infty$
- `mean_correction`: $m$
- `variance`: $v$
- `third_cumulant`: $k$

Each value is compared with the exact constant and must agree to a relative error of at most $10^{-6}$. Extra keys are ignored.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
