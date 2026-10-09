Fix $N$ and draw two independent, uniformly random permutations $\sigma$ and $\tau$ of $\{0, \dots, N-1\}$. They define a random 4-regular graph on the sites $0, \dots, N-1$: site $v$ has an $a$-bond to $\sigma(v)$ and a $b$-bond to $\tau(v)$. A spin-1/2 particle hops on this graph with the Hamiltonian on $\mathbb{C}^2 \otimes \mathbb{C}^N$ (spin $\otimes$ site)

$$H_N = B\,(\sigma_z \otimes I_N) + \sum_{\gamma \in \{a, b\}} \left( T_\gamma \otimes P_\gamma + T_\gamma^\dagger \otimes P_\gamma^{T} \right),$$

where $P_a e_v = e_{\sigma(v)}$, $P_b e_v = e_{\tau(v)}$, $T_a = \cos(\lambda)\, I_2 + i \sin(\lambda)\, \sigma_x$ and $T_b = \cos(\lambda)\, I_2 + i \sin(\lambda)\, \sigma_y$, with the Pauli matrices $\sigma_x, \sigma_y, \sigma_z$. Use $\lambda = 0.6$ and $B = 1.5$.

The subspace $\mathbb{C}^2 \otimes \operatorname{span}\{(1, \dots, 1)\}$ is invariant under $H_N$. Let $Q_N$ be the orthogonal projection onto its orthogonal complement, and consider the partition function

$$X_N = \operatorname{Tr}\left[ e^{-\beta H_N}\, Q_N \right], \qquad \beta = 0.3.$$

$X_N$ is random through $\sigma$ and $\tau$. As $N \to \infty$ its mean and variance have expansions in powers of $1/N$:

$$\mathbb{E}[X_N] = 2N z_\infty + m + \frac{m_1}{N} + O(N^{-2}), \qquad \operatorname{Var}(X_N) = v + \frac{v_1}{N} + O(N^{-2}).$$

We need the five constants $z_\infty$, $m$, $m_1$, $v$ and $v_1$.

Write `/app/output/answer.json` containing a JSON object with five numeric fields: `z_inf`, `m`, `m1`, `v` and `v1`. Each value is compared with the exact constant and must agree to a relative error of at most $10^{-6}$. Extra keys are ignored.

You have 3600 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
