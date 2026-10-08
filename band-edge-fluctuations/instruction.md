Fix N and draw two independent, uniformly random permutations σ and τ of {0, ..., N-1}. Every vertex v of the network has three outgoing bonds: an a-bond to σ(v), a b-bond to τ(v) and a c-bond to σ(τ(v)). So v, τ(v) and σ(τ(v)) always span a triangle.

A spin-1/2 electron moves on this network. Its Hamiltonian on C^2 ⊗ C^N (spin ⊗ site) is

H_N = B (σ_z ⊗ I_N) + Σ_{γ = a, b, c} ( T_γ ⊗ P_γ + T_γ^† ⊗ P_γ^T )

where P_a e_v = e_{σ(v)}, P_b e_v = e_{τ(v)}, P_c = P_a P_b, and T_γ = cos(λ) I_2 + i sin(λ) s_γ with (s_a, s_b, s_c) = (σ_x, σ_y, σ_z), the Pauli matrices. Use λ = 0.6 and B = 1.5.

The subspace C^2 ⊗ span{(1, ..., 1)} is invariant under H_N. Let Q_N be the orthogonal projection onto its orthogonal complement. Near the top of the band we look at the Lorentzian-broadened density of states

X_N = Tr[ h(H_N) Q_N ],   h(x) = η / ((x - x_0)^2 + η^2),   x_0 = 5.6,   η = 0.06.

X_N is random through σ and τ. As N → ∞:

- E[X_N] = 2N z_∞ + m + o(1)
- Var(X_N) → v
- E[(X_N - E[X_N])^3] → k

These limits exist. We need all four constants z_∞, m, v and k.

**Output.** Write `/app/output/answer.json` containing a JSON object with four numeric fields:

- `z_inf`: z_∞
- `mean_correction`: m
- `variance`: v
- `third_cumulant`: k

Each value is compared with the exact constant and must agree to a relative error of at most 1e-6. Extra keys are ignored.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
