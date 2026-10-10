"""Literal-protocol simulation for the hint-free variant (Delta = 1): eps_N^8 T_N = 1/10, no gap shift,
alpha = (0, pi/3, 2pi/3), phi = (pi/6, pi/4, pi/3).  Uses sim2.pulse and the 40-digit Floquet propagator of sim3."""
import numpy as np, sys
import sim2, sim3
sim2.gam = 0.0
sim2.alphas = [0.0, np.pi/3, 2*np.pi/3]; sim2.phis = [np.pi/6, np.pi/4, np.pi/3]
def prob(N):
    T = 2*np.pi*N; eps = (1/(10*T))**0.125
    A = np.eye(4, dtype=complex); B = np.eye(4, dtype=complex)
    for r in range(3): A = sim2.pulse(N, sim2.alphas[r], sim2.phis[r], False, sim3.floquet_mp(eps, sim2.phis[r])) @ A
    for r in (2, 1, 0): B = sim2.pulse(N, sim2.alphas[r], sim2.phis[r], True, sim3.floquet_mp(eps, -sim2.phis[r])) @ B
    w = np.trace(B.conj().T @ A @ sim2.P)/2
    return (1 - w.real)/2, eps
if __name__ == "__main__":
    pstar = 0.6867588973650293858
    for N in [int(float(v)) for v in sys.argv[1:]]:
        p, eps = prob(N); print(f"N={N:.0e} eps={eps:.5f} p={p:.10f} p-p*={p-pstar:+.3e}", flush=True)
