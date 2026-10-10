"""Literal protocol at large N without stepping every period.
Moving frame: H_phi = H_F(t) + (1/T)[-gamma pi Q + C(s)], H_F 2pi-periodic.  With U_F = exp(-i Lam), the product
prod_k U_F exp(-i <G(s_k)>_F) equals U_F^N V(1), dV/ds = -i N e^{i Lam N s} <G(s)>_F e^{-i Lam N s} V  (error O(1/N))."""
import numpy as np, scipy.linalg as sla, sys
from scipy.integrate import solve_ivp
import mpmath as mp
s3 = np.sqrt(3)
Jx = np.array([[0, s3/2, 0, 0], [s3/2, 0, 1, 0], [0, 1, 0, s3/2], [0, 0, s3/2, 0]], complex)
Jy = np.array([[0, -1j*s3/2, 0, 0], [1j*s3/2, 0, -1j, 0], [0, 1j, 0, -1j*s3/2], [0, 0, 1j*s3/2, 0]])
P = np.diag([0, 1, 1, 0]).astype(complex); Q = np.eye(4) - P
X = Jx @ Jx - (7*P + 3*Q)/4
alphas = [0, np.pi/2, np.arctan(3/4)]; phis = [np.arctan(1/2), np.arctan(2/3), np.arctan(3/4)]
gam = 23539/(8447*s3)
basis = [np.eye(4)[:, [i]] @ np.eye(4)[[j], :] for i in range(4) for j in range(4)]
def floquet(eps, ph):
    a = eps**2 + 4/3*eps**4 + 148/45*eps**6; b = eps**2 + 5/3*eps**4 + 167/36*eps**6
    HF = lambda t: 2*Q + 2*eps*np.cos(t + ph)*Jx + a*(P - Q) - 2*b*np.cos(2*(t + ph))*X
    def rhs(t, y):
        U = y[:16].reshape(4, 4); dU = -1j*HF(t) @ U
        return np.concatenate([dU.ravel()] + [(U.conj().T @ Bm @ U).ravel() for Bm in basis])
    sol = solve_ivp(rhs, (0, 2*np.pi), np.concatenate([np.eye(4).ravel(), np.zeros(256, complex)]), method='DOP853', rtol=1e-13, atol=1e-15)
    UF = sol.y[:16, -1].reshape(4, 4); L = sol.y[16:, -1].reshape(16, 4, 4)
    # enforce unitarity, Lam = i log UF (eigenphases near 0)
    u, _, vh = np.linalg.svd(UF); UF = u @ vh
    w, V = np.linalg.eig(UF); Lam = V @ np.diag(1j*np.log(w)) @ np.linalg.inv(V)
    return Lam, L
def pulse(N, alpha, phi, reverse, Fl):
    T = 2*np.pi*N
    Lam, L = Fl
    avg = lambda Y: np.einsum('k,kab->ab', Y.ravel(), L)
    Jr = Jx*np.cos(alpha) + Jy*np.sin(alpha); K = P @ Jr @ P
    D = lambda s: sla.expm(-1j*np.pi*s*K/2)
    def G(s):
        if reverse: Dm = D(1 - s); C = np.pi*Dm.conj().T @ Jr @ Dm + np.pi/2*K
        else: Dm = D(s); C = -np.pi*Dm.conj().T @ Jr @ Dm - np.pi/2*K
        return avg(-gam*np.pi*Q + C)/T
    NL = N*Lam
    def rhs(s, y):
        E = sla.expm(1j*NL*s)
        Ht = N*(E @ G(s) @ E.conj().T)
        return (-1j*Ht @ y.reshape(4, 4)).ravel()
    sol = solve_ivp(rhs, (0, 1), np.eye(4, dtype=complex).ravel(), method='DOP853', rtol=1e-11, atol=1e-13)
    U = sla.expm(-1j*NL) @ sol.y[:, -1].reshape(4, 4)
    S = sla.expm(-1j*np.pi*Jr) @ D(1)
    return S @ U if not reverse else U @ S.conj().T
def prob(N):
    T = 2*np.pi*N; eps = (1080*np.pi/(8447*s3*T))**0.125
    A = np.eye(4, dtype=complex); B = np.eye(4, dtype=complex)
    for r in range(3): A = pulse(N, alphas[r], phis[r], False, floquet(eps, phis[r])) @ A
    for r in (2, 1, 0): B = pulse(N, alphas[r], phis[r], True, floquet(eps, -phis[r])) @ B
    w = np.trace(B.conj().T @ A @ P)/2
    return (1 - w.real)/2, eps
if __name__ == "__main__":
    for N in [int(float(v)) for v in sys.argv[1:]]:
        p, eps = prob(N); print(f"N={N:.0e} eps={eps:.5f} p={p:.10f} p-p*={p-1467167323933/5157470703125:+.3e}", flush=True)
