"""Direct simulation of the literal protocol (Delta = 1) at finite N.
Moving frame W_r(s) = R_r(s) D_r(s): H_phi(t) = H_F(t) + (1/T)[-gamma*pi*Q + C(s)], H_F exactly 2pi-periodic,
C(s) = -pi D^+ J_r D - (pi/2) K_r (forward) or +pi D(1-s)^+ J_r D(1-s) + (pi/2) K_r (reversed).
One period: U_F (DOP853, rtol 1e-13) times exp(-i <C(s_k)>_F / T), <Y>_F = int_0^{2pi} U_F(t)^+ Y U_F(t) dt
(second order in 1/T per period, total error O(1/N))."""
import numpy as np, scipy.linalg as sla, sys
from scipy.integrate import solve_ivp
s3 = np.sqrt(3)
Jx = np.array([[0, s3/2, 0, 0], [s3/2, 0, 1, 0], [0, 1, 0, s3/2], [0, 0, s3/2, 0]], complex)
Jy = np.array([[0, -1j*s3/2, 0, 0], [1j*s3/2, 0, -1j, 0], [0, 1j, 0, -1j*s3/2], [0, 0, 1j*s3/2, 0]])
P = np.diag([0, 1, 1, 0]).astype(complex); Q = np.eye(4) - P
X = Jx @ Jx - (7*P + 3*Q)/4
alphas = [0, np.pi/2, np.arctan(3/4)]; phis = [np.arctan(1/2), np.arctan(2/3), np.arctan(3/4)]
gam = 23539/(8447*s3)
def pulse(N, alpha, phi, reverse):
    T = 2*np.pi*N; eps = (1080*np.pi/(8447*s3*T))**0.125
    a = eps**2 + 4/3*eps**4 + 148/45*eps**6; b = eps**2 + 5/3*eps**4 + 167/36*eps**6
    ph = -phi if reverse else phi
    HF = lambda t: 2*Q + 2*eps*np.cos(t + ph)*Jx + a*(P - Q) - 2*b*np.cos(2*(t + ph))*X
    # basis of Hermitian generators for the averaging map
    basis = [np.eye(4)[:, [i]] @ np.eye(4)[[j], :] for i in range(4) for j in range(4)]
    def rhs(t, y):
        U = y[:16].reshape(4, 4); dU = -1j*HF(t) @ U
        acc = [U.conj().T @ Bm @ U for Bm in basis]
        return np.concatenate([dU.ravel()] + [A.ravel() for A in acc])
    sol = solve_ivp(rhs, (0, 2*np.pi), np.concatenate([np.eye(4).ravel(), np.zeros(16*16, complex)]), method='DOP853', rtol=1e-13, atol=1e-14)
    UF = sol.y[:16, -1].reshape(4, 4); Lmap = sol.y[16:, -1].reshape(16, 4, 4)      # <E_ij>_F
    avg = lambda Y: np.einsum('k,kab->ab', Y.ravel(), Lmap)
    Jr = Jx*np.cos(alpha) + Jy*np.sin(alpha); K = P @ Jr @ P
    D = lambda s: sla.expm(-1j*np.pi*s*K/2); Rr = lambda s: sla.expm(-1j*np.pi*s*Jr)
    U = np.eye(4, dtype=complex)
    # C(s) is a quadratic trigonometric function of s: build <C(s)>_F from the exact P-block rotation
    for k in range(N):
        s = (k + 0.5)/N
        if reverse:
            Dm = D(1 - s); C = np.pi*Dm.conj().T @ Jr @ Dm + np.pi/2*K
        else:
            Dm = D(s); C = -np.pi*Dm.conj().T @ Jr @ Dm - np.pi/2*K
        G = avg(-gam*np.pi*Q + C)/T
        U = UF @ sla.expm(-1j*G) @ U
    S = Rr(1) @ D(1)
    return S @ U if not reverse else U @ S.conj().T
def prob(N):
    A = np.eye(4, dtype=complex); B = np.eye(4, dtype=complex)
    for r in range(3):
        A = pulse(N, alphas[r], phis[r], False) @ A
    for r in (2, 1, 0):
        B = pulse(N, alphas[r], phis[r], True) @ B
    w = np.trace(B.conj().T @ A @ P)/2
    return (1 - w.real)/2, (1080*np.pi/(8447*s3*2*np.pi*N))**0.125
if __name__ == "__main__":
    for N in map(int, sys.argv[1:]):
        p, eps = prob(N); print(N, eps, repr(p), p - 1467167323933/5157470703125, flush=True)
