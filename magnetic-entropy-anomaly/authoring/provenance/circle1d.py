"""Exact (Hermite-truncated) finite-m stationary solution of the 1D step model on the circle, piecewise constant T.
Scaled: xi = x/sqrt(m), u = sqrt(m) v; each stripe has length L = pi/sqrt(m).  u d_xi p = d_u(u p + T d_u p).
sqrt(m) S_dot_m = q(0) (1/Th - 1/Tc) with q = int u^3 p du at the line x = 0 (hot side to the right), mass normalised."""
import numpy as np, scipy.linalg as sla, sys
from layer1d import maxw
def run(K, m, Th=2.0, Tc=1.0, Tr=1.5):
    n = np.arange(K); L = np.pi / np.sqrt(m)
    V = np.sqrt(Tr) * (np.diag(np.sqrt(n[1:]), 1) + np.diag(np.sqrt(n[1:]), -1))
    def A(T):
        t = T / Tr; M = -np.diag(n).astype(float); M[2:, :-2] -= (1 - t) * np.diag(np.sqrt((n[:-2] + 1) * (n[:-2] + 2))); return M
    blocks = []
    for T in (Th, Tc):      # stripe 0: hot (xi in [0,L]), stripe 1: cold (xi in [L, 2L])
        B = np.linalg.solve(V, A(T)); lam, W = np.linalg.eig(B)
        nz = np.abs(lam) > 1e-7
        cM = maxw(K, T / Tr)
        cJ = np.linalg.lstsq(B, cM, rcond=None)[0]          # current mode: c(xi) = xi cM + cJ  (B cJ = cM)
        blocks.append((lam[nz], W[:, nz], cM, cJ))
    def mode_vals(blk, xi0, xi1, xi):
        lam, W, cM, cJ = blk
        # columns: decaying modes anchored at the end they decay from, then Maxwellian, then current mode (anchored at xi0)
        e = np.where(lam.real < 0, np.exp(lam * (xi - xi0)), np.exp(lam * (xi - xi1)))
        return np.hstack([W * e[None, :], cM[:, None], ((xi - xi0) * cM + cJ)[:, None]])
    (h, c) = blocks
    nh, nc = h[1].shape[1] + 2, c[1].shape[1] + 2
    # continuity at xi = L (hot end -> cold start) and at xi = 2L == 0 (cold end -> hot start)
    E1 = np.hstack([mode_vals(h, 0, L, L), -mode_vals(c, L, 2 * L, L)])
    E2 = np.hstack([-mode_vals(h, 0, L, 0), mode_vals(c, L, 2 * L, 2 * L)])
    Mt = np.vstack([E1, E2])
    # null vector
    U, s, Vh = np.linalg.svd(Mt); a = Vh[-1].conj()
    # mass: int rho dxi = int c0(xi) dxi over both stripes, rho component = coefficient 0
    def mass(blk, xi0, xi1, amp):
        lam, W, cM, cJ = blk
        ints = np.where(lam.real < 0, (1 - np.exp(lam * (xi1 - xi0))) / (-lam), (1 - np.exp(-lam * (xi1 - xi0))) / lam)
        Lh = xi1 - xi0
        return (W[0] * ints) @ amp[:-2] + cM[0] * Lh * amp[-2] + (cM[0] * Lh ** 2 / 2 + cJ[0] * Lh) * amp[-1]
    Mtot = mass(h, 0, L, a[:nh]) + mass(c, L, 2 * L, a[nh:])
    a = a / (Mtot * np.sqrt(m))                       # physical mass = sqrt(m) int rho dxi = 1
    c0 = mode_vals(h, 0, L, 0) @ a[:nh]
    q = (np.sqrt(6) * c0[3] + 3 * c0[1]) * Tr ** 1.5
    J = c0[1] * np.sqrt(Tr)
    return float(np.real(q * (1 / Th - 1 / Tc))), float(np.abs(J)), s[-1] / s[0], s[-2] / s[0]
if __name__ == "__main__":
    for m in (0.01, 0.0025, 0.000625):
        for K in (40, 80):
            print(m, K, run(K, m), flush=True)
