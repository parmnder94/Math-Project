"""1D kinetic layer at a temperature step, moment method with ordered real Schur subspaces (stable)."""
import numpy as np, scipy.linalg as sla, sys
from math import lgamma
def maxw(K, t):
    c = np.zeros(K)
    for k2 in range(0, K, 2):
        kk = k2 // 2
        c[k2] = 1.0 if (t == 1 and kk == 0) else (0.0 if t == 1 else np.sign(t - 1) ** kk * np.exp(kk * np.log(abs(t - 1) / 2) - lgamma(kk + 1) + 0.5 * lgamma(k2 + 1)))
    return c
def run(K, Tm, Tp, Tr):
    n = np.arange(K)
    V = np.diag(np.sqrt(n[1:]), 1) + np.diag(np.sqrt(n[1:]), -1)
    def A(T):
        t = T / Tr; M = -np.diag(n).astype(float)
        M[2:, :-2] -= (1 - t) * np.diag(np.sqrt((n[:-2] + 1) * (n[:-2] + 2)))
        return M
    bases = []
    for T, side in ((Tm, -1), (Tp, +1)):
        B = np.linalg.solve(np.sqrt(Tr) * V, A(T))
        # stable subspace for the side: right side keeps Re<0 ; left side keeps Re>0
        sel = (lambda x: x.real < -1e-9) if side > 0 else (lambda x: x.real > 1e-9)
        Tsch, Z, sdim = sla.schur(B, output='real', sort=lambda re, im: sel(complex(re, im)))
        bases.append((Z[:, :sdim], maxw(K, T / Tr)))
    (Zm, cMm), (Zp, cMp) = bases
    Mat = np.hstack([Zm, -Zp])
    rhs = -(cMm / Tm - cMp / Tp)
    sol = np.linalg.lstsq(Mat, rhs, rcond=None)[0]
    resid = np.abs(Mat @ sol - rhs).max()
    c0 = Zp @ sol[Zm.shape[1]:] + cMp / Tp
    q = (np.sqrt(6) * c0[3] + 3 * c0[1]) * Tr ** 1.5
    P = (c0[0] + np.sqrt(2) * c0[2]) * Tr
    return q / P, resid, Zm.shape[1], Zp.shape[1]
if __name__ == "__main__":
    Tm, Tp, Tr = 1.0, 2.0, 1.5
    prev = None
    for K in [40, 80, 120, 160, 240, 320, 480, 640]:
        R, res, a, b = run(K, Tm, Tp, Tr)
        print(K, repr(R), f"res={res:.1e} modes={a},{b}", f"diff={R-prev:.3e}" if prev is not None else "", flush=True); prev = R
