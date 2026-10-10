"""Layer-resolved anomalous Hall conductivity of a slab, summed over the top half. It reproduces theta/2pi (the
Chern-Simons part) and NOT the full alpha_zz: the cross-gap part is a bulk polarization that this decomposition does
not see. usage: python3 surface_hall.py <L> <Nk> [J1 J2 m]"""
import sys

import numpy as np

from model import PARAMS
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..", "solution"))
import solve  # noqa: E402


def top_half(p, L, nk):
    g = (np.arange(nk) + 0.5) * 2 * np.pi / nk
    C = np.zeros(L)
    for kx in g:
        H, DX, DY = solve.slab_blocks(p, L, np.full(nk, kx), g)
        E, V = np.linalg.eigh(H)
        n = 2 * L
        Vh = np.conj(np.swapaxes(V, -1, -2))
        X, Y = Vh @ DX @ V, Vh @ DY @ V
        den = E[:, None, n:] - E[:, :n, None]
        Xo, Yo = X[:, :n, n:] / den, Y[:, :n, n:] / den
        F = 1j * (Xo @ np.conj(np.swapaxes(Yo, -1, -2)) - Yo @ np.conj(np.swapaxes(Xo, -1, -2)))
        Vo = V[:, :, :n]
        W = np.einsum("kia,kab,kib->ki", Vo, F, np.conj(Vo)).real
        C += W.reshape(nk, L, 4).sum(axis=(0, 2))
    C *= (2 * np.pi / nk) ** 2 / (2 * np.pi)
    return C[:L // 2].sum(), C.sum()


if __name__ == "__main__":
    L, nk = int(sys.argv[1]), int(sys.argv[2])
    p = dict(PARAMS)
    if len(sys.argv) > 3:
        p.update(J1=float(sys.argv[3]), J2=float(sys.argv[4]), m=float(sys.argv[5]))
    top, tot = top_half(p, L, nk)
    print(f"L={L} Nk={nk} top-half Hall conductivity = {top:.10f} e^2/h (slab total {tot:.2e})")
