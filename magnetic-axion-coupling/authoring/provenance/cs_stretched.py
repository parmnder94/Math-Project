"""theta from the Chern-Simons definition in a Loewdin-projected trial gauge on stretched periodic grids.

usage: python3 cs_stretched.py <N list, e.g. 96,128,160> <stretch c (0 = uniform grid)> [trial: tauy | entangled]
"""
import sys

import numpy as np

from model import PARAMS, ham


def stretch(n, k0, c):
    s = np.arange(n) * 2 * np.pi / n
    return k0 + (s - k0) - c * np.sin(s - k0)


def trial(kind):
    t1 = np.array([1, -1j]) / np.sqrt(2)          # tau_y = -1
    t2 = np.array([1, 1j]) / np.sqrt(2)           # tau_y = +1
    if kind == "tauy":
        return np.stack([np.kron([1, 0], t1), np.kron([0, 1], t1)], axis=1)
    a = np.kron([1, 0], t1) + 0.4 * np.kron([0, 1], t2)
    b = np.kron([0, 1], t1) - 0.4 * np.kron([1, 0], t2)
    return np.linalg.qr(np.stack([a, b], axis=1))[0]


def theta_cs(p, n, c, kind="tauy", k0=(0.0, 0.0, np.pi)):
    ks = [stretch(n, k0[i], c) for i in range(3)]
    KX, KY, KZ = np.meshgrid(*ks, indexing="ij")
    E, V = np.linalg.eigh(ham(KX, KY, KZ, p))
    del KX, KY, KZ
    gap = (E[..., 2] - E[..., 1]).min()
    Vo = V[..., :2]
    del V, E
    phi = trial(kind)
    B = np.conj(phi.T) @ Vo
    w, u = np.linalg.eigh(B @ np.conj(np.swapaxes(B, -1, -2)))
    U = Vo @ np.conj(np.swapaxes(B, -1, -2)) @ (u @ (w[..., None] ** -0.5 * np.conj(np.swapaxes(u, -1, -2))))
    smin = w.min()
    del Vo, B, u, w
    kf = np.fft.fftfreq(n, 1 / n)

    def der(X, ax):
        sh = [1] * X.ndim
        sh[ax] = n
        return np.fft.ifft(np.fft.fft(X, axis=ax) * (1j * kf).reshape(sh), axis=ax)

    Uh = np.conj(np.swapaxes(U, -1, -2))
    A = [1j * (Uh @ der(U, a)) for a in range(3)]
    del U, Uh
    eps = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}
    tot = 0
    for (i, j, k), s in eps.items():
        tot += s * np.einsum("...ab,...ba->...", A[i], der(A[k], j)).sum()
        tot += -(2j / 3) * s * np.einsum("...ab,...bc,...ca->...", A[i], A[j], A[k]).sum()
    th = -(tot * (2 * np.pi / n) ** 3) / (4 * np.pi)
    return th.real, th.imag, smin, gap


if __name__ == "__main__":
    c = float(sys.argv[2])
    kind = sys.argv[3] if len(sys.argv) > 3 else "tauy"
    for n in [int(x) for x in sys.argv[1].split(",")]:
        th = theta_cs(PARAMS, n, c, kind)
        print(f"N={n} c={c} trial={kind} theta={th[0]:.12f} wrapped={np.angle(np.exp(1j * th[0])):.12f} "
              f"imag={th[1]:.1e} min_overlap={th[2]:.3f} grid_gap={th[3]:.4f}", flush=True)
