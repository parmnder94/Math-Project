"""Reference solution for magnetic-axion-coupling.

theta is the Chern-Simons integral of the occupied (lower two) bands in a smooth periodic gauge.

1. Gap and topology. The direct gap between bands 2 and 3 stays open (minimum about 0.050 eV, on a ring of radius
   about 0.21 around (0, 0, pi)). The occupied bundle is trivial, so a smooth periodic gauge exists.
2. Smooth gauge. Project the tau_y = -1 trial states |s> (x) |tau_y = -1> (s = up, down) onto the occupied subspace and
   Loewdin-orthonormalize: U = P Phi (Phi^+ P Phi)^(-1/2). The overlap Phi^+ P Phi has smallest eigenvalue 0.43
   everywhere, so U is smooth and periodic.
3. Grid. The small-gap ring makes the integrand sharp. Each coordinate is mapped by k = k0 + (s - k0) - c sin(s - k0)
   (c = 0.8; k0 = 0, 0, pi), a smooth periodic orientation-preserving bijection that concentrates points near the ring.
   The Chern-Simons density is a 3-form, so it is evaluated directly in the s coordinates (no Jacobian), with FFT
   derivatives on the uniform s grid. Convergence is spectral: N = 128 agrees with N = 160 to about 1e-11.
4. theta = -(1/4 pi) * integral, wrapped into (-pi, pi].
Run time about 35 s; peak memory about 1.3 GB (the eigen-decomposition is done in slabs).
"""
import json
import os

import numpy as np

P = dict(M0=0.28, A1=0.22, A2=0.40, B1=0.08, B2=0.50, m=0.03, J1=0.20, J2=-0.05)
N, C, K0 = 128, 0.8, (0.0, 0.0, np.pi)

s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1.0, -1.0]).astype(complex)
kron = np.kron  # sigma (spin) (x) tau (orbital)


def ham(kx, ky, kz, p):
    Mk = p["M0"] - 2 * p["B1"] * (1 - np.cos(kz)) - 2 * p["B2"] * (2 - np.cos(kx) - np.cos(ky))
    c = lambda a: np.asarray(a)[..., None, None]
    return (c(Mk) * kron(s0, sz) + c(p["A1"] * np.sin(kz)) * kron(sz, sx)
            + c(p["A2"] * np.sin(kx)) * kron(sx, sx) + c(p["A2"] * np.sin(ky)) * kron(sy, sx)
            + p["m"] * kron(s0, sy) + kron(sz, p["J1"] * (s0 + sz) / 2 + p["J2"] * (s0 - sz) / 2))


def stretch(n, k0, c):
    s = np.arange(n) * 2 * np.pi / n
    return k0 + (s - k0) - c * np.sin(s - k0)


def smooth_gauge(p, n, c, k0):
    t = np.array([1, -1j]) / np.sqrt(2)                       # tau_y eigenvector with eigenvalue -1
    phi = np.stack([kron([1, 0], t), kron([0, 1], t)], axis=1)  # 4 x 2 trial states
    ks = [stretch(n, k0[i], c) for i in range(3)]
    U = np.empty((n, n, n, 4, 2), dtype=complex)
    gap, smin = np.inf, np.inf
    for ix in range(n):                                       # slabs keep memory low
        KY, KZ = np.meshgrid(ks[1], ks[2], indexing="ij")
        E, V = np.linalg.eigh(ham(np.full_like(KY, ks[0][ix]), KY, KZ, p))
        gap = min(gap, (E[..., 2] - E[..., 1]).min())
        Vo = V[..., :2]
        B = phi.conj().T @ Vo                                 # <phi_a|u_n>
        w, u = np.linalg.eigh(B @ np.conj(np.swapaxes(B, -1, -2)))
        smin = min(smin, w.min())
        U[ix] = Vo @ np.conj(np.swapaxes(B, -1, -2)) @ (u @ (w[..., None] ** -0.5 * np.conj(np.swapaxes(u, -1, -2))))
    return U, gap, smin


def chern_simons_theta(U, n):
    kf = np.fft.fftfreq(n, 1 / n)

    def der(X, ax):
        shape = [1] * X.ndim
        shape[ax] = n
        return np.fft.ifft(np.fft.fft(X, axis=ax) * (1j * kf).reshape(shape), axis=ax)

    Uh = np.conj(np.swapaxes(U, -1, -2))
    A = [1j * (Uh @ der(U, a)) for a in range(3)]
    del Uh
    eps = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}
    total = 0
    for (i, j, k), s in eps.items():
        total += s * np.einsum("...ab,...ba->...", A[i], der(A[k], j)).sum()
        total += -(2j / 3) * s * np.einsum("...ab,...bc,...ca->...", A[i], A[j], A[k]).sum()
    return -(total * (2 * np.pi / n) ** 3) / (4 * np.pi)


def main():
    U, gap, smin = smooth_gauge(P, N, C, K0)
    th = chern_simons_theta(U, N)
    assert abs(th.imag) < 1e-10
    theta = float(np.angle(np.exp(1j * th.real)))
    print(f"N={N} c={C}  min direct gap on grid={gap:.4f}  min trial overlap={smin:.3f}  theta={theta:.12f}")
    os.makedirs("/app/output", exist_ok=True)
    with open("/app/output/answer.json", "w") as fh:
        json.dump({"theta": theta}, fh)


if __name__ == "__main__":
    main()
