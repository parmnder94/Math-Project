"""Independent check: theta = pi + integral of the second Chern form along a gapped path.

The path starts at the time-reversal-symmetric strong TI (m = J1 = J2 = 0, theta = pi): first m is ramped from 0 to
0.03 with J = 0, then (J1, J2) from 0 to (0.20, -0.05). d theta = -(1/4 pi) tr(F ^ F) on T^3 x path, with the
gauge-free curvature operator F_ab = i P [d_a P, d_b P] P (a, b in lambda, kx, ky, kz). Stretched grids (Jacobian
included), Gauss-Legendre in lambda on each segment.

usage: python3 path_second_chern.py <N> <Gauss nodes per segment> <stretch c>
"""
import itertools
import sys

import numpy as np

from model import PARAMS, K, dham, ham, projector_and_derivatives, s0, sy, sz

PERMS = [(p, np.linalg.det(np.eye(4)[list(p)])) for p in itertools.permutations(range(4))]


def stretch(n, k0, c):
    s = np.arange(n) * 2 * np.pi / n
    return k0 + (s - k0) - c * np.sin(s - k0), 1 - c * np.cos(s - k0)


def dtheta(q, dlam, n, c):
    (kx, fx), (ky, fy), (kz, fz) = [stretch(n, k0, c) for k0 in (0.0, 0.0, np.pi)]
    tot, gmin = 0.0, np.inf
    for ix in range(n):
        KX, KY, KZ = np.meshgrid(kx[ix:ix + 1], ky, kz, indexing="ij")
        jac = fx[ix] * fy[None, :, None] * fz[None, None, :]
        H = ham(KX, KY, KZ, q)
        P, dP, E = projector_and_derivatives(H, [np.broadcast_to(dlam, H.shape)] + dham(KX, KY, KZ, q))
        gmin = min(gmin, (E[..., 2] - E[..., 1]).min())
        F = {}
        for a in range(4):
            for b in range(a + 1, 4):
                F[(a, b)] = 1j * P @ (dP[a] @ dP[b] - dP[b] @ dP[a]) @ P
                F[(b, a)] = -F[(a, b)]
        d = sum(sg * np.einsum("...ij,...ji->...", F[(a, b)], F[(c_, d_)]) for (a, b, c_, d_), sg in PERMS)
        tot += (d * jac).sum()
    return -(tot * (2 * np.pi / n) ** 3 / 4) / (4 * np.pi), gmin


def theta_path(n, nl, c):
    x, w = np.polynomial.legendre.leggauss(nl)
    u, w = (x + 1) / 2, w / 2
    dm = PARAMS["m"] * K(s0, sy)
    dJ = K(sz, PARAMS["J1"] * (s0 + sz) / 2 + PARAMS["J2"] * (s0 - sz) / 2)
    total, gmin = 0.0, np.inf
    for seg in (0, 1):
        for ui, wi in zip(u, w):
            q = dict(PARAMS)
            if seg == 0:
                q.update(m=ui * PARAMS["m"], J1=0.0, J2=0.0)
                v, g = dtheta(q, dm, n, c)
            else:
                q.update(J1=ui * PARAMS["J1"], J2=ui * PARAMS["J2"])
                v, g = dtheta(q, dJ, n, c)
            total += wi * v
            gmin = min(gmin, g)
    return np.pi + total, gmin


if __name__ == "__main__":
    n, nl, c = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    th, g = theta_path(n, nl, c)
    print(f"N={n} nodes={nl} c={c} theta={th.real:.12f} imag={th.imag:.1e} min_gap_on_path_grid={g:.4f}", flush=True)
