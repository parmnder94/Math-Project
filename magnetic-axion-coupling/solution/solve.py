"""Reference solution for magnetic-axion-coupling: the full orbital magnetoelectric coefficient alpha_zz.

alpha_zz = dP_z/dB_z = dM_z/dE_z (Maxwell relation of the insulating ground state). It contains the topological
Chern-Simons part theta e^2/(2 pi h) and a non-topological cross-gap part; here the second is -0.036 e^2/h, so the
axion angle alone does not give the answer.

Route used here (converse effect): a slab of L layers along z (open in z, periodic in x, y) is placed in a uniform
field E_z, i.e. the on-site energy +e E_z z_l of an electron (charge -e) in layer l. Its orbital magnetization per area
follows from the multiband insulator formula
    M_z = (e/hbar) Im sum_{n occ, m unocc} (E_n + E_m) <n|d_x H|m><m|d_y H|n> / (E_m - E_n)^2   (per (2 pi)^2 d^2k),
valid because the slab's total Chern number is zero. dM_z/dE_z per area grows linearly with L once the surfaces are
converged; the slope per layer is alpha_zz (lattice constant 1). In units with e = hbar = 1 the result is converted to
e^2/h by a factor 2 pi. The slope converges exponentially in L (L = 14 -> 18 is within 1e-8 of L = 18 -> 22) and in
the k grid (48^2 points).
Run time about 15 s.
"""
import json
import os

import numpy as np

P = dict(M0=0.5, A1=0.6, A2=0.9, B1=0.3, B2=0.6, m=0.25, J1=0.35, J2=-0.25)
NK, LAYERS, DE = 48, (14, 18), 1e-3

s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1.0, -1.0]).astype(complex)
K = np.kron  # sigma (spin) (x) tau (orbital)


def hoppings(p):
    """t_0 and t_{+d} with H(k) = sum_R t_R e^{i k.R}, t_{-d} = t_{+d}^dagger."""
    t0 = (p["M0"] - 2 * p["B1"] - 4 * p["B2"]) * K(s0, sz) + p["m"] * K(s0, sy) \
        + K(sz, p["J1"] * (s0 + sz) / 2 + p["J2"] * (s0 - sz) / 2)
    tx = p["B2"] * K(s0, sz) + p["A2"] / 2j * K(sx, sx)
    ty = p["B2"] * K(s0, sz) + p["A2"] / 2j * K(sy, sx)
    tz = p["B1"] * K(s0, sz) + p["A1"] / 2j * K(sz, sx)
    return t0, tx, ty, tz


def slab_blocks(p, L, kx, ky):
    t0, tx, ty, tz = hoppings(p)
    ex, ey = np.exp(1j * kx)[:, None, None], np.exp(1j * ky)[:, None, None]
    h = t0 + tx * ex + tx.conj().T / ex + ty * ey + ty.conj().T / ey
    hx = 1j * (tx * ex - tx.conj().T / ex)
    hy = 1j * (ty * ey - ty.conj().T / ey)
    n = 4 * L
    H = np.zeros((len(kx), n, n), complex)
    DX, DY = np.zeros_like(H), np.zeros_like(H)
    for l in range(L):
        s = slice(4 * l, 4 * l + 4)
        H[:, s, s], DX[:, s, s], DY[:, s, s] = h, hx, hy
        if l + 1 < L:
            s2 = slice(4 * l + 4, 4 * l + 8)
            H[:, s, s2], H[:, s2, s] = tz, tz.conj().T
    return H, DX, DY


def magnetization(p, L, Ez):
    g = (np.arange(NK) + 0.5) * 2 * np.pi / NK
    U = np.kron(np.diag(Ez * (np.arange(L) - (L - 1) / 2)), np.eye(4))
    total, gap = 0.0, np.inf
    for kx in g:
        H, DX, DY = slab_blocks(p, L, np.full(NK, kx), g)
        E, V = np.linalg.eigh(H + U)
        n = 2 * L
        gap = min(gap, (E[:, n] - E[:, n - 1]).min())
        Vh = np.conj(np.swapaxes(V, -1, -2))
        X, Y = Vh @ DX @ V, Vh @ DY @ V
        Eo, Eu = E[:, :n, None], E[:, None, n:]
        total += np.imag(((Eo + Eu) * X[:, :n, n:] * np.swapaxes(Y[:, n:, :n], -1, -2) / (Eu - Eo) ** 2).sum())
    return total * (2 * np.pi / NK) ** 2 / (2 * np.pi) ** 2, gap


def main():
    resp = []
    for L in LAYERS:
        mp, gap = magnetization(P, L, DE)
        mm, _ = magnetization(P, L, -DE)
        resp.append((mp - mm) / (2 * DE))
    slope = (resp[1] - resp[0]) / (LAYERS[1] - LAYERS[0])
    alpha = 2 * np.pi * slope                      # units of e^2/h
    alpha = float((alpha + 0.5) % 1.0 - 0.5)
    print(f"slab gap {gap:.3f}  dM/dE per area {resp}  alpha_zz = {alpha:.10f} e^2/h")
    os.makedirs("/app/output", exist_ok=True)
    with open("/app/output/answer.json", "w") as fh:
        json.dump({"alpha_zz": alpha}, fh)


if __name__ == "__main__":
    main()
