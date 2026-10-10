"""Bloch Hamiltonian of the instruction, its k-derivatives, and the occupied-projector derivative."""
import numpy as np

PARAMS = dict(M0=0.28, A1=0.22, A2=0.40, B1=0.08, B2=0.50, m=0.03, J1=0.20, J2=-0.05)
s0 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex)
sy = np.array([[0, -1j], [1j, 0]])
sz = np.diag([1.0, -1.0]).astype(complex)
K = np.kron  # sigma (spin) (x) tau (orbital)


def _c(a):
    return np.asarray(a)[..., None, None]


def ham(kx, ky, kz, p):
    Mk = p["M0"] - 2 * p["B1"] * (1 - np.cos(kz)) - 2 * p["B2"] * (2 - np.cos(kx) - np.cos(ky))
    return (_c(Mk) * K(s0, sz) + _c(p["A1"] * np.sin(kz)) * K(sz, sx)
            + _c(p["A2"] * np.sin(kx)) * K(sx, sx) + _c(p["A2"] * np.sin(ky)) * K(sy, sx)
            + p["m"] * K(s0, sy) + K(sz, p["J1"] * (s0 + sz) / 2 + p["J2"] * (s0 - sz) / 2))


def dham(kx, ky, kz, p):
    dx = _c(-2 * p["B2"] * np.sin(kx)) * K(s0, sz) + _c(p["A2"] * np.cos(kx)) * K(sx, sx)
    dy = _c(-2 * p["B2"] * np.sin(ky)) * K(s0, sz) + _c(p["A2"] * np.cos(ky)) * K(sy, sx)
    dz = _c(-2 * p["B1"] * np.sin(kz)) * K(s0, sz) + _c(p["A1"] * np.cos(kz)) * K(sz, sx)
    return [dx, dy, dz]


def projector_and_derivatives(H, dHs, nocc=2):
    """P onto the lowest nocc bands and dP = sum_{n occ, m unocc} |m><m|dH|n><n|/(E_n - E_m) + h.c."""
    E, V = np.linalg.eigh(H)
    Vh = np.conj(np.swapaxes(V, -1, -2))
    den = E[..., nocc:, None] - E[..., None, :nocc]          # E_m - E_n  (m unocc, n occ)
    out = []
    for dH in dHs:
        M = Vh @ dH @ V
        Y = np.zeros_like(M)
        Y[..., nocc:, :nocc] = -M[..., nocc:, :nocc] / den
        Y[..., :nocc, nocc:] = np.conj(np.swapaxes(Y[..., nocc:, :nocc], -1, -2))
        out.append(V @ Y @ Vh)
    Vo = V[..., :nocc]
    return Vo @ np.conj(np.swapaxes(Vo, -1, -2)), out, E
