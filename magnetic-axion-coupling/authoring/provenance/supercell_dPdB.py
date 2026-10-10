"""Independent route: alpha_zz = dP_z/dB_z from magnetic (Hofstadter) supercells.

Landau gauge A = (0, B x, 0); flux phi = +-1/q per xy plaquette (in units of h/e) on a q x 1 x 1 supercell. With the
Peierls rule of the instruction, the block <R|H|R+y> at column x gets exp(+2 pi i phi x). P_z (charge -e, units e per
cell volume) is -Phi/(2 pi q), Phi the Berry phase of the 2q lowest bands along k_z (discrete Wilson loop, Nz points),
averaged over the magnetic Brillouin zone (k_x, k_y in [0, 2 pi/q), using the magnetic translation symmetry).
alpha_q = q (P(+1/q) - P(-1/q)) / 2 in units of e^2/h (defined mod 1/2, branch fixed by the forward difference against
the zero-field cell polarization), and alpha = lim alpha_q with corrections O(1/q^2).

usage: python3 supercell_dPdB.py <q list> <Nz>
"""
import sys

import numpy as np

from model import PARAMS, hoppings


def h_super(p, q, flux, kx, ky, kzs):
    t0, tx, ty, tz = hoppings(p)
    n = 4 * q
    H = np.zeros((len(kzs), n, n), complex)
    ez = np.exp(1j * kzs)[:, None, None]
    for x in range(q):
        s = slice(4 * x, 4 * x + 4)
        ph = np.exp(1j * (ky + 2 * np.pi * flux * x))
        H[:, s, s] = (t0 + ty * ph + ty.conj().T / ph)[None] + tz[None] * ez + tz.conj().T[None] / ez
        xn = (x + 1) % q
        f = np.exp(1j * kx * q) if x == q - 1 else 1.0
        sn = slice(4 * xn, 4 * xn + 4)
        H[:, s, sn] += tx * f
        H[:, sn, s] += (tx * f).conj().T
    return H


def berry_z(p, q, flux, kx, ky, nz):
    kzs = np.arange(nz) * 2 * np.pi / nz
    E, V = np.linalg.eigh(h_super(p, q, flux, kx, ky, kzs))
    U = V[..., :2 * q]
    M = np.conj(np.swapaxes(U, -1, -2)) @ np.roll(U, -1, axis=0)
    return -np.sum(np.angle(np.linalg.det(M))), (E[:, 2 * q] - E[:, 2 * q - 1]).min()


def mean_phase(p, q, flux, nk, nz):
    vals, gmin = [], np.inf
    for a in range(nk):
        for b in range(nk):
            v, g = berry_z(p, q, flux, (a + 0.5) * 2 * np.pi / q / nk, (b + 0.5) * 2 * np.pi / q / nk, nz)
            vals.append(v)
            gmin = min(gmin, g)
    return np.unwrap(np.array(vals)).mean(), gmin


def cell_phase(p, nk, nz):
    vals = np.array([[berry_z(p, 1, 0.0, a * 2 * np.pi / nk, b * 2 * np.pi / nk, nz)[0] for b in range(nk)]
                     for a in range(nk)])
    return np.unwrap(np.unwrap(vals, axis=0), axis=1).mean()


def alpha_q(p, q, nz, nk=2, phi0=None):
    phi0 = cell_phase(p, 24, nz) if phi0 is None else phi0
    vp, g = mean_phase(p, q, 1.0 / q, nk, nz)
    vm, _ = mean_phase(p, q, -1.0 / q, nk, nz)
    fwd = (-(vp - q * phi0) / (2 * np.pi) + 0.5) % 1 - 0.5
    d = -(vp - vm) / (2 * np.pi)
    d += np.round(2 * fwd - d)
    return d / 2, g


if __name__ == "__main__":
    nz = int(sys.argv[2])
    phi0 = cell_phase(PARAMS, 24, nz)
    for q in [int(x) for x in sys.argv[1].split(",")]:
        a, g = alpha_q(PARAMS, q, nz, phi0=phi0)
        print(f"q={q} Nz={nz} alpha_q={a:.8f} gap={g:.3f}", flush=True)
