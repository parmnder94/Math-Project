"""Direct and indirect gap of the instruction's model.  usage: python3 gap_scan.py"""
import json

import numpy as np
from scipy.optimize import minimize

from model import PARAMS, ham


def gap_at(k, p):
    E = np.linalg.eigvalsh(ham(*[np.array(x) for x in k], p))
    return E[2] - E[1]


def min_gap(p, n=48):
    g = np.arange(n) * 2 * np.pi / n - np.pi
    KX, KY, KZ = np.meshgrid(g, g, g, indexing="ij")
    E = np.linalg.eigvalsh(ham(KX, KY, KZ, p))
    G = E[..., 2] - E[..., 1]
    starts = np.argsort(G, axis=None)[:8]
    best = None
    for s in starts:
        i = np.unravel_index(s, G.shape)
        r = minimize(gap_at, [g[i[0]], g[i[1]], g[i[2]]], args=(p,), method="Nelder-Mead",
                     options=dict(xatol=1e-10, fatol=1e-13, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return float(best.fun), [float(x) for x in best.x]


if __name__ == "__main__":
    g, k = min_gap(PARAMS)
    E = np.linalg.eigvalsh(ham(*np.meshgrid(*(np.arange(40) * 2 * np.pi / 40,) * 3, indexing="ij"), PARAMS))
    out = {"min direct gap (eV) and location": [g, k],
           "indirect gap on a 40^3 grid (eV)": float(E[..., 2].min() - E[..., 1].max())}
    print(json.dumps(out, indent=1))
