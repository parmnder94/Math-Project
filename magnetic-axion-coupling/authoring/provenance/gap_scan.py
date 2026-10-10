"""Direct gap of the instruction's model: global minimum and its location, the gap along the path used in
path_second_chern.py, and the gap for parameter changes used in the ablations.  usage: python3 gap_scan.py"""
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
    out = {}
    g, k = min_gap(PARAMS)
    out["model: min direct gap (eV) and location"] = [g, k, "radius in (kx,ky) = %.4f" % np.hypot(k[0], k[1])]
    path = []
    for lam in np.linspace(0, 1, 21):
        q = dict(PARAMS)
        a, b = min(1, 2 * lam), max(0, 2 * lam - 1)
        q.update(m=a * PARAMS["m"], J1=b * PARAMS["J1"], J2=b * PARAMS["J2"])
        path.append([round(float(lam), 3), round(min_gap(q, 40)[0], 5)])
    out["gap along the path (lambda, gap): m ramped on [0,1/2], exchange on [1/2,1]"] = path
    jp = (PARAMS["J1"] + PARAMS["J2"]) / 2
    out["J1 = J2 = (J1+J2)/2: min gap"] = min_gap(dict(PARAMS, J1=jp, J2=jp))[0]
    out["exchange doubled: min gap"] = min_gap(dict(PARAMS, J1=2 * PARAMS["J1"], J2=2 * PARAMS["J2"]))[0]
    print(json.dumps(out, indent=1))
