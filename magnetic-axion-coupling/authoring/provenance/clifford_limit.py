"""Normalization and sign check: in the Clifford limit J1 = J2 = 0 (Neel mass only) the cross-gap part vanishes, so
the slab dM/dE result must equal theta/2pi from the Chern-Simons code. usage: python3 clifford_limit.py"""
import os
import sys

import numpy as np

from cs_stretched import theta_cs
from model import PARAMS
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "solution"))
import solve  # noqa: E402

if __name__ == "__main__":
    p = dict(PARAMS, J1=0.0, J2=0.0)
    th = theta_cs(p, 64, 0.0)[0]
    resp = []
    for L in (14, 18):
        resp.append((solve.magnetization(p, L, 1e-3)[0] - solve.magnetization(p, L, -1e-3)[0]) / 2e-3)
    alpha = 2 * np.pi * (resp[1] - resp[0]) / 4
    print(f"Clifford limit (J1=J2=0): theta/2pi = {th / (2 * np.pi):.10f}   slab dM/dE alpha_zz = {alpha:.10f}   "
          f"difference = {alpha - th / (2 * np.pi):.1e}")
