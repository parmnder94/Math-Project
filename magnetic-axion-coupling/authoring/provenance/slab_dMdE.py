"""Truth route: alpha_zz = dM_z/dE_z from the orbital magnetization of slabs in a uniform field E_z (same method as
solution/solve.py, with more layers and k points). usage: python3 slab_dMdE.py <Nk> <comma-separated L list>"""
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..", "solution"))
import solve  # noqa: E402

if __name__ == "__main__":
    solve.NK = int(sys.argv[1])
    Ls = [int(x) for x in sys.argv[2].split(",")]
    resp = []
    for L in Ls:
        mp, gap = solve.magnetization(solve.P, L, 1e-3)
        mm, _ = solve.magnetization(solve.P, L, -1e-3)
        resp.append((mp - mm) / 2e-3)
    for (L1, a1), (L2, a2) in zip(zip(Ls, resp), zip(Ls[1:], resp[1:])):
        print(f"Nk={solve.NK} L{L1}->L{L2}: alpha_zz = {2 * np.pi * (a2 - a1) / (L2 - L1):.12f} e^2/h", flush=True)
