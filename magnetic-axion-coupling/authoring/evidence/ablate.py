"""Values produced by plausible shortcuts. usage: python3 ablate.py <out.json>"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "provenance"))
from cs_stretched import theta_cs  # noqa: E402
from model import PARAMS  # noqa: E402
from surface_hall import top_half  # noqa: E402

TRUTH = float(json.load(open(os.path.join(HERE, "..", "..", "tests", "truth.json")))["alpha_zz"])


def main():
    th = theta_cs(PARAMS, 64, 0.0)[0]
    out = {
        "reference alpha_zz (e^2/h)": TRUTH,
        "axion angle only, theta/2pi (Chern-Simons part, what a theta calculation gives)": th / (2 * np.pi),
        "surface (top-half layer-resolved) Hall conductivity of a slab": top_half(PARAMS, 24, 64)[0],
        "sign or field orientation reversed": -TRUTH,
        "half-quantized value assumed": 0.5,
        "Hofstadter supercell q = 32, forward difference, not extrapolated": 0.15491451,
    }
    json.dump(out, open(sys.argv[1], "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
