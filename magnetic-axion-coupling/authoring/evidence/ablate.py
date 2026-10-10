"""Values produced by plausible shortcuts. usage: python3 ablate.py <out.json>"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "provenance"))
from cs_stretched import theta_cs  # noqa: E402
from gap_scan import min_gap  # noqa: E402
from model import PARAMS  # noqa: E402

TRUTH = float(json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tests", "truth.json")))["theta"])


def wrap(x):
    return float(np.angle(np.exp(1j * x)))


def main():
    out = {"reference theta": TRUTH,
           "quantized value pi assumed (time-reversal-symmetric TI)": np.pi,
           "sign convention reversed (A = -i<u|du> or left-handed k basis)": -TRUTH}
    out["exchange dropped (Neel mass only, J1 = J2 = 0)"] = wrap(theta_cs(dict(PARAMS, J1=0.0, J2=0.0), 128, 0.8)[0])
    out["orbital dependence of the exchange dropped (J1 = J2 = average)"] = (
        "gapless: min direct gap %.1e, theta undefined" % min_gap(dict(PARAMS, J1=0.075, J2=0.075))[0])
    for n in (64, 96, 128):
        out[f"uniform {n}^3 grid, same smooth gauge"] = wrap(theta_cs(PARAMS, n, 0.0)[0])
    json.dump(out, open(sys.argv[1], "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
