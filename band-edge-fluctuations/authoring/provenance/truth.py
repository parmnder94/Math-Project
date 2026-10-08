"""Ground truth for band-edge-fluctuations, with cross-checks.

Run from this directory:  python3 truth.py <out_dir>
  A. Primary: 4x4-linearization engine (lib/engine4c.py), every term of total power <= 4 resummed exactly by
     marked transfer operators (+ Moebius inversion), remaining terms summed over primitive classes up to
     length 13 (convergence shown with lengths 11, 12, 13).
  B. Independent linearization: 6x6 engine (lib/edge_assemble.py, a different self-adjoint linearization and
     different code for the marked transfers) against the 4x4 engine with the same exact/enumerated split.
  C. End-to-end: the exact assembly against plain enumeration (lib/cross_assemble.py logic) at a strongly
     convergent point, where enumeration alone is accurate.
Per-length brute-force validations of every transfer construction are in validation/ (logs in
../evidence/validation_logs.txt). Monte Carlo of the finite networks: ../evidence/mc_summary.json.
"""
import json, sys, time, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "lib"))
import numpy as np
from engine4c import assemble4
from edge_assemble import assemble as assemble6

B, LAM, X0, ETA = 1.5, 0.6, 5.6, 0.06
out_dir = sys.argv[1] if len(sys.argv) > 1 else "."
report = {}
t0 = time.time()
conv = {}
for lmax in (11, 12, 13):
    conv[lmax] = assemble4(B, LAM, X0, ETA, lmax=lmax, lnmax=6 * lmax, exact_max=4)
    print("A", lmax, conv[lmax], round(time.time() - t0), flush=True)
truth = conv[13]
report["A_convergence"] = {str(l): v for l, v in conv.items()}
report["A_rel_change_12_to_13"] = {k: abs(conv[13][k] - conv[12][k]) / abs(conv[13][k]) for k in truth}

r6, _ = assemble6(B, LAM, X0, ETA, lmax=11, lnmax=66)
r4 = assemble4(B, LAM, X0, ETA, lmax=11, lnmax=66, exact_max=3)
report["B_K6_vs_K4"] = {k: {"K6": r6[k], "K4": r4[k], "rel_diff": abs(r6[k] - r4[k]) / abs(r4[k])} for k in r4}
print("B", report["B_K6_vs_K4"], flush=True)

from cross_assemble import TwoPoint
zc = complex(3.0, 7.0); Sx = TwoPoint(B, LAM, zc)
tr = np.einsum('n,nii->', Sx.w, Sx.G[:, :2, :2]).real
acc = Sx.run2(11, 66, verbose=False)
enum = {"z_inf": tr / 2, "mean_correction": acc['M'] - tr, "variance": acc['V'], "third_cumulant": acc['K3']}
ex = assemble4(B, LAM, zc.real, zc.imag, lmax=8, lnmax=48, exact_max=4)
report["C_convergent_point"] = {"z": str(zc), "rows": {k: {"enumeration": enum[k], "exact": ex[k], "rel_diff": abs(enum[k] - ex[k]) / abs(ex[k])} for k in ex}}
print("C", report["C_convergent_point"], flush=True)

report["seconds"] = time.time() - t0
json.dump(truth, open(os.path.join(out_dir, "truth.json"), "w"), indent=2)
json.dump(report, open(os.path.join(out_dir, "truth_report.json"), "w"), indent=2)
print("truth", truth)
