"""Plausible wrong routes, each changing one step of the reference computation (60-digit arithmetic).
Usage: python3 ablate.py <solution_dir> <out.json>"""
import importlib.util, json, os, sys
import mpmath as mp
import sympy as sp
spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
bq, bx = sp.Rational(23539, 2160), sp.Rational(8447, 540)
P = lambda **kw: S.probability(bq, bx, **kw)[0]
C = lambda **kw: S.first_correction(bq, bx, **kw)
ref_p, ref_c = P(), C()
routes = {
    "reference": (ref_p, ref_c),
    # routes for p_inf
    "endpoint frame D_r(1) set to the identity": (P(use_D=False), C(use_D=False)),
    "slow connection omitted": (P(connection=0), C(connection=0)),
    "connection from R_r only (-pi K_r)": (P(connection=1), C(connection=1)),
    "reversed pulse keeps the carrier phase": (P(flip_phase=False), C(flip_phase=False)),
    "(Q-P) part of the eighth-order block omitted": (S.probability(sp.Integer(0), bx)[0], S.first_correction(sp.Integer(0), bx)),
    # routes for c
    "stroboscopic gauge: no kick at the pulse boundaries (c = 0)": (ref_p, mp.mpf(0)),
    "kick only at the start of each pulse": (ref_p, C(kick_ends="start")),
    "kick only at the end of each pulse": (ref_p, C(kick_ends="end")),
    "kick denominators without the 2Q gap (-l only)": (ref_p, C(kick_den=lambda Ea, Eb, l: -l)),
}
out = {}
for k, (p, c) in routes.items():
    out[k] = {"p": float(p), "c": float(c), "abs_error_p": float(abs(p - ref_p)), "abs_error_c": float(abs(c - ref_c))}
    print(f"{k:62s} p = {float(p):.10f} ({out[k]['abs_error_p']:.1e})   c = {float(c):+.10f} ({out[k]['abs_error_c']:.1e})")
out["naive numerical slope (p_N - p_inf)/eps_N at N = 1e40 (c2 about -15)"] = {"c_estimate": float(ref_c) - 15*6.1e-6, "abs_error_c": 15*6.1e-6}
json.dump(out, open(sys.argv[2], "w"), indent=2)
