"""Plausible wrong routes, each changing one step of the reference computation (40-digit arithmetic).
Usage: python3 ablate.py <solution_dir> <out.json>"""
import importlib.util, json, os, sys
import mpmath as mp
import sympy as sp
spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
bq, bx = sp.Rational(23539, 2160), sp.Rational(8447, 540)
ref = S.probability(bq, bx)[0]
routes = {
    "reference": S.probability(bq, bx)[0],
    "endpoint frame D_r(1) set to the identity": S.probability(bq, bx, use_D=False)[0],
    "slow connection omitted": S.probability(bq, bx, connection=0)[0],
    "connection from R_r only (-pi K_r, D_r frame term dropped)": S.probability(bq, bx, connection=1)[0],
    "reversed pulse keeps the carrier phase (no phi -> -phi)": S.probability(bq, bx, flip_phase=False)[0],
    "carrier phases ignored in the effective coupling": S.probability(bq, bx, carrier=False)[0],
    "(Q-P) part of the eighth-order block omitted": S.probability(sp.Integer(0), bx)[0],
    "eighth-order coupling omitted": S.probability(bq, sp.Integer(0))[0],
    "B_8 off-diagonal coefficient 1e-6 too large": S.probability(bq, bx*(1 + sp.Rational(1, 10**6)))[0],
}
out = {k: {"p": float(v), "abs_error": float(abs(v - ref))} for k, v in routes.items()}
for k, v in out.items(): print(f"{k:62s} p = {v['p']:.12f}   |error| = {v['abs_error']:.2e}")
json.dump(out, open(sys.argv[2], "w"), indent=2)
