"""Plausible wrong routes, each changing one step of the reference computation (exact arithmetic).
Usage: python3 ablate.py <solution_dir> <out.json>"""
import importlib.util, json, os, sys
import sympy as sp
spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
bq, bx = sp.Rational(23539, 2160), sp.Rational(8447, 540)
routes = {
    "reference": S.probability(bq, bx),
    "endpoint frame D_r(1) set to the identity": S.probability(bq, bx, use_D=False),
    "slow connection omitted": S.probability(bq, bx, connection=0),
    "connection from R_r only (-pi K_r, D_r frame term dropped)": S.probability(bq, bx, connection=1),
    "reversed pulse keeps the carrier phase (no phi -> -phi)": S.probability(bq, bx, flip_phase=False),
    "carrier phases ignored in the effective coupling": S.probability(bq, bx, carrier=False),
}
out = {k: {"p": (str(v[0]) if len(str(v[0])) < 60 else None), "decimal": float(v[0])} for k, v in routes.items()}
for k, v in out.items(): print(f"{k:62s} {v['decimal']:.12f}  {v['p']}")
json.dump(out, open(sys.argv[2], "w"), indent=2)
