"""Wrong-but-plausible routes, each changing one modelling choice in the exact reference computation.
Usage: python3 ablate.py <solution_dir> <out.json>"""
import importlib.util, json, os, sys
import sympy as sp
spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py")); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
half = sp.Rational(3, 2)
routes = {
    "reference": S.delta(),
    "time reversal also flips the field (or: uniform density taken as reversible), I0 = 0": S.delta(field_in_reversal=False),
    "no entropy anomaly: heat rate equals the overdamped path irreversibility": (S.delta()[1], S.delta()[1], sp.Integer(0)),
    "Lorentz force dropped from the fast velocity process (kept in the current)": (S.delta(b=0)[0], S.delta()[1], S.delta(b=0)[0] - S.delta()[1]),
    "isotropic drag Gamma -> (3/2) I": S.delta(g1=half, g2=half),
}
out = {k: {"lim_Sdot": str(v[0]), "I0": str(v[1]), "Delta": str(v[2])} for k, v in routes.items()}
for k, v in out.items(): print(f"{k:85s} Delta = {v['Delta']}")
json.dump(out, open(sys.argv[2], "w"), indent=2)
