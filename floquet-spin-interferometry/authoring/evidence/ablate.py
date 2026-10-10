"""Plausible wrong routes, each changing one step of the reference computation (60-digit arithmetic).
Usage: python3 ablate.py <solution_dir> <out.json>   (the c2, c3 routes need python-flint)"""
import importlib.util, json, os, sys
import mpmath as mp
import sympy as sp
sys.path.insert(0, sys.argv[1])
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

# routes for c2 and c3 (reference: Taylor coefficients of the smooth interpolant F, solution/smooth.py)
ref_c2, ref_c3, _ = S.taylor_coefficients(ref_p, ref_c)
def coeffs_of(f, h, m=8):
    pts = [(k*h, f(k*h)) for k in list(range(-m, 0)) + list(range(1, m + 1))]
    return mp.lu_solve(mp.matrix([[e**i for i in range(2*m)] for e, _ in pts]), mp.matrix([v for _, v in pts]))
mp.mp.dps = 60
kick = coeffs_of(lambda e: P(eps=e), mp.mpf("1e-3"))
_, _, und = S.taylor_coefficients(ref_p, ref_c, dressed=False)
c_routes = {
    "reference (smooth interpolant F)": (ref_c2, ref_c3),
    "static boundary kicks around exp(-i G_r) (first-order K1 exponentiated, B_8 only)": (kick[2], kick[3]),
    "exact one-period Floquet Hamiltonian, slow connection not dressed by the micromotion": (und[2], und[3]),
    "remainder (p_N - p_inf - c1 eps_N)/eps_N^2 at N = 1e40 taken as c2": (ref_c2 + ref_c3*mp.mpf("6.03e-6"), mp.nan),
}
for k, (c2, c3) in c_routes.items():
    out[k] = {"c2": float(c2), "abs_error_c2": float(abs(c2 - ref_c2))}
    if not mp.isnan(c3): out[k].update({"c3": float(c3), "abs_error_c3": float(abs(c3 - ref_c3))})
    print(f"{k:62s} c2 = {float(c2):+.10f} ({out[k]['abs_error_c2']:.1e})   c3 = {float(c3):+.10f}")
out["exact one-period Floquet Hamiltonian, slow connection not dressed by the micromotion"].update(
    {"p": float(und[0]), "c1": float(und[1]), "abs_error_c1": float(abs(und[1] - ref_c))})
json.dump(out, open(sys.argv[2], "w"), indent=2)
