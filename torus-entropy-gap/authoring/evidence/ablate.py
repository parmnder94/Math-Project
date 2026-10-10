"""Plausible wrong routes for torus-entropy-gap, each changing one step of the reference (exact arithmetic).
Usage: python3 ablate.py <solution_dir> <out.json>"""
import importlib.util, json, os, sys
import sympy as sp
spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
R_ = sp.Rational
base = dict(a=R_(4, 5), lam=R_(1, 2), gam=(1, 3), tau=(1, 2), b=R_(1, 2), k=2, alpha=R_(1, 3), beta=R_(1, 2))
m = S.Model(**base)
ref = S.analyse(m)
routes = {"reference": (ref["heat"], ref["sigma"])}
r = S.analyse(m, kl_velocity="jkin"); routes["KL rate from the kinetic current rho h instead of j"] = (ref["heat"], r["sigma"])
# x treated as a one-dimensional Markov diffusion with y-averaged coefficients: rate = J^2 int dx / (rho_x Dbar)
jx_avg = S.torus_average(m, ref["j"][0]) / S.torus_average(m, m.rho)          # = J / (2 pi * normalised) up to norm
# exact: rho_x Dbar = T^2/<T^2> * T * tr(D0)/2 / (2 pi);  J = 2 pi <j_x>/<T^2> / (4 pi^2) ... use normalised torus units
MC0 = [[sum(m.M0[i][l]*m.C0[l][j] for l in range(2)) for j in range(2)] for i in range(2)]
trD0 = sp.Rational(int((MC0[0][0] + MC0[1][1]).numerator), int((MC0[0][0] + MC0[1][1]).denominator))
a = R_(4, 5); T2 = S.Tmoment(a, 2)
jbar = S.torus_average(m, ref["j"][0])/T2          # <j_x>/<T^2> with rho normalised to mean 1 on the torus
# with rho normalised to mean 1, the 1D rate is jbar^2 * < 1/(rho_x Dbar) >, rho_x = T^2/<T^2>, Dbar = T trD0/2
sig1d = jbar**2*T2*S.Tmoment(a, -3)/(trD0/2)
routes["x treated as a 1D Markov diffusion with y-averaged coefficients"] = (ref["heat"], sp.nsimplify(sig1d))
r = S.analyse(S.Model(**base, rho_pow=1)); routes["stationary density taken proportional to T instead of T^2"] = (r["heat"], r["sigma"])
routes["finite heat without the current part B h of the force"] = (S.analyse(m, heat_force=[m.F[i] - m.Bh[i] for i in range(2)])["heat"], ref["sigma"])
r = S.analyse(S.Model(**base, wtau=(1, 1))); routes["entropy weight gam_a/T for both reservoirs (T_2 = 2T ignored)"] = (r["heat"], ref["sigma"])
routes["observed x path taken as reversible (KL rate 0)"] = (ref["heat"], sp.Integer(0))
out = {}
for k, (E, I) in routes.items():
    D = sp.nsimplify(E - I)
    out[k] = {"finite_heat": str(E), "kl_rate": str(I), "Delta": str(D), "Delta_float": float(D),
              "abs_error": float(abs(D - ref["Delta"]))}
    print(f"{k:70s} Delta = {float(D):+.10f}  (error {out[k]['abs_error']:.3e})")
g = S.analyse(S.Model(a=R_(3, 5), lam=R_(2, 3), gam=(1, 2), tau=(1, 2), b=1, k=1, alpha=R_(2, 5), beta=R_(2, 3)))
out["original problem (user-supplied parameters) reproduced by the same pipeline"] = {
    "finite_heat": str(g["heat"]), "kl_rate": str(g["sigma"]), "Delta": str(g["Delta"]),
    "independent_derivation_value": "-141231719/146966400", "match": str(g["Delta"]) == "-141231719/146966400"}
print("original problem:", g["Delta"], g["heat"], g["sigma"])
json.dump(out, open(sys.argv[2], "w"), indent=2)
