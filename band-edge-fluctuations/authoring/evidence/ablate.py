"""Wrong-but-plausible routes at the shipped parameters, versus the reference values.
Usage: python3 ablate.py <solution_dir> <out.json>
Uses the reference machinery (solution/solve.py); each route changes one modelling or numerical choice."""
import json, sys, itertools, importlib.util
from math import gcd
import numpy as np
spec = importlib.util.spec_from_file_location("solve", sys.argv[1] + "/solve.py"); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
z = complex(S.X0, S.ETA); PT = {'+': S.Point(z), '-': S.Point(np.conj(z))}; A = PT['+']
trce = float(-np.imag(A.trce)); TA = S.class_table(A); IT = TA.imag; LN = S.LNMAX
lens = []
for l in range(1, S.LMAX + 1):
    for u in S.cyc_reduced(l):
        if not S.is_power(u) and min(S.rotations(u) + S.rotations(S.winv(u))) == u: lens.append(l)
lens = np.array(lens)
memo = {}
def Pi(ns, signs, exact_max):
    key = (ns, signs, exact_max)
    if key not in memo:
        if sum(ns) > exact_max:
            memo[key] = np.sum(np.prod([TA[:, k] if s == '+' else np.conj(TA[:, k]) for k, s in zip(ns, signs)], axis=0))
        else:
            psi = sum(S.psi_marked([(PT[signs[0]], False, ns[0])] + [(PT[s], iv, k) for k, s, iv in zip(ns[1:], signs[1:], invs)])
                      for invs in itertools.product([False, True], repeat=len(ns) - 1))
            corr = sum(q ** (len(ns) - 1) * Pi(tuple(q * k for k in ns), signs, exact_max) for q in range(2, LN + 1) if q * max(ns) <= LN)
            memo[key] = psi - corr
    return memo[key]
def ImIm(n1, n2, em, mask=None):
    if n1 + n2 > em: return float(np.sum((IT[:, n1] * IT[:, n2])[mask] if mask is not None else IT[:, n1] * IT[:, n2]))
    return 0.5 * np.real(Pi((n1, n2), ('+', '-'), em) - Pi((n1, n2), ('+', '+'), em))
def ImImIm(ns, em, mask=None):
    if sum(ns) > em:
        p = IT[:, ns[0]] * IT[:, ns[1]] * IT[:, ns[2]]
        return float(np.sum(p[mask] if mask is not None else p))
    return float(np.real(sum((-1) ** s.count('-') * Pi(ns, s, em) for s in itertools.product('+-', repeat=3)) / (-8j)))
Md_exact = {d: S.M_d(A, d) for d in range(2, 5)}
def Md(d, em, mask=None):
    if d <= em: return Md_exact[d]
    cols = TA[:, [q * d for q in range(1, LN // d + 1)]]
    return np.sum(cols[mask] if mask is not None else cols)
def route(em=4, wv=lambda g: sum(S.divisors(g)), wk=lambda g: sum(d * d for d in S.divisors(g)), unit_mean=False, mask=None):
    if unit_mean:
        mob = lambda q: 0 if any(q % (p * p) == 0 for p in range(2, q + 1)) else (-1) ** sum(1 for p in range(2, q + 1) if q % p == 0 and all(p % r for r in range(2, p)))
        Sn = lambda n: sum(mob(q) * Md(q * n, em, mask) for q in range(1, LN // n + 1))
        m = sum(-np.imag(Sn(n)) for n in range(2, LN + 1)) - trce
    else:
        m = sum(-np.imag(Md(d, em, mask)) for d in range(2, LN + 1)) - trce
    v = sum(wv(gcd(a, b)) * ImIm(a, b, em, mask) for a in range(1, LN) for b in range(1, LN - a + 1))
    k = sum(wk(gcd(gcd(a, b), c)) * -ImImIm((a, b, c), em, mask) for a in range(1, LN) for b in range(1, LN - a + 1) for c in range(1, LN - a - b + 1))
    return {"z_inf": trce / 2, "mean_correction": float(m), "variance": float(v), "third_cumulant": float(k)}
ref = route()
routes = {"reference": ref,
          "gaussian fluctuations (k = 0)": dict(ref, third_cumulant=0.0),
          "power cycles ignored (j = 1 only)": route(wv=lambda g: 1, wk=lambda g: 1),
          "E fix(u^d) = 2 for every proper power": route(unit_mean=True),
          "E X_N = 2(N-1) z + m convention": dict(ref, mean_correction=ref["mean_correction"] + 2 * ref["z_inf"]),
          "exact resummation only to total power 3": route(em=3)}
for L in (6, 8, 10):
    routes[f"no resummation: primitive classes up to length {L}"] = route(em=0, mask=lens <= L)
out = {}
for name, r in routes.items():
    out[name] = {f: {"value": r[f], "rel_err": abs(r[f] - ref[f]) / abs(ref[f])} for f in ref}
    print(f"{name:52s} " + " ".join(f"{f}={r[f]:+.10f}({abs(r[f]-ref[f])/abs(ref[f]):.1e})" for f in ref), flush=True)
json.dump(out, open(sys.argv[2], "w"), indent=2)
