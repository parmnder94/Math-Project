"""Ground truth for heat-trace-fluctuations, with independent cross-checks.

A. Primary: class sums over primitive conjugacy classes (mod inversion) up to length LMAX_TRUTH, with the
   conjugator sums resummed exactly (solution/solve.py machinery, deeper truncation than the solver).
B. Independent check of every class sum T_kappa(n) for |u|<=3: cactus-geodesic coefficients (no
   linearization) with explicit enumeration of conjugators g, |g| <= GMAX.
C. Independent word-level enumeration with cactus coefficients, |w| <= LW (no class resummation).
Usage: python3 truth.py <solution_dir> <out_dir>
"""
import json, sys, time, importlib.util
import numpy as np
sys.path.insert(0, __file__.rsplit("/", 1)[0])
from cactus_coef import CactusCoef
spec = importlib.util.spec_from_file_location("solve", sys.argv[1] + "/solve.py"); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
OUT = sys.argv[2]; LMAX_TRUTH, LNMAX_TRUTH = 13, 26; GMAX = 9; LW = 10
t0 = time.time(); report = {}

# ---------- A
wq, Phi, G = S.contour_data(S.Tree()); tr_ce = np.einsum('n,nii->', wq, G[:, :2, :2]).real; Om = S.omegas(Phi, G)
msum = v = k3 = 0.0; per_len = {}; Tstore = {}
for u in S.primitive_classes(LMAX_TRUTH):
    nmax = max(1, LNMAX_TRUTH // len(u)); T = S.class_T(u, nmax, S_wq := wq, Phi, Om)
    if len(u) <= 3: Tstore[u] = T
    dm = sum((S.ndiv(n) - 1) * T[n] for n in range(2, nmax + 1)); dv = dk = 0.0
    for j in range(1, nmax + 1):
        Sj = sum(T[j * q] for q in range(1, nmax // j + 1)); dv += j * Sj ** 2; dk += j * j * Sj ** 3
    msum += dm; v += dv; k3 += dk
    pl = per_len.setdefault(len(u), [0.0, 0.0, 0.0]); pl[0] += dm; pl[1] += dv; pl[2] += dk
truth = {"z_inf": tr_ce / 2, "mean_correction": msum - tr_ce, "variance": v, "third_cumulant": k3}
report["A_class_sums"] = {"lmax": LMAX_TRUTH, "values": truth,
                          "per_class_length": {l: {"m": x[0], "v": x[1], "k": x[2]} for l, x in sorted(per_len.items())},
                          "seconds": time.time() - t0}
print("A", truth, flush=True)

# ---------- B
cc = CactusCoef(S.BFIELD, S.LAM, S.BETA)
report["B_zinf_cactus"] = float(np.einsum('n,nii->', cc.w, cc.G).real / 2)
INV = S.INV
def conjugators(first, last, gmax):
    """reduced g = g_1..g_k (k<=gmax) with g y g^-1 reduced for y with given first/last letter"""
    out = [()]; frontier = [()]
    for k in range(1, gmax + 1):
        new = []
        for g in frontier:
            for x in range(4):
                if k == 1:
                    if x == INV[first] or x == last: continue
                    new.append((x,))
                else:
                    if x == INV[g[0]]: continue
                    new.append((x,) + g)
        frontier = new; out += new
    return out
worst = 0.0; rows = []
for u, T in Tstore.items():
    for n in range(1, min(len(T) - 1, 4) + 1):
        tot = 0.0
        for word in (u, S.winv(u)):
            for r in S.rotations(word):
                y = r * n
                for g in conjugators(y[0], y[-1], GMAX):
                    tot += cc.trace_coef(g + y + S.winv(g)).real
        err = abs(tot - T[n]); worst = max(worst, err)
        rows.append({"u": "".join("aAbB"[x] for x in u), "n": n, "omega_resummed": T[n], "cactus_enumerated": tot, "abs_diff": err})
report["B_class_sum_check"] = {"gmax": GMAX, "worst_abs_diff": worst, "rows": rows}
print("B worst", worst, time.time() - t0, flush=True)

# ---------- C
def cycred(w):
    i, j = 0, len(w) - 1
    while i < j and w[i] == INV[w[j]]: i += 1; j -= 1
    return w[i:j + 1]
def root(w):
    for p in range(1, len(w) + 1):
        if len(w) % p == 0 and w == w[:p] * (len(w) // p): return w[:p], len(w) // p
def canon(w): return min(S.rotations(w) + S.rotations(S.winv(w)))
conv = []
Sacc = {}; Macc = 0.0; frontier = [()]
for L in range(1, LW + 1):
    new = []
    for w in frontier:
        for x in range(4):
            if w and x == INV[w[0]]: continue
            new.append((x,) + w)
    frontier = new
    for w in frontier:
        t = cc.trace_coef(w).real; u, d = root(cycred(w)); kap = canon(u)
        Macc += (S.ndiv(d) - 1) * t
        for j in range(1, d + 1):
            if d % j == 0: Sacc[(kap, j)] = Sacc.get((kap, j), 0.0) + t
    if L >= 6:
        vv = sum(j * s * s for (k, j), s in Sacc.items()); kk = sum(j * j * s ** 3 for (k, j), s in Sacc.items())
        conv.append({"max_word_length": L, "mean_correction": Macc - tr_ce, "variance": vv, "third_cumulant": kk})
        print("C", conv[-1], time.time() - t0, flush=True)
report["C_word_enumeration"] = conv
json.dump(truth, open(OUT + "/truth.json", "w"), indent=2)
json.dump(report, open(OUT + "/truth_report.json", "w"), indent=2, default=str)
print("done", time.time() - t0)
