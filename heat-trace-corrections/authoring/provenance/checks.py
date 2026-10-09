"""Independent checks of the reference machinery (solution/solve.py).
Usage: python3 checks.py <out.json>
(1) Heat kernel: tr c_e from the contour integral of the tree resolvent against the exact moment series
    sum_n (-beta)^n tau(P^n) / n!, with tau(P^n) from walks on reduced words of F_2.
(2) Word maps: the quotient expansion E[prod fix_{w_i}](N) = sum_G (N)_V / ((N)_{E_a} (N)_{E_b}) (only
    quotients with V <= N) against exact averages over all of S_N x S_N, single words and pairs, N = 5.
(3) Finite N: E X_N and Var X_N computed exactly by brute force over S_N x S_N (sigma over conjugacy-class
    representatives, weighted) against the class sums combined with the finite-N quotient expansion, N = 4..7.
(4) Coefficient extraction: m1 and v1 from the pruned chi >= -1 enumeration against Richardson extrapolation
    of the full (unpruned) quotient expansion evaluated in exact rational arithmetic at N = 2e4 ... 1.6e5,
    with the same class truncation (words and pairs up to length 8)."""
import importlib.util, itertools, json, os, sys, time
from collections import Counter
from fractions import Fraction
from math import factorial
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("solve", os.path.join(HERE, "..", "..", "solution", "solve.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
INV = S.INV; rep = {}
Ta, Tb = S.hop(S.sx), S.hop(S.sy)


def class_sums(L):
    S.LMAX = S.LNMAX = L; TA = None; trce = 0.0
    for k in range(S.NODES):
        z = S.RADIUS * np.exp(2j * np.pi * (k + 0.5) / S.NODES); pt = S.Point(z); w = np.exp(-S.BETA * z) * z / S.NODES
        T = S.class_table(pt) * w; TA = T if TA is None else TA + T; trce = trce + pt.trce * w
    words = [tuple(int(x) for x in r) for l in range(1, L + 1) for r in S.class_words(l)]
    items = [(u * d, float(TA.real[i, d])) for i, u in enumerate(words) for d in range(1, L // len(u) + 1)]
    return float(trce.real), items


def ff(N, k):
    r = 1
    for i in range(k): r *= N - i
    return r


def expect(ws, N, exact=False):
    tot = Fraction(0) if exact else 0.0
    for V, Ea, Eb, _ in S.quotients(list(ws), max_V=N if N < 10 ** 3 else None):
        tot += Fraction(ff(N, V), ff(N, Ea) * ff(N, Eb)) if exact else ff(N, V) / (ff(N, Ea) * ff(N, Eb))
    return tot


# (1) heat kernel against exact moments
terms = [(S.BFIELD * S.sz, ()), (Ta, (0,)), (Ta.conj().T, (1,)), (Tb, (2,)), (Tb.conj().T, (3,))]
def mul(g, w):
    out = list(w)
    for x in reversed(g):
        if out and out[0] == INV[x]: out.pop(0)
        else: out.insert(0, x)
    return tuple(out)
NM = 26; v = {(): np.eye(2, dtype=complex)}; mom = []
for n in range(NM + 1):
    mom.append(np.trace(v.get((), np.zeros((2, 2)))))
    nv = {}
    for w, cw in v.items():
        for T, g in terms:
            k = mul(g, w); nv[k] = nv.get(k, 0) + T @ cw
    v = {k: c for k, c in nv.items() if len(k) <= NM - n}
series = float(sum((-S.BETA) ** n / factorial(n) * mom[n] for n in range(NM + 1)).real)
trce, _ = class_sums(1)
rep["heat_kernel_trace"] = {"contour": trce, "moment_series": series, "abs_diff": abs(trce - series)}
print("(1)", rep["heat_kernel_trace"], flush=True)

# (2) word-map expansion against exact averages over S_5 x S_5
def brute_words(ws, N):
    perms = [np.array(p) for p in itertools.permutations(range(N))]; tot = 0
    for s in perms:
        for t in perms:
            P = {0: s, 1: np.argsort(s), 2: t, 3: np.argsort(t)}; prod = 1
            for w in ws:
                x = np.arange(N)
                for c in reversed(w): x = P[c][x]
                prod *= int(np.sum(x == np.arange(N)))
            tot += prod
    return Fraction(tot, len(perms) ** 2)
cases = [[(0, 2, 1, 3)], [(0, 0, 2, 2)], [(0, 2, 0, 3)], [(0, 0)], [(0, 2, 0, 2)], [(0, 0, 2, 1, 3)], [(0,), (2,)], [(0,), (0,)],
         [(0, 2), (0, 2)], [(0, 2), (3, 1)], [(0, 2, 1, 3), (0,)], [(0, 0, 2), (2, 1)], [(0, 2, 1, 3), (0, 2, 1, 3)],
         [(0, 0, 2, 2), (0, 2)], [(0, 2, 1, 3), (2, 0, 3, 1)], [(0, 0, 0), (0, 2, 2)]]
res2 = []
for ws in cases:
    a, b = expect(ws, 5, exact=True), brute_words(ws, 5)
    res2.append({"words": str(ws), "expansion": str(a), "brute_force": str(b), "equal": a == b})
rep["word_maps_N5"] = res2; rep["word_maps_all_equal"] = all(r["equal"] for r in res2)
print("(2) all equal:", rep["word_maps_all_equal"], flush=True)

# (3) exact finite-N moments
hop_triv = S.BFIELD * S.sz + Ta + Ta.conj().T + Tb + Tb.conj().T
Xtriv = float(np.sum(np.exp(-S.BETA * np.linalg.eigvalsh(hop_triv))))
def brute_moments(N):
    def parts(n, mx=None):
        if n == 0: yield []; return
        for k in range(min(n, mx or n), 0, -1):
            for p in parts(n - k, k): yield [k] + p
    taus = [np.array(t) for t in itertools.permutations(range(N))]; s1 = s2 = W = 0.0
    for lamb in parts(N):
        sig = []; st = 0
        for k in lamb: sig += [st + (i + 1) % k for i in range(k)]; st += k
        sig = np.array(sig); size = factorial(N)
        for k, c in Counter(lamb).items(): size //= (k ** c) * factorial(c)
        Pa = np.zeros((N, N)); Pa[sig, np.arange(N)] = 1; Hs = []
        for t in taus:
            Pb = np.zeros((N, N)); Pb[t, np.arange(N)] = 1
            Hs.append(np.kron(S.BFIELD * S.sz, np.eye(N)) + np.kron(Ta, Pa) + np.kron(Ta.conj().T, Pa.T) + np.kron(Tb, Pb) + np.kron(Tb.conj().T, Pb.T))
        X = np.exp(-S.BETA * np.linalg.eigvalsh(np.array(Hs))).sum(axis=1) - Xtriv
        s1 += size * X.sum(); s2 += size * (X ** 2).sum(); W += size * len(taus)
    return s1 / W, s2 / W - (s1 / W) ** 2
trce10, items10 = class_sums(10)
res3 = []
for N in (4, 5, 6, 7):
    bm, bv = brute_moments(N)
    E1 = {w: expect([w], N) for w, _ in items10}
    fm = (N - 1) * trce10 + sum(T * (E1[w] - 1) for w, T in items10)
    fv = 0.0
    for i, (w1, T1) in enumerate(items10):
        for j in range(i, len(items10)):
            w2, T2 = items10[j]
            if len(w1) + len(w2) > 10: continue
            fv += (1 if i == j else 2) * T1 * T2 * (expect([w1, w2], N) - E1[w1] * E1[w2])
    res3.append({"N": N, "brute_mean": bm, "formula_mean": fm, "brute_var": bv, "formula_var": fv,
                 "mean_rel_diff": abs(bm - fm) / abs(bm), "var_rel_diff": abs(bv - fv) / abs(bv)})
    print("(3)", res3[-1], flush=True)
rep["finite_N_exact"] = res3

# (4) Richardson extraction of m1, v1 from the full expansion
trce8, items8 = class_sums(8)
pairs = [(i, j) for i in range(len(items8)) for j in range(i, len(items8)) if len(items8[i][0]) + len(items8[j][0]) <= 8]
Q1 = {w: S.quotients([w]) for w, _ in items8}; Q2 = {p: S.quotients([items8[p[0]][0], items8[p[1]][0]]) for p in pairs}
def stats(N):
    E1 = {w: sum(Fraction(ff(N, V), ff(N, Ea) * ff(N, Eb)) for V, Ea, Eb, _ in Q1[w]) for w, _ in items8}
    mc = sum(T * float(E1[w] - 1) for w, T in items8); var = 0.0
    for (i, j) in pairs:
        w1, T1 = items8[i]; w2, T2 = items8[j]
        cov = sum(Fraction(ff(N, V), ff(N, Ea) * ff(N, Eb)) for V, Ea, Eb, _ in Q2[(i, j)]) - E1[w1] * E1[w2]
        var += (1 if i == j else 2) * T1 * T2 * float(cov)
    return mc, var
Ns = [20000, 40000, 80000, 160000]; vals = [stats(N) for N in Ns]
A = np.array([[1, 1 / N, 1 / N ** 2, 1 / N ** 3] for N in Ns])
cm = np.linalg.solve(A, np.array([x[0] for x in vals])); cv = np.linalg.solve(A, np.array([x[1] for x in vals]))
out8 = os.path.join(os.path.dirname(os.path.abspath(sys.argv[1])), "pruned_8.json")
import subprocess
subprocess.run([sys.executable, os.path.join(HERE, "..", "..", "solution", "solve.py")], check=True, stdout=subprocess.DEVNULL,
               env=dict(os.environ, LMAX="8", LPAIR="8", ANSWER_PATH=out8))
pruned = json.load(open(out8)); os.remove(out8)
rep["coefficient_extraction"] = {"truncation": "classes and pairs up to length 8",
                                 "richardson_m1": cm[1], "pruned_m1": pruned["m1"], "m1_rel_diff": abs(cm[1] - pruned["m1"]) / abs(pruned["m1"]),
                                 "richardson_v1": cv[1], "pruned_v1": pruned["v1"], "v1_rel_diff": abs(cv[1] - pruned["v1"]) / abs(pruned["v1"])}
print("(4)", rep["coefficient_extraction"], flush=True)
json.dump(rep, open(sys.argv[1], "w"), indent=2, default=str)
