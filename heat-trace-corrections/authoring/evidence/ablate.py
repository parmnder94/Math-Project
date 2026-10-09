"""Wrong-but-plausible routes at the shipped parameters, versus the reference values.
Usage: python3 ablate.py <solution_dir> <out.json>
Uses the reference machinery (solution/solve.py); each route changes one modelling or numerical choice
(items to length 10, pairs to total length 10, except where a route says otherwise)."""
import importlib.util, json, os, sys
import numpy as np

spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py")); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
ref = json.load(open(sys.argv[3])) if len(sys.argv) > 3 else None
L = 10
TA = None; trce = 0.0
for k in range(S.NODES):
    z = S.RADIUS * np.exp(2j * np.pi * (k + 0.5) / S.NODES); pt = S.Point(z); w = np.exp(-S.BETA * z) * z / S.NODES
    S.LMAX = S.LNMAX = L
    T = S.class_table(pt) * w; TA = T if TA is None else TA + T; trce = trce + pt.trce * w
TA = TA.real; trce = float(trce.real)
words = [tuple(int(x) for x in r) for l in range(1, L + 1) for r in S.class_words(l)]
items = [(u * d, float(TA[i, d]), d) for i, u in enumerate(words) for d in range(1, L // len(u) + 1)]
q0 = {w: S.chi0_quotients(w) for w, _, _ in items}
def parts1(w):
    q1 = sum(1 for V, Ea, Eb, _ in S.quotients([w], max_excess=1) if Ea + Eb - V == 1)
    return q1, sum(a * b for a, b in q0[w])
def parts2(w1, w2):
    c0 = conn1 = conn0ab = 0
    for V, Ea, Eb, _ in S.quotients([w1, w2], max_excess=1, connected_only=True):
        if Ea + Eb - V == 0: c0 += 1; conn0ab += Ea * Eb
        else: conn1 += 1
    cross = sum(x1 * y2 + y1 * x2 for (x1, y1) in q0[w1] for (x2, y2) in q0[w2])
    return c0, conn1, conn0ab, cross
P1 = {w: parts1(w) for w, _, _ in items}
P2 = {}
for i, (w1, T1, d1) in enumerate(items):
    for j in range(i, len(items)):
        w2, T2, d2 = items[j]
        if len(w1) + len(w2) <= L: P2[(i, j)] = parts2(w1, w2)
def route(a1=lambda q1, ab, d: q1 - ab, c1=lambda c0, n1, n0ab, cr: n1 - n0ab - cr):
    m = sum(T * (len(q0[w]) - 1) for w, T, d in items) - trce
    m1 = sum(T * a1(*P1[w], d) for w, T, d in items)
    v = v1 = 0.0
    for (i, j), p in P2.items():
        f = (1 if i == j else 2) * items[i][1] * items[j][1]
        v += f * p[0]; v1 += f * c1(*p)
    return {"z_inf": trce / 2, "m": m, "m1": m1, "v": v, "v1": v1}
R = {"truncation: items and pairs to length 10 (reference machinery)": route(),
     "1/N of E fix: chi = -1 quotient count only (no -E_a E_b from the falling factorials)": route(a1=lambda q1, ab, d: q1),
     "1/N of E fix: proper powers given no 1/N term": route(a1=lambda q1, ab, d: (q1 - ab) if d == 1 else 0),
     "1/N of Cov: no disconnected cross term": route(c1=lambda c0, n1, n0ab, cr: n1 - n0ab),
     "1/N of Cov: no -E_a E_b on connected chi = 0 quotients": route(c1=lambda c0, n1, n0ab, cr: n1 - cr),
     "leading order only (m1 = v1 = 0)": dict(route(), m1=0.0, v1=0.0)}
if ref is None: ref = R["truncation: items and pairs to length 10 (reference machinery)"]
out = {}
for name, r in R.items():
    out[name] = {k: {"value": r[k], "rel_err": abs(r[k] - ref[k]) / abs(ref[k])} for k in ref}
    print(f"{name:85s}", " ".join(f"{k}={out[name][k]['rel_err']:.1e}" for k in ref), flush=True)
json.dump(out, open(sys.argv[2], "w"), indent=2)
