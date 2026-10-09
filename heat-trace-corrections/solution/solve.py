"""Reference solution: 1/N corrections to the heat-trace statistics of spin-orbit random 4-regular graphs.

Model. H_N = P(U_a, U_b) with P = B sz + Ta l(a) + Tb l(b) + h.c. in M_2(C[F_2]), U_a, U_b uniform random
permutation matrices, X_N = Tr[exp(-beta H_N) Q_N].

1. Exact trace expansion (every N):  X_N = (N-1) tr c_e + sum_{w != e} tr(c_w) (fix(w(sigma,tau)) - 1),
   c_w = <delta_w, exp(-beta P) delta_e>.  fix(w) depends only on the conjugacy class y of w, so all sums run
   over classes with weights F(y) = sum of tr c_w over the class, and T_u(d) = F(u^d) + F(u^-d).
2. Class sums.  On the 4-regular Cayley tree G_z(w,e) = Phi_{s_1}...Phi_{s_k} G_z(e,e) with Phi_s = H_s A_s
   (branch resolvents H_s by Newton's method); conjugator sums are resummed exactly, F_z(y) = tr(Om Phi_y).
   The heat kernel is the contour integral exp(-beta P) = (1/2 pi i) oint exp(-beta z) (z - P)^-1 dz
   (trapezoid rule on |z| = R, exponentially convergent).
3. Word maps.  E[prod_i fix_{w_i}](N) = sum over quotients G of the disjoint union of the cycles C_{w_i}
   (folded labelled graphs, Linial-Puder) of (N)_V / ((N)_{E_a} (N)_{E_b}), and (N)_k = N^k (1 - k(k-1)/2N + ...).
   Hence, with chi = V - E:
     E fix_w      = tau(d) + a1(w)/N + O(N^-2),
        a1(w)     = #{quotients with chi = -1} - sum_{chi = 0 quotients} E_a E_b          (w = u^d)
     Cov(fix_w1, fix_w2) = C0 + C1/N + O(N^-2),
        C0        = #{connected quotients with chi = 0}                                   (= sigma(gcd) or 0)
        C1        = #{connected, chi = -1} - sum_{connected, chi = 0} E_a E_b
                    - sum_{G1 in Q0(w1), G2 in Q0(w2)} (E_a(G1) E_b(G2) + E_b(G1) E_a(G2)),
   the last term coming from (N)_{p+q} / ((N)_p (N)_q) = 1 - pq/N + ... for disconnected quotients.
4. Constants:
     z_inf = tr c_e / 2,  m = sum_y F(y)(tau(d) - 1) - tr c_e,  m1 = sum_y F(y) a1(y),
     v = sum_{y1,y2} F F C0,  v1 = sum_{y1,y2} F F C1.
   Class contributions fall about 10x per unit of word length; classes up to length 11 and pairs of classes up
   to total length 12 give the constants to about 1e-8 relative.  The quotient counts are invariant under the
   eight signed permutations of the letters, so each orbit of words (or pairs) is enumerated once.
"""
import json, os, time
import numpy as np

LAM, BFIELD, BETA = 0.6, 1.5, 0.3
LMAX = int(os.environ.get("LMAX", 11)); LNMAX = LMAX          # class words (and powers) up to this length
LPAIR = int(os.environ.get("LPAIR", 12))                       # pairs of classes up to this total length
RADIUS, NODES = 12.0, 48
OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")
INV = [1, 0, 3, 2]                                  # letters a, a^-1, b, b^-1
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def hop(s): return np.cos(LAM) * I2 + 1j * np.sin(LAM) * s


def linearization():
    """P is already linear in l(a), l(b): K = 2, coefficient A_s of l(s)."""
    Ta, Tb = hop(sx), hop(sy)
    return 2, BFIELD * sz, [Ta, Ta.conj().T, Tb, Tb.conj().T]


class Point:
    """Resolvent data at a non-real z: Phi_s, G(e,e) and Om[(f,l)] with F(y) = tr(Om Phi_y)."""
    def __init__(self, z):
        K, b0, A = linearization(); self.K = K
        Lam = np.zeros((K, K), complex)

        def F(H, w):
            Lam[0, 0] = Lam[1, 1] = w
            return [np.linalg.inv(Lam - b0 - sum(A[t].conj().T @ H[t] @ A[t] for t in range(4) if t != INV[a])) for a in range(4)]

        def newton(H, w):
            for _ in range(80):
                Fs = F(H, w); r = np.concatenate([(Fs[a] - H[a]).ravel() for a in range(4)])
                if np.abs(r).max() < 1e-14: return H
                n = 4 * K * K; J = -np.eye(n, dtype=complex)
                for t in range(4):
                    for j in range(K * K):
                        E = np.zeros(K * K, complex); E[j] = 1; X = A[t].conj().T @ E.reshape(K, K) @ A[t]
                        J[:, t * K * K + j] += np.concatenate([(Fs[a] @ X @ Fs[a]).ravel() if t != INV[a] else np.zeros(K * K) for a in range(4)])
                dv = np.linalg.solve(J, -r)
                H = [H[a] + dv[a * K * K:(a + 1) * K * K].reshape(K, K) for a in range(4)]
            raise RuntimeError("Newton did not converge")

        sgn = 1.0 if z.imag > 0 else -1.0
        w0 = complex(z.real, sgn * 0.5); H = [-1j * sgn * np.eye(K) for _ in range(4)]
        for _ in range(5000):
            Hn = F(H, w0); H = [0.5 * (H[a] + Hn[a]) for a in range(4)]
        for e in [0.5, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02, 0.01]:
            if abs(z.imag) < e: H = newton(H, complex(z.real, sgn * e))
        H = newton(H, z)
        Lam[0, 0] = Lam[1, 1] = z
        self.Phi = np.array([H[a] @ A[a] for a in range(4)])
        self.G = np.linalg.inv(Lam - b0 - sum(A[t].conj().T @ H[t] @ A[t] for t in range(4)))
        self.trce = np.trace(self.G[:2, :2])
        KK = K * K; Ak = [np.kron(self.Phi[a], self.Phi[INV[a]].T) for a in range(4)]
        T = np.zeros((4 * KK, 4 * KK), complex)
        for a in range(4):
            for b in range(4):
                if a != INV[b]: T[a * KK:(a + 1) * KK, b * KK:(b + 1) * KK] = Ak[a]
        assert max(abs(np.linalg.eigvals(T))) < 1, "conjugator series diverges"
        g = np.zeros((K, K), complex); g[:2, :] = self.G[:, :2].T; gv = g.ravel()
        adj = np.linalg.solve((np.eye(4 * KK) - T).T, np.concatenate([gv] * 4))
        self.Om = {}
        for f in range(4):
            for l in range(4):
                if l == INV[f]: continue
                om = gv.copy()
                for t in range(4):
                    if t != INV[f] and t != l: om = om + Ak[t].T @ adj[t * KK:(t + 1) * KK]
                self.Om[(f, l)] = om.reshape(K, K).T


def winv(w): return tuple(INV[x] for x in reversed(w))
def rotations(w): return [w[k:] + w[:k] for k in range(len(w))]
def is_power(w): return any(len(w) % p == 0 and w == w[:p] * (len(w) // p) for p in range(1, len(w)))
def divisors(n): return [d for d in range(1, n + 1) if n % d == 0]


def cyc_reduced(l):
    out = []
    def rec(w):
        if len(w) == l:
            if w[-1] != INV[w[0]]: out.append(tuple(w))
            return
        for x in range(4):
            if not w or x != INV[w[-1]]: rec(w + [x])
    rec([]); return out


def class_words(l):
    """Primitive cyclically reduced words of length l that are minimal among their rotations and the rotations
    of their inverse (one per class mod inversion), as rows of an int array, in lexicographic order."""
    inv = np.array(INV)
    W = np.arange(4)[:, None]
    for _ in range(l - 1):
        W = np.concatenate([np.column_stack([W, np.full(len(W), x)]) for x in range(4)])
        W = W[W[:, -1] != inv[W[:, -2]]]
    W = W[W[:, 0] != inv[W[:, -1]]]
    pw = 4 ** np.arange(l - 1, -1, -1)
    code = lambda A: A @ pw
    c0 = code(W); cmin = c0.copy(); prim = np.ones(len(W), bool)
    Wi = inv[W[:, ::-1]]
    for k in range(l):
        ck = code(np.roll(W, -k, axis=1)); cmin = np.minimum(cmin, np.minimum(ck, code(np.roll(Wi, -k, axis=1))))
        if k: prim &= ck != c0
    keep = prim & (c0 == cmin)
    return W[keep][np.argsort(c0[keep])]


def class_table(pt):
    """rows: primitive classes (mod rotation and inversion) up to LMAX; T[n] = sum over conjugates of u^{+-n} of tr c_w."""
    K = pt.K; Phi = np.array(pt.Phi); Om = np.zeros((4, 4, K, K), complex)
    for (f, l), M in pt.Om.items(): Om[f, l] = M
    inv = np.array(INV); rows = []
    for l in range(1, LMAX + 1):
        U = class_words(l); nmax = LNMAX // l
        T = np.zeros((len(U), LNMAX + 1), complex)
        for Wd in (U, inv[U[:, ::-1]]):
            for k in range(l):
                R = np.roll(Wd, -k, axis=1)
                M = Phi[R[:, 0]]
                for j in range(1, l): M = M @ Phi[R[:, j]]
                O = Om[R[:, 0], R[:, -1]]; Mn = M
                for n in range(1, nmax + 1):
                    if n > 1: Mn = Mn @ M
                    T[:, n] += np.einsum('cij,cji->c', O, Mn)
        rows.append(T)
    return np.concatenate(rows)



def quotients(words, max_excess=None, connected_only=False, max_V=None):
    """Yield (V, Ea, Eb, ncomp) for every quotient of the disjoint union of the cycles C_w.
    max_excess bounds E - V (i.e. -chi) of the final graph."""
    nc = len(words)
    out = []          # out[c][x] -> class reached from class c by letter x, or -1
    comp = []         # component id of each class (union-find lite: components only merge)
    state = {"V": 0, "E": 0, "Ea": 0}
    parent = []
    def find(c):
        while parent[c] != c: c = parent[c]
        return c
    results = []
    def new_class():
        out.append([-1, -1, -1, -1]); parent.append(len(parent)); state["V"] += 1; return len(out) - 1
    def add_edge(c, x, d):
        out[c][x] = d; out[d][INV[x]] = c; state["E"] += 1
        if x < 2: state["Ea"] += 1
    def del_edge(c, x, d):
        out[c][x] = -1; out[d][INV[x]] = -1; state["E"] -= 1
        if x < 2: state["Ea"] -= 1
    def excess(): return state["E"] - state["V"]
    def walk(k, i, cur, base, merged):
        w = words[k]; L = len(w)
        if max_excess is not None and excess() - (nc - 1 - k) > max_excess: return
        if i == L:
            if cur != base: return
            if k + 1 == nc:
                comps = len({find(c) for c in range(len(out))})
                if connected_only and comps > 1: return
                results.append((state["V"], state["Ea"], state["E"] - state["Ea"], comps)); return
            start_cycle(k + 1); return
        x = w[i]; nxt = out[cur][x]
        if nxt >= 0:
            walk(k, i + 1, nxt, base, merged); return
        last = (i == L - 1)
        if not last and (max_V is None or len(out) < max_V):
            d = new_class(); add_edge(cur, x, d); parent[d] = find(cur)
            walk(k, i + 1, d, base, merged)
            del_edge(cur, x, d); out.pop(); parent.pop(); state["V"] -= 1
        cands = [base] if last else range(len(out))
        for d in cands:
            if out[d][INV[x]] >= 0: continue
            add_edge(cur, x, d)
            rc, rd = find(cur), find(d); saved = None
            if rc != rd: saved = (rc, parent[rc]); parent[rc] = rd
            walk(k, i + 1, d, base, merged)
            if saved: parent[saved[0]] = saved[1]
            del_edge(cur, x, d)
    def start_cycle(k):
        # base vertex of cycle k: an existing class or a new one
        nV = len(out)
        for c in range(nV):
            walk(k, 0, c, c, True)
        if max_V is None or nV < max_V:
            c = new_class(); walk(k, 0, c, c, False); out.pop(); parent.pop(); state["V"] -= 1
    c0 = new_class(); walk(0, 0, c0, c0, False); out.pop(); parent.pop(); state["V"] -= 1
    return results


def tau(d): return len(divisors(d))


def a1_coef(w):
    s = 0
    for V, Ea, Eb, _ in quotients([w], max_excess=1):
        s += 1 if Ea + Eb - V == 1 else -Ea * Eb
    return s


def chi0_quotients(w):
    return [(Ea, Eb) for V, Ea, Eb, _ in quotients([w], max_excess=0)]


def cov_coefs(w1, w2, q1, q2):
    c0 = c1 = 0
    for V, Ea, Eb, _ in quotients([w1, w2], max_excess=1, connected_only=True):
        if Ea + Eb - V == 0: c0 += 1; c1 -= Ea * Eb
        else: c1 += 1
    for (x1, y1) in q1:
        for (x2, y2) in q2: c1 -= x1 * y2 + y1 * x2
    return c0, c1


SYM = [(0, 1, 2, 3), (1, 0, 2, 3), (0, 1, 3, 2), (1, 0, 3, 2), (2, 3, 0, 1), (3, 2, 0, 1), (2, 3, 1, 0), (3, 2, 1, 0)]


def canon(w): return min(rotations(w) + rotations(winv(w)))


def orbit_key(ws):
    return min(tuple(sorted(canon(tuple(g[x] for x in w)) for w in ws)) for g in SYM)


def main():
    t0 = time.time()
    TA = None; trce = 0.0
    for k in range(NODES):
        z = RADIUS * np.exp(2j * np.pi * (k + 0.5) / NODES)
        pt = Point(z); wgt = np.exp(-BETA * z) * z / NODES
        T = class_table(pt) * wgt
        TA = T if TA is None else TA + T
        trce = trce + pt.trce * wgt
    TA = TA.real; trce = float(trce.real)
    words = [tuple(int(x) for x in r) for l in range(1, LMAX + 1) for r in class_words(l)]
    items = [(u * d, float(TA[i, d])) for i, u in enumerate(words) for d in range(1, LMAX // len(u) + 1)]
    q0 = {w: chi0_quotients(w) for w, _ in items}
    m = sum(T * (len(q0[w]) - 1) for w, T in items) - trce          # len(q0[w]) = tau(d) for w = u^d
    cache = {}
    def cached(key, fn):
        if key not in cache: cache[key] = fn()
        return cache[key]
    m1 = sum(T * cached(orbit_key([w]), lambda: a1_coef(w)) for w, T in items)
    v = v1 = 0.0
    for i, (w1, T1) in enumerate(items):
        for j in range(i, len(items)):
            w2, T2 = items[j]
            if len(w1) + len(w2) > LPAIR: continue
            c0, c1 = cached(orbit_key([w1, w2]), lambda: cov_coefs(w1, w2, q0[w1], q0[w2]))
            f = (1 if i == j else 2) * T1 * T2
            v += f * c0; v1 += f * c1
    ans = {"z_inf": trce / 2, "m": float(m), "m1": float(m1), "v": float(v), "v1": float(v1)}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh: json.dump(ans, fh, indent=2)
    print(json.dumps(ans, indent=2), f"\n[{time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
