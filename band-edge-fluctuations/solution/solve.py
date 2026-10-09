"""Reference solution: band-edge DOS fluctuations on random spin-orbit triangle networks.

Model. H_N = P(U_a, U_b) with P = B sz + Ta l(a) + Tb l(b) + Tc l(a)l(b) + h.c. in M_2(C[F_2]),
U_a, U_b uniform random permutation matrices.  X_N = Tr[h(H_N) Q_N] with the Lorentzian
h(x) = eta/((x-x0)^2+eta^2) = (i/2)[R_z(x) - R_zbar(x)],  R_z(x) = 1/(z-x),  z = x0 + i eta.

1. Exact trace expansion (every N):
       X_N = (N-1) tr c_e + sum_{w != e} tr(c_w) (fix(w(sigma,tau)) - 1),   c_w = <delta_w, h(P) delta_e>.
2. Random word maps (Nica; Linial-Puder): if w = g u^d g^-1 with u primitive, fix(w) -> sum_{j|d} j C_j(kappa),
   C_j(kappa) independent Poisson(1/j), one family per primitive class kappa (up to conjugation AND inversion).
   With T_kappa(n) = sum of tr c_w over all conjugates of u^n and u^-n:
       z_inf = tr c_e / 2
       m     = sum_{d>=2} sum_{y cyclically reduced} F(y^d) - tr c_e
       v     = sum_{n,n'} sigma(gcd(n,n')) sum_kappa T(n)T(n')
       k     = sum_{n,n',n''} J2(gcd(n,n',n'')) sum_kappa T(n)T(n')T(n'')      (J2(g) = sum_{j|g} j^2)
3. Resolvent coefficients from a 4x4 self-adjoint linearization of P (Y = 1 l(a^-1) + Tc l(b), Q = -1,
   P = L0 + Y*Y) on the 4-regular Cayley tree:  G(w,e) = Phi_{s_1}...Phi_{s_k} G(e,e), Phi_s = H_s A_s.
   Conjugator sums F(y) = sum_g tr c_{g y g^-1} = tr(Om_(first,last) Phi_y) are resummed exactly.
4. At this z the class sums decay only like q^l with q ~ 0.94, so every term of total power <= 4 is resummed
   EXACTLY: sums over all cyclically reduced words of products of rotated traces are transfer operators on
   tensor powers of Phi, with Om inserted at each rotation's cut (inverse rotations: transposed inverse
   letters); Moebius inversion then restricts to primitive classes.  Terms of total power >= 5 converge like
   (0.17)^l and are summed over primitive classes up to length LMAX.
"""
import itertools, json, os, time
from math import gcd
import numpy as np

LAM, BFIELD, X0, ETA = 0.6, 1.5, 5.6, 0.06
LMAX, EXACT = 12, 4
LNMAX = 6 * LMAX
OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")
INV = [1, 0, 3, 2]                                  # letters a, a^-1, b, b^-1
PAIRS = [(p, n) for p in range(4) for n in range(4) if n != INV[p]]
sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0, -1.0]).astype(complex)
I2 = np.eye(2, dtype=complex)


def hop(s): return np.cos(LAM) * I2 + 1j * np.sin(LAM) * s


def linearization():
    Ta, Tb, Tc = hop(sx), hop(sy), hop(sz)
    K = 4; m, x = slice(0, 2), slice(2, 4)
    b0 = np.zeros((K, K), complex); A = [np.zeros((K, K), complex) for _ in range(4)]
    b0[m, m] = BFIELD * sz - 2 * I2; b0[x, x] = -I2
    A[0][m, m] += Ta; A[1][m, m] += Ta.conj().T; A[2][m, m] += Tb; A[3][m, m] += Tb.conj().T
    A[0][m, x] += I2; A[1][x, m] += I2; A[2][x, m] += Tc; A[3][m, x] += Tc.conj().T
    return K, b0, A


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


def kron_all(ms):
    out = ms[0]
    for m in ms[1:]: out = np.kron(out, m)
    return out


_CC = {}
def cyc_closure(p, K):
    """C on (C^K)^(x)p with tr[(A_1 x ... x A_p) C] = tr(A_1 ... A_p)."""
    if (p, K) not in _CC:
        D = K ** p; C = np.zeros((D, D))
        for idx in itertools.product(range(K), repeat=p):
            C[np.ravel_multi_index(idx[1:] + idx[:1], (K,) * p), np.ravel_multi_index(idx, (K,) * p)] = 1
        _CC[(p, K)] = C
    return _CC[(p, K)]


def transfer(slots):
    """X[a][b] = sum over reduced words a...b of the tensor product (over slots) of the per-letter matrices."""
    Dl = [kron_all([sl[a] for sl in slots]) for a in range(4)]; D = Dl[0].shape[0]
    T = np.zeros((4 * D, 4 * D), complex)
    for a in range(4):
        for b in range(4):
            if b != INV[a]: T[a * D:(a + 1) * D, b * D:(b + 1) * D] = Dl[b]
    W = np.linalg.inv(np.eye(4 * D) - T)
    return [[Dl[a] @ W[a * D:(a + 1) * D, b * D:(b + 1) * D] for b in range(4)] for a in range(4)]


def M_d(pt, d):
    """sum over cyclically reduced y of F(y^d)."""
    K = pt.K; X = transfer([list(pt.Phi)] * d); C = cyc_closure(d, K); I = np.eye(K ** (d - 1))
    return sum(np.trace(np.kron(pt.Om[(f, l)], I) @ X[f][l] @ C) for f in range(4) for l in range(4) if l != INV[f])


def slot_mats(P, inv): return [P.Phi[INV[a]].T for a in range(4)] if inv else list(P.Phi)
def cut_om(P, inv, prev, nxt): return P.Om[(INV[prev], INV[nxt])].T if inv else P.Om[(nxt, prev)]


def ordered_partitions(items):
    if not items: yield []; return
    for k in range(1, len(items) + 1):
        for first in itertools.combinations(items, k):
            rem = tuple(x for x in items if x not in first)
            for tail in ordered_partitions(rem): yield [first] + tail


def psi_marked(slots):
    """sum_{r cyclically reduced} F_0(r^p0) prod_{i>=1} sum_{k} F_i(rot_k(r or r^-1)^{p_i});  slots = [(Point, inv, p)]."""
    K = slots[0][0].K; n = len(slots)
    X = transfer(sum([[slot_mats(P, inv)] * p for (P, inv, p) in slots], []))
    dims = [K ** p for (_, _, p) in slots]; Close = kron_all([cyc_closure(p, K) for (_, _, p) in slots])
    op = lambda i, M: kron_all([np.kron(M, np.eye(K ** (p - 1))) if j == i else np.eye(dims[j]) for j, (P, inv, p) in enumerate(slots)])
    cache = {}

    def Ins(block, pn):
        if (block, pn) not in cache:
            M = None
            for i in block:
                O = op(i, cut_om(slots[i][0], slots[i][1], *pn)); M = O if M is None else M @ O
            cache[(block, pn)] = M
        return cache[(block, pn)]
    kernels = {}

    def kernel(blocks):
        key = tuple(blocks)
        if key not in kernels:
            if not blocks: kernels[key] = X
            else:
                L = kernel(blocks[:-1])
                kernels[key] = [[sum(L[a][p] @ Ins(blocks[-1], (p, q)) @ X[q][b] for (p, q) in PAIRS) for b in range(4)] for a in range(4)]
        return kernels[key]
    rest = tuple(range(1, n)); tot = 0
    for mask in range(1 << len(rest)):
        at0 = tuple(rest[i] for i in range(len(rest)) if mask >> i & 1); inner = tuple(x for x in rest if x not in at0)
        for blocks in ordered_partitions(inner):
            Kb = kernel(blocks)
            for f in range(4):
                for l in range(4):
                    if l == INV[f]: continue
                    S = op(0, slots[0][0].Om[(f, l)])
                    for i in at0: S = S @ op(i, cut_om(slots[i][0], slots[i][1], l, f))
                    tot += np.trace(Close @ S @ Kb[f][l])
    return tot


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


def class_table(pt):
    """rows: primitive classes (mod rotation and inversion) up to LMAX; T[n] = sum over conjugates of u^{+-n} of tr c_w."""
    K = pt.K; rows = []
    for l in range(1, LMAX + 1):
        for u in cyc_reduced(l):
            if is_power(u) or min(rotations(u) + rotations(winv(u))) != u: continue
            T = np.zeros(LNMAX + 1, complex); nmax = LNMAX // l
            for word in (u, winv(u)):
                for r in rotations(word):
                    M = np.eye(K, dtype=complex)
                    for x in r: M = M @ pt.Phi[x]
                    Om = pt.Om[(r[0], r[-1])]; Mn = M
                    for k in range(1, nmax + 1):
                        if k > 1: Mn = Mn @ M
                        T[k] += np.trace(Om @ Mn)
            rows.append(T)
    return np.array(rows)


def main():
    t0 = time.time(); z = complex(X0, ETA)
    PT = {'+': Point(z), '-': Point(np.conj(z))}
    A = PT['+']; trce_h = float(-np.imag(A.trce))
    TA = class_table(A); IT = TA.imag
    memo = {}

    def Pi(ns, signs):          # sum_kappa prod_i T^{(s_i)}(n_i);  T^{(-)} = conj T = T at zbar
        # exact symmetries: the slots commute, and flipping every sign conjugates the sum
        pairs = sorted(zip(ns, signs)); flip = sorted((k, '+' if c == '-' else '-') for k, c in pairs)
        if flip < pairs: return np.conj(Pi(tuple(k for k, _ in flip), tuple(c for _, c in flip)))
        if [k for k, _ in pairs] != list(ns) or [c for _, c in pairs] != list(signs):
            return Pi(tuple(k for k, _ in pairs), tuple(c for _, c in pairs))
        key = (ns, signs)
        if key not in memo:
            if sum(ns) > EXACT:
                memo[key] = np.sum(np.prod([TA[:, k] if s == '+' else np.conj(TA[:, k]) for k, s in zip(ns, signs)], axis=0))
            else:
                psi = sum(psi_marked([(PT[signs[0]], False, ns[0])] + [(PT[s], iv, k) for k, s, iv in zip(ns[1:], signs[1:], invs)])
                          for invs in itertools.product([False, True], repeat=len(ns) - 1))
                corr = sum(q ** (len(ns) - 1) * Pi(tuple(q * k for k in ns), signs) for q in range(2, LNMAX + 1) if q * max(ns) <= LNMAX)
                memo[key] = psi - corr
        return memo[key]

    def ImIm(n1, n2):          # sum_kappa Im T(n1) Im T(n2)
        if n1 + n2 > EXACT: return float(np.sum(IT[:, n1] * IT[:, n2]))
        return 0.5 * np.real(Pi((n1, n2), ('+', '-')) - Pi((n1, n2), ('+', '+')))

    def ImImIm(ns):            # sum_kappa Im T Im T Im T
        if sum(ns) > EXACT: return float(np.sum(IT[:, ns[0]] * IT[:, ns[1]] * IT[:, ns[2]]))
        return float(np.real(sum((-1) ** s.count('-') * Pi(ns, s) for s in itertools.product('+-', repeat=3)) / (-8j)))

    # T(h) = -Im T(z) for the Lorentzian
    m = sum(-np.imag(M_d(A, d) if d <= EXACT else np.sum(TA[:, [q * d for q in range(1, LNMAX // d + 1)]])) for d in range(2, LNMAX + 1)) - trce_h
    # Enumerated class sums via sigma(gcd) = sum_{d | gcd} d and J2(gcd) = sum_{d | gcd} d^2:
    # sum_{n1,n2} sigma(gcd) IT_n1 IT_n2 = sum_d d S_d^2, with S_d = sum_q IT_{dq}; same for k with d^2 S_d^3.
    S = [None] + [IT[:, d::d].sum(axis=1) for d in range(1, LNMAX + 1)]
    small2 = [(a, b) for a in range(1, EXACT) for b in range(1, EXACT - a + 1)]
    small3 = [(a, b, c) for a in range(1, EXACT) for b in range(1, EXACT - a + 1) for c in range(1, EXACT - a - b + 1)]
    v = sum(d * float(np.sum(S[d] ** 2)) for d in range(1, LNMAX + 1))
    v += sum(sum(divisors(gcd(a, b))) * (ImIm(a, b) - float(np.sum(IT[:, a] * IT[:, b]))) for a, b in small2)
    k = -sum(d * d * float(np.sum(S[d] ** 3)) for d in range(1, LNMAX + 1))
    k += sum(sum(d * d for d in divisors(gcd(gcd(a, b), c))) * -(ImImIm((a, b, c)) - float(np.sum(IT[:, a] * IT[:, b] * IT[:, c])))
             for a, b, c in small3)
    ans = {"z_inf": trce_h / 2, "mean_correction": float(m), "variance": float(v), "third_cumulant": float(k)}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh: json.dump(ans, fh, indent=2)
    print(json.dumps(ans, indent=2), f"\n[{time.time() - t0:.0f} s]")


if __name__ == "__main__":
    main()
