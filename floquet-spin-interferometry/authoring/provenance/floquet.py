"""Exact Floquet-space perturbation theory at quasienergy 0 (Brillouin-Wigner/Feshbach recurrence)."""
import sympy as sp
s3 = sp.sqrt(3); R = sp.Rational
Jx = sp.Matrix([[0, s3/2, 0, 0], [s3/2, 0, 1, 0], [0, 1, 0, s3/2], [0, 0, s3/2, 0]])
P = sp.diag(0, 1, 1, 0); Q = sp.eye(4) - P
X = (Jx**2 - (7*P + 3*Q)/4).applyfunc(sp.nsimplify)
q = [1, 0, 0, 1]
ret_n = [-2, 0, 0, -2]                      # photon index of retained state for each spin state
c = {2: 1, 4: R(4, 3), 6: R(148, 45)}; d = {2: 1, 4: R(5, 3), 6: R(167, 36)}
def vertices(o):
    if o == 1: return {1: Jx, -1: Jx}
    return {0: c[o]*(P - Q), 2: -d[o]*X, -2: -d[o]*X}
def apply_V(o, W):
    out = {}
    for n, M in W.items():
        for l, Vm in vertices(o).items():
            out[n + l] = out.get(n + l, sp.zeros(4)) + Vm*M
    return out
def resolvent(W):
    out = {}
    for n, M in W.items():
        M2 = M.copy()
        for r in range(4):
            e = 2*q[r] + n
            if n == ret_n[r]: M2[r, :] = sp.zeros(1, 4)      # retained state excluded
            else: M2[r, :] = M2[r, :] * R(-1, e)
        out[n] = M2
    return out
def project(W):
    B = sp.zeros(4)
    for r in range(4):
        if ret_n[r] in W: B[r, :] = W[ret_n[r]][r, :]
    return B
orders = [1, 2, 4, 6]
W = {0: {n: sp.zeros(4) for n in [0, -2]}}
W[0][0] = P; W[0][-2] = Q
Bs = {}
for n in range(1, 9):
    acc = {}
    for o in orders:
        if n - o >= 0:
            for k, M in apply_V(o, W[n - o]).items(): acc[k] = acc.get(k, sp.zeros(4)) + M
    acc = {k: M.applyfunc(sp.nsimplify) for k, M in acc.items()}
    Bs[n] = project(acc).applyfunc(sp.nsimplify)
    W[n] = resolvent(acc)
    print(n, Bs[n].tolist())
target = R(23539, 2160)*(Q - P) + R(8447, 540)*X
print('B8 matches golden:', (Bs[8] - target).applyfunc(sp.nsimplify) == sp.zeros(4))
