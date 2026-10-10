# Independent degree check: count preimages of the regular values +-e_i (i=1..4) of d/|d| on the l-torus.
# For v = s*e_i: the d_k (k != i, k<=4) vanish at l_k = pi n_k/3, then d5 must vanish as a function of l_i,
# with s*d_i > 0. Local degree = s*(-1)^(i-1)*sign det d(d_k (k != i), d_5)/d(l_1..l_4).
import itertools, math, cmath
import numpy as np
from scipy.optimize import brentq
w = cmath.exp(2j * math.pi / 3)


def d5(l):
    p = [math.cos(3 * x) for x in l]; u = [cmath.exp(2j * x) for x in l]
    chi = lambda a, b: (1 + a + b - a * b) / 2
    S = lambda a, v: (1 + a + (1 - a) * v) / 2
    Xi = lambda a, b: ((1 + b + b * b) + a * (1 + w * w * b + w * b * b) + a * a * (1 + w * b + w * w * b * b)) / 3
    pc = p + [p[0]]
    Z1 = Xi(u[0], u[1]) * Xi(u[2], u[3]) * S(p[0], u[0]) * S(p[1], u[1]) * S(p[2], u[2]) * S(p[3], u[3])
    Z2 = Xi(u[0], u[2]) * Xi(u[1], u[3]) * S(pc[1], u[0]) * S(pc[2], u[1]) * S(pc[3], u[2]) * S(pc[4], u[3])
    return chi(pc[0], pc[1]) * chi(pc[1], pc[2]) * chi(pc[2], pc[3]) * chi(pc[3], pc[4]) * (0.2 + Z1.real + 0.6 * Z2.real)


def dvec(l):
    return [math.sin(3 * x) for x in l] + [d5(l)]


def jac(l, h=1e-6):
    J = np.zeros((5, 4))
    for k in range(4):
        lp = list(l); lm = list(l); lp[k] += h; lm[k] -= h
        J[:, k] = (np.array(dvec(lp)) - np.array(dvec(lm))) / (2 * h)
    return J


G = np.linspace(0, 2 * math.pi, 3001)
for i in range(4):
    deg = {1: 0, -1: 0}; nroots = 0; mindet = 1e9; minsep = 1e9
    for n in itertools.product(range(6), repeat=3):
        def f(t):
            l = [0.0] * 4; it = iter(n)
            for k in range(4):
                l[k] = t if k == i else math.pi * next(it) / 3
            return l
        vals = [d5(f(t)) for t in G]
        roots = []
        for a, b, fa, fb in zip(G[:-1], G[1:], vals[:-1], vals[1:]):
            if fa * fb < 0:
                roots.append(brentq(lambda t: d5(f(t)), a, b, xtol=1e-14))
            elif fa == 0:
                raise SystemExit("grid hit a root exactly; shift grid")
        for t in roots:
            l = f(t); J = jac(l); rows = [k for k in range(4) if k != i] + [4]
            det = np.linalg.det(J[rows, :]); mindet = min(mindet, abs(det)); nroots += 1
            s = 1 if math.sin(3 * t) > 0 else -1
            minsep = min(minsep, abs(math.sin(3 * t)))
            deg[s] += s * (-1) ** i * int(np.sign(det))
    print(f"i={i+1}: roots {nroots}, deg at +e{i+1} = {deg[1]}, deg at -e{i+1} = {deg[-1]}, "
          f"min|det| {mindet:.3g}, min|d_i| at roots {minsep:.3g}", flush=True)
