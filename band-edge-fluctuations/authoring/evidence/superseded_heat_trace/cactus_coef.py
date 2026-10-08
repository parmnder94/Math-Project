"""Independent coefficients c_w = <delta_w, exp(-beta P) delta_e> from the triangle-cactus cavity solution.
The Cayley graph of F_2 w.r.t. {a, b, ab} is a tree of triangles {x, bx, abx}; G(w,e) factorizes over the
triangles crossed by the cactus geodesic: G(w,e) = Theta_m ... Theta_1 G(e,e), Theta = (M_r^-1 C_r)_exit.
No linearization is used."""
import numpy as np
from cactus import modelA, I2, Z2
a_, A_, b_, B_ = 0, 1, 2, 3
SINGLE = {b_: (0, 1), B_: (1, 0), a_: (1, 2), A_: (2, 1)}


def steps(w):
    path = list(reversed(w)); out = []; i = 0
    while i < len(path):
        s = path[i]; nxt = path[i + 1] if i + 1 < len(path) else None
        if s == b_ and nxt == a_: out.append((0, 2)); i += 2; continue      # x -> bx -> abx is one c-bond
        if s == A_ and nxt == B_: out.append((2, 0)); i += 2; continue
        out.append(SINGLE[s]); i += 1
    return out


class CactusCoef:
    def __init__(self, B, lam, beta, center=0.5, ra=8.0, rb=3.5, nodes=160):
        c = modelA(B, lam)
        ts = np.pi / 2 + 2 * np.pi * np.arange(nodes) / nodes
        zs = center + ra * np.cos(ts) + 1j * rb * np.sin(ts)
        dz = (-ra * np.sin(ts) + 1j * rb * np.cos(ts)) * (2 * np.pi / nodes) / (2j * np.pi)
        self.w = np.exp(-beta * zs) * dz
        Sig, _ = c.solve(zs[0], iters=200000, tol=1e-15)
        Th = {(r, e): [] for r in range(3) for e in c.others[r]}; G = []
        for z in zs:
            Sig, res = c.newton(z, Sig); assert res < 1e-12, res
            St = sum(Sig); base = z * I2 - c.D - St
            for r in range(3):
                r1, r2 = c.others[r]
                M = np.block([[base + Sig[r1], Z2], [Z2, base + Sig[r2]]]) - c.hp[r]
                X = np.linalg.solve(M, c.C[r]); Th[(r, r1)].append(X[:2]); Th[(r, r2)].append(X[2:])
            G.append(c.G(Sig, z))
        self.Theta = {k: np.array(v) for k, v in Th.items()}; self.G = np.array(G)

    def trace_coef(self, w):
        M = self.G
        for st in steps(w): M = self.Theta[st] @ M
        return np.einsum('n,nii->', self.w, M).real if False else np.einsum('n,nii->', self.w, M)
