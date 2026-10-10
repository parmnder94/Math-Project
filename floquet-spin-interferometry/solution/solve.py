"""Reference solution for floquet-spin-interferometry (exact arithmetic, SymPy).

1. Fast Floquet problem.  With tau = Delta t and the carrier phase set to zero, H/Delta = 2Q + eps V1 + eps^2 V2
   + eps^4 V4 + eps^6 V6 with V1 = 2 Jx cos tau and V_o = c_o (P - Q) - 2 d_o X cos 2tau.  In Floquet space the
   states |q, -2> (q in Q) and |p, 0> (p in P) are degenerate at quasienergy 0.  The Feshbach (Brillouin-Wigner)
   matrix at energy 0 is generated exactly by the recurrence W_n = R sum_o V_o W_{n-o}, B_n = Pi sum_o V_o W_{n-o},
   R = (1-Pi)(-F0)^-1(1-Pi).  B_1..B_7 vanish identically (the counterterms a_N, b_N are tuned for this), and the
   first non-zero block B_8 = beta_q (Q - P) + beta_x X is computed here.  Because E = O(eps^8), energy dependence
   enters only at order eps^10.
2. Slow frames.  W_r = R_r D_r gives the connection -(pi/T)D^+ J_r D - (pi/(2T)) K_r, whose secular part on the
   quasi-degenerate cluster is -(3 pi/(2T)) K_r.  With Delta T eps^8 = 1080 pi/(8447 sqrt3) the forward generator
   integrates to pi[-(gamma/2) I + (2/sqrt3) X_phi - (3/2) K_r]; X_phi = e^{-2i phi} QXP + h.c.  The reversed pulse
   has phi -> -phi and the connection sign flipped; U_r = S_r exp(-i pi M_r), Utilde_r = exp(-i pi Mtilde_r) S_r^+,
   S_r = R_r(1) D_r(1).  The common scalar cancels between the arms.
3. Exact algebra.  exp(-i pi M) = (4M^2 - I)/15 + (8i/15)(M^3 - 4M) on the spectrum {+-2, +-1/2} (checked);
   R_r(1) = -(7i/3) J_r + (4i/3) J_r^3 and D_r(1) = Q - i K_r.  p(-) = (1 - Re w)/2, w = Tr(B^+ A P)/2.
"""
import json, os
from fractions import Fraction
import sympy as sp

OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")
I, s3, R = sp.I, sp.sqrt(3), sp.Rational
Jx = sp.Matrix([[0, s3/2, 0, 0], [s3/2, 0, 1, 0], [0, 1, 0, s3/2], [0, 0, s3/2, 0]])     # basis m = 3/2, 1/2, -1/2, -3/2
Jy = sp.Matrix([[0, -I*s3/2, 0, 0], [I*s3/2, 0, -I, 0], [0, I, 0, -I*s3/2], [0, 0, I*s3/2, 0]])
P = sp.diag(0, 1, 1, 0); Q = sp.eye(4) - P
X = (Jx**2 - (7*P + 3*Q)/4).applyfunc(sp.nsimplify)
C_DIAG = {2: 1, 4: R(4, 3), 6: R(148, 45)}
D_OFF = {2: 1, 4: R(5, 3), 6: R(167, 36)}
qflag, ret_photon = [1, 0, 0, 1], [-2, 0, 0, -2]
simp = lambda M: M.applyfunc(lambda e: sp.nsimplify(sp.expand(e)))


def floquet_B(nmax=8, orders=(1, 2, 4, 6), c=C_DIAG, d=D_OFF):
    def vertex(o):
        return {1: Jx, -1: Jx} if o == 1 else {0: c[o]*(P - Q), 2: -d[o]*X, -2: -d[o]*X}

    def apply(o, W):
        out = {}
        for n, M in W.items():
            for l, V in vertex(o).items():
                out[n + l] = out.get(n + l, sp.zeros(4)) + V*M
        return out

    W = {0: {0: P, -2: Q}}
    B = {}
    for n in range(1, nmax + 1):
        acc = {}
        for o in orders:
            if n - o >= 0:
                for k, M in apply(o, W[n - o]).items():
                    acc[k] = acc.get(k, sp.zeros(4)) + M
        acc = {k: simp(M) for k, M in acc.items()}
        Bn = sp.zeros(4)
        for r in range(4):
            if ret_photon[r] in acc:
                Bn[r, :] = acc[ret_photon[r]][r, :]
        B[n] = simp(Bn)
        Wn = {}
        for k, M in acc.items():
            M2 = M.copy()
            for r in range(4):
                M2[r, :] = sp.zeros(1, 4) if k == ret_photon[r] else M2[r, :]*R(-1, 2*qflag[r] + k)
            Wn[k] = M2
        W[n] = Wn
    return B


def probability(beta_q, beta_x, connection=R(3, 2), flip_phase=True, use_D=True, carrier=True):
    kappa = 1080/(8447*s3)                              # Delta T eps^8 / pi
    gamma = 23539/(8447*s3)
    # (Q - P) part of the eighth-order term plus the explicit gap correction -gamma pi Q must leave a pure scalar
    scalar = sp.nsimplify(kappa*beta_q - gamma/2)       # zero: then pi[kappa beta_q (Q-P) - gamma Q] = -(gamma/2) pi I
    cx = sp.nsimplify(kappa*beta_x)                     # coefficient of X_phi (2/sqrt3 expected)
    exp_ia = [1, I, (4 + 3*I)/5]                        # e^{i alpha_r}
    e_iphi = [(2 + I)/sp.sqrt(5), (3 + 2*I)/sp.sqrt(13), (4 + 3*I)/5]
    def Xphi(z):                                        # z = e^{-2 i phi}
        return z*Q*X*P + sp.conjugate(z)*P*X*Q
    def expm_mpi(M):
        M2 = simp(M*M)
        if simp((M2 - 4*sp.eye(4))*(M2 - sp.eye(4)/4)) != sp.zeros(4):   # only for altered (ablation) generators
            import mpmath as mp
            mp.mp.dps = 40
            E = mp.expm(-1j*mp.pi*mp.matrix(M.evalf(45).tolist()))
            return sp.Matrix(4, 4, lambda i, j: sp.Float(mp.re(E[i, j]), 40) + I*sp.Float(mp.im(E[i, j]), 40))
        return simp((4*M2 - sp.eye(4))/15 + R(8, 15)*I*(M2*M - 4*M))
    A, B = sp.eye(4), sp.eye(4)
    Us, Uts = [], []
    for r in range(3):
        ea = exp_ia[r]
        Jr = simp(Jx*sp.re(ea) + Jy*sp.im(ea))
        K = simp(P*Jr*P)
        S = simp((-R(7, 3)*I*Jr + R(4, 3)*I*Jr**3)*((Q - I*K) if use_D else sp.eye(4)))
        z = sp.nsimplify(sp.expand(sp.conjugate(e_iphi[r])**2)) if carrier else 1
        zt = sp.conjugate(z) if flip_phase else z
        M = cx*Xphi(z) - connection*K
        Mt = cx*Xphi(zt) + connection*K
        Us.append(simp(S*expm_mpi(M)))
        Uts.append(simp(expm_mpi(Mt)*S.H))
    A = simp(Us[2]*Us[1]*Us[0]); B = simp(Uts[0]*Uts[1]*Uts[2])
    w = sp.nsimplify(sp.expand((B.H*A*P).trace()/2))
    return sp.nsimplify((1 - sp.re(w))/2), w, scalar, cx


def main():
    B = floquet_B()
    for n in range(1, 8):
        assert B[n] == sp.zeros(4), n
    beta_q = sp.nsimplify(B[8][0, 0])                   # B8 = beta_q (Q - P) + beta_x X
    beta_x = sp.nsimplify(B[8][0, 2]/X[0, 2])
    assert simp(B[8] - (beta_q*(Q - P) + beta_x*X)) == sp.zeros(4)
    p, w, scalar, cx = probability(beta_q, beta_x)
    assert scalar == 0 and sp.simplify(cx - 2/s3) == 0
    frac = Fraction(int(sp.numer(p)), int(sp.denom(p)))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump({"p": f"{frac.numerator}/{frac.denominator}"}, fh)
    print(f"B8 = {beta_q} (Q-P) + {beta_x} X ;  X_phi coefficient {cx} ;  w = {w}")
    print(f"p(-) = {frac} = {float(frac):.17f}")


if __name__ == "__main__":
    main()
