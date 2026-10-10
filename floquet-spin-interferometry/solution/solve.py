"""Reference solution for floquet-spin-interferometry.

1. Fast Floquet problem.  With tau = Delta t and the carrier phase set to zero, H/Delta = 2Q + eps V1 + eps^2 V2
   + eps^4 V4 + eps^6 V6 with V1 = 2 Jx cos tau and V_o = c_o (P - Q) - 2 d_o X cos 2tau.  In Floquet space the
   states |q, -2> (q in Q) and |p, 0> (p in P) are degenerate at quasienergy 0.  The Feshbach (Brillouin-Wigner)
   matrix at energy 0 is generated exactly by the recurrence W_n = R sum_o V_o W_{n-o}, B_n = Pi sum_o V_o W_{n-o},
   R = (1-Pi)(-F0)^-1(1-Pi).  B_1..B_7 vanish identically (a_N and b_N cancel orders 2, 4, 6; odd orders vanish by
   photon parity) and the first non-zero block is B_8 = beta_q (Q - P) + beta_x X, computed here in exact arithmetic.
   Because the quasienergy is O(eps^8), energy dependence enters only at order eps^10.
2. Carrier phase and slow frames.  A Fourier vertex of photon shift l picks up e^{i l phi}, so X -> X_phi =
   e^{-2i phi} QXP + h.c.  The frame W_r = R_r D_r gives the connection -(pi/T)D^+ J_r D - (pi/(2T)) K_r, whose
   secular part on the cluster is -(3 pi/(2T)) K_r.  With Delta T_N eps_N^8 = LAM = 1/10 the forward pulse is
   U_r = S_r exp(-i G_r), G_r = LAM [beta_q (Q-P) + beta_x X_phi] - (3 pi/2) K_r, S_r = R_r(1) D_r(1).
   The literal reversal H(T_N - t) has phi -> -phi and the connection sign flipped: Utilde_r = exp(-i Gtilde_r) S_r^+,
   Gtilde_r = LAM [beta_q (Q-P) + beta_x X_{-phi}] + (3 pi/2) K_r.  Since 2 Delta T_N = 4 pi N, exp(-2i Delta T_N Q) = I.
3. Evaluation.  R_r(1) = -(7i/3) J_r + (4i/3) J_r^3 and D_r(1) = Q - i K_r exactly; the exponentials and the overlap
   w = Tr(B^+ A P)/2 are evaluated in 40-digit arithmetic, and p(-) = (1 - Re w)/2.
"""
import json, os
import mpmath as mp
import sympy as sp

OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")
LAM = sp.Rational(1, 10)                                 # Delta T_N eps_N^8
ALPHAS = [0, sp.pi/3, 2*sp.pi/3]
PHIS = [sp.pi/6, sp.pi/4, sp.pi/3]
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


def probability(beta_q, beta_x, lam=LAM, alphas=ALPHAS, phis=PHIS, connection=R(3, 2),
                flip_phase=True, use_D=True, carrier=True, dps=40):
    mp.mp.dps = dps
    tomp = lambda M: mp.matrix([[mp.mpc(sp.re(e).evalf(dps + 10), sp.im(e).evalf(dps + 10)) for e in M.row(i)]
                                for i in range(4)])
    Pm, Qm, Xm = tomp(P), tomp(Q), tomp(X)
    lam, bq, bx, conn = mp.mpf(sp.Rational(lam).p)/sp.Rational(lam).q, mp.mpf(beta_q.p)/beta_q.q, \
        mp.mpf(beta_x.p)/beta_x.q, mp.mpf(sp.Rational(connection).p)/sp.Rational(connection).q
    j = mp.mpc(0, 1)
    dag = lambda A: A.transpose_conj()
    Us, Uts = [], []
    for a, ph in zip(alphas, phis):
        Jr = tomp(Jx*sp.cos(a) + Jy*sp.sin(a))
        K = Pm*Jr*Pm
        S = (-j*mp.mpf(7)/3*Jr + j*mp.mpf(4)/3*Jr**3)*((Qm - j*K) if use_D else mp.eye(4))
        z = mp.exp(-2*j*mp.mpf(sp.N(ph, dps + 10))) if carrier else mp.mpf(1)
        zt = mp.conj(z) if flip_phase else z
        Xphi = lambda zz: zz*Qm*Xm*Pm + mp.conj(zz)*Pm*Xm*Qm
        G = lam*(bq*(Qm - Pm) + bx*Xphi(z)) - conn*mp.pi*K
        Gt = lam*(bq*(Qm - Pm) + bx*Xphi(zt)) + conn*mp.pi*K
        Us.append(S*mp.expm(-j*G))
        Uts.append(mp.expm(-j*Gt)*dag(S))
    A = Us[2]*Us[1]*Us[0]
    B = Uts[0]*Uts[1]*Uts[2]
    w = sum((dag(B)*A*Pm)[i, i] for i in range(4))/2
    return (1 - mp.re(w))/2, w


def main():
    B = floquet_B()
    for n in range(1, 8):
        assert B[n] == sp.zeros(4), n
    beta_q = sp.nsimplify(B[8][0, 0])                   # B_8 = beta_q (Q - P) + beta_x X
    beta_x = sp.nsimplify(B[8][0, 2]/X[0, 2])
    assert simp(B[8] - (beta_q*(Q - P) + beta_x*X)) == sp.zeros(4)
    p, w = probability(beta_q, beta_x)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump({"p": float(p)}, fh)
    print(f"B8 = {beta_q} (Q-P) + {beta_x} X ;  w = {mp.nstr(w, 20)}")
    print(f"p(-) = {mp.nstr(p, 30)}")


if __name__ == "__main__":
    main()
