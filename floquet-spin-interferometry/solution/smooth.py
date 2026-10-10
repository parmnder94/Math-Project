"""Smooth interpolant F(eps) of p_N, p_N = F(eps_N) + O(eps_N^8), in python-flint arithmetic (Delta = 1).

Moving frame phi = (R_r D_r)^+ psi: H_phi(t) = H0(t) + C(s)/T with H0 2pi-periodic and C(s) = -pi D^+ J_r D - (pi/2) K_r
(forward) or +pi D(1-s)^+ J_r D(1-s) + (pi/2) K_r with carrier phase -phi_r (reversed).  Floquet: U0(t) = P0(t) e^{-i Lam t/2pi},
P0 periodic, Lam = i log U0(2pi).  Averaging the O(1/T) term over each period (error O(1/N) = O(eps^8)) leaves
    i d eta/ds = (N Lam + <C(s)>) eta,   <C> = period average of P0^+ C P0  (trapezoid rule, spectrally exact),
and U_phi(T) = eta(1) since P0(T) = I.  N enters only through N Lam = Lam/(20 pi eps^8), so F is analytic in eps and its
Taylor coefficients up to eps^7 are those of p_N.  One period: Taylor integrator; slow equation: Taylor-series integrator.
"""
from flint import acb_mat, acb, arb, ctx



def setup(dps):
    ctx.dps = dps
    global j, pi, Jx, Jy, P, Q, X, I4, Z4, ALPHAS, PHIS
    j = acb(0, 1); pi = arb.pi(); s3 = arb(3).sqrt()
    Jx = acb_mat([[0, s3/2, 0, 0], [s3/2, 0, 1, 0], [0, 1, 0, s3/2], [0, 0, s3/2, 0]])
    Jy = acb_mat([[0, -j*s3/2, 0, 0], [j*s3/2, 0, -j, 0], [0, j, 0, -j*s3/2], [0, 0, j*s3/2, 0]])
    P = acb_mat([[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]])
    I4 = acb_mat([[1 if a == b else 0 for b in range(4)] for a in range(4)]); Z4 = acb_mat(4, 4)
    Q = I4 - P
    X = Jx*Jx - (7*P + 3*Q)/4
    ALPHAS = [arb(0), pi/3, 2*pi/3]
    PHIS = [pi/6, pi/4, pi/3]


dag = lambda A: A.conjugate().transpose()


def one_period(eps, ph, n=64, order=30):
    a = eps**2 + arb(4)/3*eps**4 + arb(148)/45*eps**6
    b = eps**2 + arb(5)/3*eps**4 + arb(167)/36*eps**6
    H00 = 2*Q + a*(P - Q)
    h = 2*pi/n
    U = I4; out = [U]
    fact = [arb.fac_ui(k) for k in range(order + 1)]
    half = pi/2
    for st in range(n):
        t0 = st*h
        c1 = [2*eps*(t0 + ph + k*half).cos()/fact[k] for k in range(order + 1)]
        c2 = [-2*b*arb(2)**k*(2*(t0 + ph) + k*half).cos()/fact[k] for k in range(order + 1)]
        Uk = [U]
        for k in range(order):
            S1 = Z4; S2 = Z4
            for i in range(k + 1):
                S1 = S1 + c1[i]*Uk[k - i]; S2 = S2 + c2[i]*Uk[k - i]
            Uk.append(-j*(H00*Uk[k] + Jx*S1 + X*S2)/(k + 1))
        U = Uk[order]
        for k in range(order - 1, -1, -1):
            U = U*h + Uk[k]
        out.append(U)
    return out


def logm_near_identity(U):
    E, R = U.mid().eig(right=True, algorithm="approx")
    Rinv = R.inv()
    D = acb_mat([[j*E[a].log() if a == b else 0 for b in range(4)] for a in range(4)])
    return R*D*Rinv


def expm(A):
    return A.exp()


def pulse(eps, alpha, ph, reverse, n=64, order=30, steps=16, torder=36, dressed=True):
    Us = one_period(eps, -ph if reverse else ph, n, order)
    Lam = logm_near_identity(Us[-1].mid())
    P0 = [Us[k].mid()*expm(j*Lam*k/n) for k in range(n)]
    def avg(M):
        if not dressed: return P*M*P + Q*M*Q          # ablation: secular part without micromotion dressing
        S = Z4
        for A in P0: S = S + dag(A)*M*A
        return S/n
    Jr = Jx*alpha.cos() + Jy*alpha.sin(); K = P*Jr*P
    def C(s):
        if reverse:
            Dm = expm(-j*pi*(1 - s)*K/2); return pi*dag(Dm)*Jr*Dm + pi/2*K
        Dm = expm(-j*pi*s*K/2); return -pi*dag(Dm)*Jr*Dm - pi/2*K
    basis = lambda s: [arb(1), (pi*s/2).cos(), (pi*s/2).sin(), (pi*s).cos(), (pi*s).sin()]
    nodes = [arb(k)/5 for k in range(5)]
    Binv = acb_mat([basis(s) for s in nodes]).inv()
    Cs = [avg(C(s)) for s in nodes]
    W = []
    for i in range(5):
        S = Z4
        for k in range(5): S = S + Binv[i, k]*Cs[k]
        W.append(S.mid())
    st = arb(37)/41
    chk = avg(C(st))
    for bv, Wi in zip(basis(st), W): chk = chk - bv*Wi
    assert max(abs(complex(chk[a, b].mid())) for a in range(4) for b in range(4)) < 1e-25
    NL = Lam/(20*pi*eps**8)
    hs = arb(1)/steps; eta = I4
    for st_ in range(steps):
        s0 = st_*hs
        tc = []
        for k in range(torder + 1):
            fk = arb.fac_ui(k)
            w = [arb(0 if k else 1),
                 (pi/2)**k*(pi*s0/2 + k*pi/2).cos()/fk, (pi/2)**k*(pi*s0/2 + k*pi/2).sin()/fk,
                 pi**k*(pi*s0 + k*pi/2).cos()/fk, pi**k*(pi*s0 + k*pi/2).sin()/fk]
            A = Z4
            for wi, Wi in zip(w, W): A = A + wi*Wi
            if k == 0: A = A + NL
            tc.append((-j*A).mid())
        ek = [eta]
        for k in range(torder):
            S = Z4
            for i in range(k + 1): S = S + tc[i]*ek[k - i]
            ek.append((S/(k + 1)).mid())
        eta = ek[torder]
        for k in range(torder - 1, -1, -1):
            eta = eta*hs + ek[k]
        eta = eta.mid()
    S = expm(-j*pi*Jr)*expm(-j*pi*K/2)
    return S*eta if not reverse else eta*dag(S)


def F(eps, **kw):
    eps = arb(eps)
    A = I4; B = I4
    for r in range(3): A = pulse(eps, ALPHAS[r], PHIS[r], False, **kw)*A
    for r in (2, 1, 0): B = pulse(eps, ALPHAS[r], PHIS[r], True, **kw)*B
    M = dag(B)*A*P
    w = (M[0, 0] + M[1, 1] + M[2, 2] + M[3, 3])/2
    return ((1 - w.real)/2).mid()


