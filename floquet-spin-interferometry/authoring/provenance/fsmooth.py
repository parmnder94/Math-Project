"""Smooth high-precision interpolant F(eps) of p_N (Delta = 1), with p_N = F(eps_N) + O(eps_N^8).

Moving frame phi = (R_r D_r)^+ psi:  H_phi(t) = H0(t) + C(s)/T with H0 2pi-periodic and
C(s) = -pi D^+ J_r D - (pi/2) K_r (forward), +pi D(1-s)^+ J_r D(1-s) + (pi/2) K_r and phase -phi_r (reversed).
Floquet: U0(t) = P0(t) exp(-i Lam t/2pi), P0 periodic, Lam = i log U0(2pi) (all four eigenphases near 0).
Averaging the O(1/T) term over each period (error O(1/N) = O(eps^8)) gives, with s = t/T and N Lam = Lam/(20 pi eps^8),
    i d eta/ds = (N Lam + <C(s)>) eta,   <C> = (1/2pi) int_0^{2pi} P0^+ C P0 dt  (trapezoid rule: P0 is periodic),
and U_phi(T) = eta(1) because P0(T) = I.  Everything below eps^8 is exact; N is treated as the continuous
variable N = 1/(20 pi eps^8), so F is analytic in eps and its Taylor coefficients are those of p_N.
"""
import mpmath as mp

mp.mp.dps = 60
j = mp.mpc(0, 1)
s3 = mp.sqrt(3)
Jx = mp.matrix([[0, s3/2, 0, 0], [s3/2, 0, 1, 0], [0, 1, 0, s3/2], [0, 0, s3/2, 0]])
Jy = mp.matrix([[0, -j*s3/2, 0, 0], [j*s3/2, 0, -j, 0], [0, j, 0, -j*s3/2], [0, 0, j*s3/2, 0]])
P = mp.diag([0, 1, 1, 0]); Q = mp.eye(4) - P
X = Jx*Jx - (7*P + 3*Q)/4
ALPHAS = [mp.mpf(0), mp.pi/3, 2*mp.pi/3]
PHIS = [mp.pi/6, mp.pi/4, mp.pi/3]
dag = lambda A: A.transpose_conj()


def one_period(eps, ph, n=96, order=40):
    """U0 at t_k = 2 pi k/n, k = 0..n, for H0(t) = 2Q + 2 eps cos(t+ph) Jx + a(P-Q) - 2b cos(2(t+ph)) X."""
    a = eps**2 + mp.mpf(4)/3*eps**4 + mp.mpf(148)/45*eps**6
    b = eps**2 + mp.mpf(5)/3*eps**4 + mp.mpf(167)/36*eps**6
    H00 = 2*Q + a*(P - Q)
    h = 2*mp.pi/n
    U = mp.eye(4); out = [U]
    fact = [mp.factorial(k) for k in range(order + 1)]
    for st in range(n):
        t0 = st*h
        c1 = [2*eps*mp.cos(t0 + ph + k*mp.pi/2)/fact[k] for k in range(order + 1)]          # coefficient of Jx
        c2 = [-2*b*2**k*mp.cos(2*(t0 + ph) + k*mp.pi/2)/fact[k] for k in range(order + 1)]   # coefficient of X
        Uk = [U]
        for k in range(order):
            S1 = mp.zeros(4); S2 = mp.zeros(4)
            for i in range(k + 1):
                S1 += c1[i]*Uk[k - i]; S2 += c2[i]*Uk[k - i]
            Uk.append(-j*(H00*Uk[k] + Jx*S1 + X*S2)/(k + 1))
        U = Uk[order]
        for k in range(order - 1, -1, -1):
            U = U*h + Uk[k]
        out.append(U)
    return out


def logm_near_identity(U):
    E, V = mp.eig(U)
    return V*mp.diag([j*mp.log(e) for e in E])*mp.inverse(V)          # Lam = i log U


def pulse(eps, alpha, ph, reverse, n=96, order=40, steps=24, torder=48):
    Us = one_period(eps, -ph if reverse else ph, n, order)
    Lam = logm_near_identity(Us[-1])
    P0 = [Us[k]*mp.expm(j*Lam*k/n) for k in range(n)]
    avg = lambda M: sum((dag(A)*M*A for A in P0), mp.zeros(4))/n
    Jr = Jx*mp.cos(alpha) + Jy*mp.sin(alpha); K = P*Jr*P
    def C(s):
        if reverse:
            Dm = mp.expm(-j*mp.pi*(1 - s)*K/2); return mp.pi*dag(Dm)*Jr*Dm + mp.pi/2*K
        Dm = mp.expm(-j*mp.pi*s*K/2); return -mp.pi*dag(Dm)*Jr*Dm - mp.pi/2*K
    # C(s) = W0 + W1 cos(pi s/2) + W2 sin(pi s/2) + W3 cos(pi s) + W4 sin(pi s); fit on 5 nodes, check on a 6th
    nodes = [mp.mpf(k)/5 for k in range(5)]
    basis = lambda s: [1, mp.cos(mp.pi*s/2), mp.sin(mp.pi*s/2), mp.cos(mp.pi*s), mp.sin(mp.pi*s)]
    Binv = mp.inverse(mp.matrix([basis(s) for s in nodes]))
    Cs = [avg(C(s)) for s in nodes]
    W = [sum((Binv[i, k]*Cs[k] for k in range(5)), mp.zeros(4)) for i in range(5)]
    st = mp.mpf(37)/41
    chk = avg(C(st)) - sum((bv*Wi for bv, Wi in zip(basis(st), W)), mp.zeros(4))
    assert mp.mnorm(chk, 1) < mp.mpf(10)**(-mp.mp.dps + 10), chk
    NL = Lam/(20*mp.pi*eps**8)
    # Taylor-series integration of i eta' = (NL + <C>(s)) eta
    hs = mp.mpf(1)/steps; eta = mp.eye(4)
    for st_ in range(steps):
        s0 = st_*hs
        tc = []
        for k in range(torder + 1):
            fk = mp.factorial(k)
            w = [0 if k else 1,
                 (mp.pi/2)**k*mp.cos(mp.pi*s0/2 + k*mp.pi/2)/fk, (mp.pi/2)**k*mp.sin(mp.pi*s0/2 + k*mp.pi/2)/fk,
                 mp.pi**k*mp.cos(mp.pi*s0 + k*mp.pi/2)/fk, mp.pi**k*mp.sin(mp.pi*s0 + k*mp.pi/2)/fk]
            A = sum((wi*Wi for wi, Wi in zip(w, W)), mp.zeros(4))
            if k == 0: A = A + NL
            tc.append(-j*A)
        ek = [eta]
        for k in range(torder):
            ek.append(sum((tc[i]*ek[k - i] for i in range(k + 1)), mp.zeros(4))/(k + 1))
        eta = ek[torder]
        for k in range(torder - 1, -1, -1):
            eta = eta*hs + ek[k]
    S = mp.expm(-j*mp.pi*Jr)*mp.expm(-j*mp.pi*K/2)
    return S*eta if not reverse else eta*dag(S)


def F(eps, **kw):
    eps = mp.mpf(eps)
    A = mp.eye(4); B = mp.eye(4)
    for r in range(3): A = pulse(eps, ALPHAS[r], PHIS[r], False, **kw)*A
    for r in (2, 1, 0): B = pulse(eps, ALPHAS[r], PHIS[r], True, **kw)*B
    w = sum((dag(B)*A*P)[i, i] for i in range(4))/2
    return (1 - mp.re(w))/2


if __name__ == "__main__":
    import sys, time
    for e in sys.argv[1:]:
        t = time.time(); print(e, mp.nstr(F(e), 40), f"{time.time()-t:.0f}s", flush=True)
