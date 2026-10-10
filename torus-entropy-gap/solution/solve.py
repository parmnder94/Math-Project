"""Reference solution for torus-entropy-gap (exact rational arithmetic, SymPy only).

Model (instruction): T = 1 + (4/5)cos x, theta = x + y + (1/2)cos x, R = R(theta), reservoirs along r_a = R e_a with
friction gam = (1, 3) and temperatures (T, 2T), magnetic term (1/2) J V, explicit force F.  In the reservoir frame the
fast matrices are constant: B0 = diag(1, 3) - (1/2) J, Q0 = diag(1, 6), entropy weights W0 = diag(1, 3/2).

1. Fast covariance.  B0 C0 + C0 B0^T = 2 Q0 gives C0 = [[55, -6], [-6, 103]]/52; the scaled velocity w = sqrt(m) V is
   locally Gaussian with covariance Sigma = T R C0 R^T.  The given force is F = div Sigma + Sigma grad(2 log T) + B h with
   B = R B0 R^T and h = (1/(3T^2), T'/(2T)) (checked symbolically below, not assumed).
2. Small-mass limit.  The generator of (X, w) is L/m + A/sqrt(m), L = -(Bw).grad_w + Q:grad_w grad_w, Q = T R Q0 R^T,
   A = w.grad_X + F.grad_w.  The limiting generator is L_eff = -<A L^-1 A>, computed here directly on coordinate
   functions: drift b = M F + (d_l M) Sigma_l, diffusion D = sym(M Sigma), M = B^-1.  The stationary density is
   rho = T^2/<T^2> (div j = 0 is verified), with current j = rho b - div(D rho).
3. Observed irreversibility.  The quadratic variation of the x path gives D_xx = T a(theta) and its covariation with
   a(theta_t) gives a'(theta); together they fix theta mod pi, hence y mod pi.  The dynamics is invariant under
   y -> y + pi, so the two lifts carry no relative entropy and the x-path rate equals the full-position rate
   sigma = int j^T D^-1 j / rho.
4. Finite heat.  The Ito form of the entropy flow is Sdot_m = (1/m) E[w^T W w/T - sum gam], W = R W0 R^T.  With
   Phi = w^T Z w, Z = R Z0 R^T/T, B0^T Z0 + Z0 B0 = W0, one has L Phi = -w^T W w/T + 2 tr(Q0 Z0), so stationarity gives
   Sdot_m - kappa/m = m^(-1/2) E_m[A Phi] exactly, kappa = 2 tr(Q0 Z0) - sum gam = 3/104.  Expanding the stationary
   density f = f0 + sqrt(m) f1 + ..., f0 = rho N(0, Sigma), L* f1 = -A* f0, gives
   Sdot_0 = lim (Sdot_m - kappa/m) = int A Phi f1 = -int f0 A L^-1 A Phi.
5. Exact averages.  All integrands are polynomials in cos x, sin x, cos theta, sin theta, 1/T.  At fixed x, theta is
   uniform in y, so the theta-average is done monomial by monomial; odd powers of sin x integrate to zero; the rest
   is a Laurent polynomial in T, whose moments <T^n> are exact (negative n via Legendre polynomials, sqrt(1 - a^2) = 3/5).
The answer is Delta = Sdot_0 - sigma as an exact fraction.
"""
import json
import os
from functools import lru_cache

import sympy as sp
from sympy import QQ
from sympy.polys.rings import ring

OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")
Rg, Cx, Sx, c, s, iT, w1, w2 = ring("Cx Sx c s iT w1 w2", QQ)
GENS = (Cx, Sx, c, s, iT, w1, w2)
WV = (w1, w2)


def toQQ(v):
    v = sp.Rational(v)
    return QQ(int(v.p), int(v.q))


def mm(A, B): return [[sum(A[i][l]*B[l][j] for l in range(2)) for j in range(2)] for i in range(2)]
def scal(f, A): return [[f*A[i][j] for j in range(2)] for i in range(2)]


class Model:
    """T = 1 + a cos x, theta = x + y + lam cos x; fast matrices from (gam, temperature factors tau, magnetic b).
    The force is built as F = div Sigma + Sigma grad(k log T) + B h, h = (alpha/T^k, beta T' T^(1-k)), and is compared
    with the explicit force of the instruction in check_force()."""

    def __init__(self, a, lam, gam, tau, b, k, alpha, beta, rho_pow=None, wtau=None):
        self.a, self.lam, self.k = toQQ(a), toQQ(lam), k
        self.gam = gam
        self.T = 1 + self.a*Cx
        self.Tp = -self.a*Sx
        self.R = [[c, -s], [s, c]]
        self.RT = [[c, s], [-s, c]]
        B0 = sp.Matrix([[gam[0], sp.Rational(b)], [-sp.Rational(b), gam[1]]])
        Q0 = sp.diag(gam[0]*tau[0], gam[1]*tau[1])
        wt = tau if wtau is None else wtau
        W0 = sp.diag(sp.Rational(gam[0], wt[0]), sp.Rational(gam[1], wt[1]))
        cc = sp.symbols('c0:3')
        C0 = sp.Matrix([[cc[0], cc[1]], [cc[1], cc[2]]])
        C0 = C0.subs(sp.solve(list(B0*C0 + C0*B0.T - 2*Q0), cc))
        zz = sp.symbols('z0:3')
        Z0 = sp.Matrix([[zz[0], zz[1]], [zz[1], zz[2]]])
        Z0 = Z0.subs(sp.solve(list(B0.T*Z0 + Z0*B0 - W0), zz))
        self.sym = dict(B0=B0, Q0=Q0, W0=W0, C0=C0, Z0=Z0)
        q = lambda A: [[toQQ(A[i, j]) for j in range(2)] for i in range(2)]
        self.B0, self.Q0, self.C0, self.Z0, self.M0 = q(B0), q(Q0), q(C0), q(Z0), q(B0.inv())
        rot = lambda A0: mm(mm(self.R, A0), self.RT)
        self.B, self.M = rot(self.B0), rot(self.M0)
        self.Sig = scal(self.T, rot(self.C0))
        self.thx = 1 - self.lam*Sx
        gpsi = [k*self.Tp*iT, Rg(0)]
        self.h = [toQQ(alpha)*iT**k, toQQ(beta)*self.Tp*(iT**(k - 1) if k >= 1 else self.T**(1 - k))]
        divS = [self.dx(self.Sig[i][0]) + self.dy(self.Sig[i][1]) for i in range(2)]
        self.Bh = [sum(self.B[i][j]*self.h[j] for j in range(2)) for i in range(2)]
        self.F = [divS[i] + sum(self.Sig[i][j]*gpsi[j] for j in range(2)) + self.Bh[i] for i in range(2)]
        self.rho_pow = k if rho_pow is None else rho_pow
        self.rho = self.T**self.rho_pow

    def dth(self, P): return P.diff(c)*(-s) + P.diff(s)*c
    def dx(self, P): return P.diff(Cx)*(-Sx) + P.diff(Sx)*Cx + self.dth(P)*self.thx + P.diff(iT)*(self.a*Sx*iT**2)
    def dy(self, P): return self.dth(P)


def subs_w(P, A):
    return P.compose([(w1, A[0][0]*w1 + A[0][1]*w2), (w2, A[1][0]*w1 + A[1][1]*w2)])


def wdeg_part(P, d):
    return Rg({mo: cf for mo, cf in P.terms() if mo[5] + mo[6] == d})


def Aop(m, P, F=None):
    F = m.F if F is None else F
    return w1*m.dx(P) + w2*m.dy(P) + F[0]*P.diff(w1) + F[1]*P.diff(w2)


def solve_DB(m, P, d):
    """solve -(B0 w').grad' u = P for P homogeneous of degree d in w' (constant rational matrix on monomials)"""
    mons = [(d - i, i) for i in range(d + 1)]
    Mt = sp.zeros(d + 1, d + 1)
    for col, (e1, e2) in enumerate(mons):
        mon = w1**e1*w2**e2
        img = -((m.B0[0][0]*w1 + m.B0[0][1]*w2)*mon.diff(w1) + (m.B0[1][0]*w1 + m.B0[1][1]*w2)*mon.diff(w2))
        for mo, cf in img.terms():
            Mt[mons.index((mo[5], mo[6])), col] += sp.Rational(int(cf.numerator), int(cf.denominator))
    Minv = Mt.inv()
    parts = {mn: Rg(0) for mn in mons}
    for mo, cf in P.terms():
        rest = list(mo); rest[5] = rest[6] = 0
        parts[(mo[5], mo[6])] += Rg({tuple(rest): cf})
    u = Rg(0)
    for i, (e1, e2) in enumerate(mons):
        cf = sum((parts[mons[j]]*toQQ(Minv[i, j]) for j in range(d + 1)), Rg(0))
        u += cf*w1**e1*w2**e2
    return u


def Linv_odd(m, P):
    """L u = P for P cubic + linear in w, solved in the reservoir frame w = R w' where L has constant coefficients"""
    Pr = subs_w(P, m.R)
    p3, p1 = wdeg_part(Pr, 3), wdeg_part(Pr, 1)
    assert Pr == p3 + p1
    u3 = solve_DB(m, p3, 3)
    lap = m.T*sum((m.Q0[i][j]*u3.diff(WV[i]).diff(WV[j]) for i in range(2) for j in range(2)), Rg(0))
    u1 = solve_DB(m, p1 - lap, 1)
    return subs_w(u3 + u1, m.RT)


def gauss_avg(m, P):
    """<P(w)> with w = R w', w' ~ N(0, T C0) (Isserlis)"""
    Pr = subs_w(P, m.R)

    @lru_cache(None)
    def mom(i, j):
        idx = (0,)*i + (1,)*j
        if len(idx) % 2:
            return QQ(0), 0

        def pairings(l):
            if not l:
                yield []
                return
            for t in range(1, len(l)):
                for rest in pairings(l[1:t] + l[t + 1:]):
                    yield [(l[0], l[t])] + rest
        tot = QQ(0)
        for pr in pairings(list(idx)):
            v = QQ(1)
            for u, w in pr:
                v *= m.C0[u][w]
            tot += v
        return tot, len(idx)//2
    out = Rg(0)
    for mo, cf in Pr.terms():
        val, npair = mom(mo[5], mo[6])
        if val == 0:
            continue
        rest = list(mo); rest[5] = rest[6] = 0
        out += Rg({tuple(rest): cf*val})*m.T**npair
    return out


@lru_cache(None)
def Tmoment(a, n):
    """<(1 + a cos x)^n> over a period, exact for |a| < 1 (a a sympy Rational)"""
    if n >= 0:
        return sum(sp.binomial(n, j)*a**j*(sp.factorial2(j - 1)/sp.factorial2(j) if j % 2 == 0 else 0)
                   for j in range(n + 1))
    r = sp.sqrt(1 - a**2)
    return sp.nsimplify(r**n*sp.legendre(-n - 1, 1/r))


def torus_average(m, P):
    """exact average over [0, 2pi)^2 with uniform weight"""
    a = sp.Rational(int(m.a.numerator), int(m.a.denominator))
    acc = {}
    for mo, cf in P.terms():
        e_Cx, e_Sx, e_c, e_s, e_iT, e_w1, e_w2 = mo
        assert e_w1 == 0 and e_w2 == 0
        if e_c % 2 or e_s % 2 or e_Sx % 2:
            continue
        th = sp.factorial2(e_c - 1)*sp.factorial2(e_s - 1)/sp.factorial2(e_c + e_s) if (e_c or e_s) else 1
        key = (e_Cx, e_Sx, e_iT)
        acc[key] = acc.get(key, 0) + sp.Rational(int(cf.numerator), int(cf.denominator))*th
    tt = sp.symbols('tt')
    lau = {}
    for (p, q2, n), cf in acc.items():
        expr = sp.expand(((tt - 1)/a)**p*(1 - ((tt - 1)/a)**2)**(q2//2))
        for (deg,), co in sp.Poly(expr, tt).terms():
            lau[deg - n] = lau.get(deg - n, 0) + cf*co
    return sp.nsimplify(sum(co*Tmoment(a, d) for d, co in lau.items() if co != 0))


def homogenized(m):
    """drift and diffusion of the m -> 0 limit from L_eff = -<A L^-1 A>, and the closed forms"""
    b, D = [], [[None]*2 for _ in range(2)]
    U = [Linv_odd(m, WV[i]) for i in range(2)]            # L^-1 (A X_i), A X_i = w_i
    for i in range(2):
        b.append(-gauss_avg(m, Aop(m, U[i])))
    for i in range(2):
        for j in range(2):
            D[i][j] = -(gauss_avg(m, WV[i]*U[j]) + gauss_avg(m, WV[j]*U[i]))/2
    return b, D


def analyse(m, kl_velocity="j", heat_force=None):
    """returns dict with kappa, sigma (KL rate), heat (finite part), Delta; options produce the ablation routes"""
    b, D = homogenized(m)
    rho = m.rho
    jv = [rho*b[i] - (m.dx(D[i][0]*rho) + m.dy(D[i][1]*rho)) for i in range(2)]
    if kl_velocity == "jkin":
        jv = [rho*m.h[i] for i in range(2)]
    MC0 = [[sum(m.M0[i][l]*m.C0[l][j] for l in range(2)) for j in range(2)] for i in range(2)]
    D0 = [[(MC0[i][j] + MC0[j][i])/2 for j in range(2)] for i in range(2)]
    detD0 = D0[0][0]*D0[1][1] - D0[0][1]*D0[1][0]
    adj = [[D[1][1], -D[0][1]], [-D[1][0], D[0][0]]]
    qf = sum(jv[i]*adj[i][j]*jv[j] for i in range(2) for j in range(2))
    Iint = qf*iT**(2 + m.rho_pow)*(1/detD0)            # j^T D^-1 j / rho, det D = T^2 det D0
    Zm = scal(iT, mm(mm(m.R, m.Z0), m.RT))
    Phi = sum(WV[i]*Zm[i][j]*WV[j] for i in range(2) for j in range(2))
    kappa = 2*(m.Q0[0][0]*m.Z0[0][0] + m.Q0[1][1]*m.Z0[1][1]) - sum(m.gam)
    F = m.F if heat_force is None else heat_force
    u = Linv_odd(m, Aop(m, Phi, F))
    Eint = -rho*gauss_avg(m, Aop(m, u, F))
    norm = torus_average(m, rho)
    sigma = torus_average(m, Iint)/norm
    heat = torus_average(m, Eint)/norm
    return dict(kappa=sp.Rational(int(kappa.numerator), int(kappa.denominator)), sigma=sp.nsimplify(sigma),
                heat=sp.nsimplify(heat), Delta=sp.nsimplify(heat - sigma), b=b, D=D, j=jv)


def instruction_model():
    R_ = sp.Rational
    return Model(a=R_(4, 5), lam=R_(1, 2), gam=(1, 3), tau=(1, 2), b=R_(1, 2), k=2, alpha=R_(1, 3), beta=R_(1, 2))


def check_force(m):
    """the constructed F equals the explicit force of the instruction (exact check in the ring, modulo the
    identities cos^2 + sin^2 = 1 and T * (1/T) = 1, tested at rational-angle points in high precision)"""
    import mpmath as mp
    mp.mp.dps = 40
    a, lam = mp.mpf(4)/5, mp.mpf(1)/2
    for xv, yv in [(mp.mpf('0.3'), mp.mpf('1.7')), (mp.mpf('2.2'), mp.mpf('5.1')), (mp.mpf('4.0'), mp.mpf('0.4'))]:
        Tv = 1 + a*mp.cos(xv); th = xv + yv + lam*mp.cos(xv); sx = mp.sin(xv)
        A1 = -mp.mpf(237)/65*sx + mp.mpf(2)/(3*Tv**2) - sx/(5*Tv)
        A2 = -1/(6*Tv**2) - 4*sx/(5*Tv)
        U = mp.mpf(72)/65*sx - mp.mpf(9)/13*Tv - mp.mpf(3)/26*Tv*sx - 1/(3*Tv**2)
        V = -mp.mpf(18)/65*sx + mp.mpf(15)/13*Tv - mp.mpf(6)/13*Tv*sx + 2*sx/(5*Tv)
        Fe = [A1 + U*mp.cos(2*th) + V*mp.sin(2*th), A2 + U*mp.sin(2*th) - V*mp.cos(2*th)]
        vals = {Cx: mp.cos(xv), Sx: sx, c: mp.cos(th), s: mp.sin(th), iT: 1/Tv}
        for i in range(2):
            Fc = sum(mp.mpf(int(cf.numerator))/int(cf.denominator)*mp.fprod(vals[g]**e for g, e in zip(GENS[:5], mo[:5]))
                     for mo, cf in m.F[i].terms())
            assert abs(Fc - Fe[i]) < mp.mpf(10)**-30, (i, Fc, Fe[i])


def main():
    m = instruction_model()
    check_force(m)
    assert m.sym["C0"] == sp.Matrix([[55, -6], [-6, 103]])/52
    res = analyse(m)
    # the homogenized drift/diffusion agree with the closed forms b = M F + (d_l M) Sigma_l, D = sym(M Sigma)
    bcf = [sum(m.M[i][j]*m.F[j] for j in range(2)) +
           sum((m.dx(m.M[i][j]) if l == 0 else m.dy(m.M[i][j]))*m.Sig[j][l] for j in range(2) for l in range(2))
           for i in range(2)]
    assert torus_average(m, (res["b"][0] - bcf[0])**2 + (res["b"][1] - bcf[1])**2) == 0
    divj = m.dx(res["j"][0]) + m.dy(res["j"][1])
    assert torus_average(m, divj**2) == 0, "rho = T^2/<T^2> must be stationary"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump({"Delta": str(res["Delta"])}, fh)
    print(f"kappa = {res['kappa']}   finite heat = {res['heat']}   KL rate = {res['sigma']}")
    print(f"Delta = {res['Delta']} = {float(res['Delta']):.15f}")


if __name__ == "__main__":
    main()
