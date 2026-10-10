"""Monte Carlo of the literal underdamped SDE at small m, for the instruction model ('task') or the author's original
problem ('original'), estimating Sdot_m - kappa/m by (i) the definition and (ii) the exact identity m^(-1/2) E[A Phi].
Step: coefficients frozen at X over dt; the pair (w', integral of w' dt/sqrt(m)) in the reservoir frame is sampled
exactly (Van Loan matrix exponentials), so the only discretization error is the freezing.
Usage: mc_gen.py task|original m dt_over_m nparticles tau seed"""
import json, sys, time
import numpy as np, sympy as sp
from scipy.linalg import expm, cholesky
case, m, r, P, tau, seed = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
R_ = sp.Rational
PAR = {"task": dict(a=R_(4, 5), lam=R_(1, 2), gam=(1, 3), tau=(1, 2), b=R_(1, 2), k=2, alpha=R_(1, 3), beta=R_(1, 2)),
       "original": dict(a=R_(3, 5), lam=R_(2, 3), gam=(1, 2), tau=(1, 2), b=1, k=1, alpha=R_(2, 5), beta=R_(2, 3))}[case]
x_, y_ = sp.symbols('x y', real=True)
a, lam, gam, temps, b, k = PAR['a'], PAR['lam'], PAR['gam'], PAR['tau'], PAR['b'], PAR['k']
T_ = 1 + a*sp.cos(x_); th_ = x_ + y_ + lam*sp.cos(x_)
Rm = sp.Matrix([[sp.cos(th_), -sp.sin(th_)], [sp.sin(th_), sp.cos(th_)]])
B0s = sp.Matrix([[gam[0], b], [-b, gam[1]]]); Q0s = sp.diag(gam[0]*temps[0], gam[1]*temps[1])
W0s = sp.diag(R_(gam[0], temps[0]), R_(gam[1], temps[1]))
cc = sp.symbols('c0:3'); C0s = sp.Matrix([[cc[0], cc[1]], [cc[1], cc[2]]]); C0s = C0s.subs(sp.solve(list(B0s*C0s + C0s*B0s.T - 2*Q0s), cc))
zz = sp.symbols('z0:3'); Z0s = sp.Matrix([[zz[0], zz[1]], [zz[1], zz[2]]]); Z0s = Z0s.subs(sp.solve(list(B0s.T*Z0s + Z0s*B0s - W0s), zz))
Sig = T_*Rm*C0s*Rm.T; Bm = Rm*B0s*Rm.T
h = sp.Matrix([PAR['alpha']/T_**k, PAR['beta']*T_**(1 - k)*sp.diff(T_, x_)])
F_ = sp.Matrix([sp.diff(Sig[i, 0], x_) + sp.diff(Sig[i, 1], y_) for i in range(2)]) + Sig*sp.Matrix([k*sp.diff(T_, x_)/T_, 0]) + Bm*h
Ff = sp.lambdify((x_, y_), [F_[0], F_[1]], 'numpy')
B0, Q0, W0 = [np.array(M.tolist(), float) for M in (B0s, Q0s, W0s)]
C0, Z0 = np.array(C0s.tolist(), float), np.array(Z0s.tolist(), float)
kappa = 2*np.trace(Q0@Z0) - sum(gam)
a, lam = float(a), float(lam)
dt = r*m; rng = np.random.default_rng(seed)
A = np.zeros((4, 4)); A[:2, :2] = -B0/m; A[2:, :2] = np.eye(2)/np.sqrt(m)
Phi = expm(A*dt)
Aug = np.zeros((8, 8)); Aug[:4, :4] = A; Aug[:4, 4:] = np.eye(4); Psi = expm(Aug*dt)[:4, 4:]
Gq = np.zeros((4, 4)); Gq[:2, :2] = 2*Q0/m
VL = np.zeros((8, 8)); VL[:4, :4] = -A; VL[:4, 4:] = Gq; VL[4:, 4:] = A.T
E = expm(VL*dt); Kc = E[4:, 4:].T @ E[:4, 4:]; Lc = cholesky((Kc + Kc.T)/2, lower=True)
J = np.array([[0, -1], [1, 0]]); JZ = J@Z0 - Z0@J
def to_frame(c, s, v): return np.stack([c*v[0] + s*v[1], -s*v[0] + c*v[1]])
def from_frame(c, s, v): return np.stack([c*v[0] - s*v[1], s*v[0] + c*v[1]])
x = np.empty(0)
while x.size < P:
    cand = rng.uniform(0, 2*np.pi, 2*P); keep = rng.uniform(0, (1 + a)**k, 2*P) < (1 + a*np.cos(cand))**k
    x = np.concatenate([x, cand[keep]])
x = x[:P]; y = rng.uniform(0, 2*np.pi, P)
T = 1 + a*np.cos(x); th = x + y + lam*np.cos(x)
w = from_frame(np.cos(th), np.sin(th), np.sqrt(T)*(np.linalg.cholesky(C0)@rng.standard_normal((2, P))))
burn = int(2.0/dt); nsteps = burn + int(tau/dt)
accA = np.zeros(P); accH = np.zeros(P); t0 = time.time()
for it in range(nsteps):
    T = 1 + a*np.cos(x); th = x + y + lam*np.cos(x); c, s = np.cos(th), np.sin(th)
    F = np.array(Ff(x, y))
    wp = to_frame(c, s, w); Fp = to_frame(c, s, F)
    if it >= burn:
        q0 = np.einsum('in,ij,jn->n', wp, Z0, wp); qJ = np.einsum('in,ij,jn->n', wp, JZ, wp)
        nx = 1 - lam*np.sin(x); Tx = -a*np.sin(x)
        Zw = from_frame(c, s, Z0@wp)/T
        accA += w[0]*(qJ*nx/T - q0*Tx/T**2) + w[1]*qJ/T + 2*(F[0]*Zw[0] + F[1]*Zw[1])
        accH += (np.einsum('in,ij,jn->n', wp, W0, wp)/T - sum(gam))/m
    z = Phi@np.concatenate([wp, np.zeros((2, P))]) + Psi@np.concatenate([Fp/np.sqrt(m), np.zeros((2, P))]) \
        + np.sqrt(T)*(Lc@rng.standard_normal((4, P)))
    w = from_frame(c, s, z[:2]); dX = from_frame(c, s, z[2:])
    x = (x + dX[0]) % (2*np.pi); y = (y + dX[1]) % (2*np.pi)
n = nsteps - burn
eA = accA/n/np.sqrt(m); eH = accH/n - kappa/m
print(json.dumps(dict(case=case, m=m, dt_over_m=r, P=P, tau=tau, seed=seed, kappa=kappa,
                      E_identity=eA.mean(), se_identity=eA.std()/np.sqrt(P), E_definition=eH.mean(),
                      se_definition=eH.std()/np.sqrt(P), secs=round(time.time() - t0))), flush=True)
