"""Monte Carlo of the literal underdamped SDE at small m (fast numpy version, generic parameters).
Cases: 'task' (instruction), 'mild' (instruction model with T = 1 + (3/5)cos x), 'original' (author-supplied problem).
The force is built from F = div Sigma + Sigma grad(k log T) + B h and checked against an independent SymPy evaluation
at start-up (for 'task' this is the explicit force of the instruction, verified in solution/solve.py).
Estimates Sdot_m - kappa/m by (i) the definition and (ii) the exact identity m^(-1/2) E[A Phi].
Step: coefficients frozen at X over dt; (w', integral of w' dt/sqrt(m)) sampled exactly in the reservoir frame.
With scheme 'mid' the step is repeated with the coefficients evaluated at the predicted midpoint X + dX/2 (same
random numbers), removing the first-order freezing error, which otherwise biases the O(sqrt(m)) correlations that
carry the finite heat by O(dt/m) even as m -> 0.
Usage: mc_fast.py case m dt_over_m nparticles tau seed [left|mid]"""
import json, sys, time
import numpy as np
from scipy.linalg import expm, cholesky, solve_continuous_lyapunov
case, m, r, P, tau, seed = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])
scheme = sys.argv[7] if len(sys.argv) > 7 else 'left'
PAR = {"task": dict(a=0.8, lam=0.5, gam=(1, 3), temps=(1, 2), b=0.5, k=2, alpha=1/3, beta=0.5),
       "mild": dict(a=0.6, lam=0.5, gam=(1, 3), temps=(1, 2), b=0.5, k=2, alpha=1/3, beta=0.5),
       "original": dict(a=0.6, lam=2/3, gam=(1, 2), temps=(1, 2), b=1.0, k=1, alpha=0.4, beta=2/3)}[case]
a, lam, b, k, alpha, beta = (PAR[q] for q in ("a", "lam", "b", "k", "alpha", "beta"))
gam, temps = np.array(PAR["gam"], float), np.array(PAR["temps"], float)
B0 = np.array([[gam[0], b], [-b, gam[1]]]); Q0 = np.diag(gam*temps); W0 = np.diag(gam/temps)
C0 = solve_continuous_lyapunov(B0, 2*Q0)                      # B0 C0 + C0 B0^T = 2 Q0
Z0 = solve_continuous_lyapunov(B0.T, W0)                      # B0^T Z0 + Z0 B0 = W0
kappa = 2*np.trace(Q0@Z0) - gam.sum()
cbar, d, q = (C0[0, 0] + C0[1, 1])/2, (C0[0, 0] - C0[1, 1])/2, C0[0, 1]
J = np.array([[0, -1], [1, 0]]); JZ = J@Z0 - Z0@J
def to_frame(c, s, v): return np.stack([c*v[0] + s*v[1], -s*v[0] + c*v[1]])
def from_frame(c, s, v): return np.stack([c*v[0] - s*v[1], s*v[0] + c*v[1]])
def force(x, y):
    T = 1 + a*np.cos(x); Tp = -a*np.sin(x); th = x + y + lam*np.cos(x); thx = 1 - lam*np.sin(x)
    c2, s2 = np.cos(2*th), np.sin(2*th)
    Pp, Ss = d*c2 - q*s2, d*s2 + q*c2
    div1 = Tp*(cbar + Pp) - 2*T*Ss*thx + 2*T*Pp
    div2 = Tp*Ss + 2*T*Pp*thx + 2*T*Ss
    g1, g2 = k*Tp*(cbar + Pp), k*Tp*Ss
    h = np.stack([alpha/T**k, beta*Tp*T**(1 - k)])
    c, s = np.cos(th), np.sin(th)
    Bh = from_frame(c, s, B0@to_frame(c, s, h))
    return np.stack([div1 + g1 + Bh[0], div2 + g2 + Bh[1]]), T, th
# start-up check of the force against an independent SymPy construction
import sympy as sp
xs, ys = sp.symbols('x y'); Ts = 1 + sp.nsimplify(a)*sp.cos(xs); ths = xs + ys + sp.nsimplify(lam)*sp.cos(xs)
Rs = sp.Matrix([[sp.cos(ths), -sp.sin(ths)], [sp.sin(ths), sp.cos(ths)]])
Sig = Ts*Rs*sp.Matrix(C0.tolist())*Rs.T; Bs = Rs*sp.Matrix(B0.tolist())*Rs.T
hs = sp.Matrix([sp.nsimplify(alpha)/Ts**k, sp.nsimplify(beta)*sp.diff(Ts, xs)*Ts**(1 - k)])
Fs = sp.Matrix([sp.diff(Sig[i, 0], xs) + sp.diff(Sig[i, 1], ys) for i in range(2)]) + Sig*sp.Matrix([k*sp.diff(Ts, xs)/Ts, 0]) + Bs*hs
Ff = sp.lambdify((xs, ys), [Fs[0], Fs[1]], "numpy")
tx, ty = np.array([0.3, 2.2, 4.0]), np.array([1.7, 5.1, 0.4])
assert np.allclose(np.array(Ff(tx, ty)), force(tx, ty)[0], atol=1e-11)
dt = r*m; rng = np.random.default_rng(seed)
A = np.zeros((4, 4)); A[:2, :2] = -B0/m; A[2:, :2] = np.eye(2)/np.sqrt(m)
Phi = expm(A*dt)
Aug = np.zeros((8, 8)); Aug[:4, :4] = A; Aug[:4, 4:] = np.eye(4); Psi = expm(Aug*dt)[:4, 4:]
Gq = np.zeros((4, 4)); Gq[:2, :2] = 2*Q0/m
VL = np.zeros((8, 8)); VL[:4, :4] = -A; VL[:4, 4:] = Gq; VL[4:, 4:] = A.T
E = expm(VL*dt); Kc = E[4:, 4:].T @ E[:4, 4:]; Lc = cholesky((Kc + Kc.T)/2, lower=True)
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
    F, T, th = force(x, y); c, s = np.cos(th), np.sin(th)
    wp = to_frame(c, s, w); Fp = to_frame(c, s, F)
    if it >= burn:
        q0 = np.einsum('in,ij,jn->n', wp, Z0, wp); qJ = np.einsum('in,ij,jn->n', wp, JZ, wp)
        nx = 1 - lam*np.sin(x); Tx = -a*np.sin(x)
        Zw = from_frame(c, s, Z0@wp)/T
        accA += w[0]*(qJ*nx/T - q0*Tx/T**2) + w[1]*qJ/T + 2*(F[0]*Zw[0] + F[1]*Zw[1])
        accH += (np.einsum('in,ij,jn->n', wp, W0, wp)/T - gam.sum())/m
    xi = rng.standard_normal((4, P))
    z = Phi@np.concatenate([wp, np.zeros((2, P))]) + Psi@np.concatenate([Fp/np.sqrt(m), np.zeros((2, P))]) + np.sqrt(T)*(Lc@xi)
    dX = from_frame(c, s, z[2:])
    if scheme == 'mid':
        xm, ym = x + dX[0]/2, y + dX[1]/2
        Fm, Tm, thm = force(xm, ym); cm, sm = np.cos(thm), np.sin(thm)
        z = Phi@np.concatenate([to_frame(cm, sm, w), np.zeros((2, P))]) + Psi@np.concatenate([to_frame(cm, sm, Fm)/np.sqrt(m), np.zeros((2, P))]) \
            + np.sqrt(Tm)*(Lc@xi)
        c, s = cm, sm
        dX = from_frame(c, s, z[2:])
    w = from_frame(c, s, z[:2])
    x = (x + dX[0]) % (2*np.pi); y = (y + dX[1]) % (2*np.pi)
n = nsteps - burn
eA = accA/n/np.sqrt(m); eH = accH/n - kappa/m
print(json.dumps(dict(case=case, scheme=scheme, m=m, dt_over_m=r, P=P, tau=tau, seed=seed, kappa=kappa,
                      E_identity=eA.mean(), se_identity=eA.std()/np.sqrt(P), E_definition=eH.mean(),
                      se_definition=eH.std()/np.sqrt(P), secs=round(time.time() - t0))), flush=True)
