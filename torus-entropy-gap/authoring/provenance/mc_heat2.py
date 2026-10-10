"""Monte Carlo of the literal underdamped SDE (instruction parameters) at small m.
Estimates  Sdot_m - kappa/m  two ways:
  (i)  the definition, (1/m) E[w^T W w / T - sum(gam)] - kappa/m        (W = R diag(gam/tau) R^T, w = sqrt(m) V)
  (ii) the exact stationary identity m^(-1/2) E[A Phi], Phi = w^T Z w, Z = R Z0 R^T / T, B0^T Z0 + Z0 B0 = W0
Step: with coefficients frozen at X over dt, the pair (w', integral of w' dt/sqrt(m)) in the reservoir frame is a
Gaussian linear SDE; it is sampled exactly (Van Loan matrix exponentials), so the only error is the freezing (O(dt)).
Usage: mc_heat.py m dt_over_m nparticles tau seed"""
import numpy as np, sys, json, time
from scipy.linalg import expm, cholesky
m, r, P, tau, seed = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5])
dt = r*m; rng = np.random.default_rng(seed)
a, lam, b, gam, temps = 0.8, 0.5, 0.5, np.array([1.0, 3.0]), np.array([1.0, 2.0])
B0 = np.array([[gam[0], b], [-b, gam[1]]]); Q0 = np.diag(gam*temps); W0 = np.diag(gam/temps)
C0 = np.array([[55/52, -3/26], [-3/26, 103/52]]); assert np.allclose(B0@C0 + C0@B0.T, 2*Q0)
Z0 = np.array([[101/208, -3/104], [-3/104, 53/208]]); assert np.allclose(B0.T@Z0 + Z0@B0, W0)
kappa = 2*np.trace(Q0@Z0) - gam.sum(); assert abs(kappa - 3/104) < 1e-14
M0 = np.linalg.inv(B0)
# augmented linear system z = (w', xi): dz = (A z + f) dt + noise, A = [[-B0/m, 0], [I/sqrt(m), 0]], noise cov (2/m) T Q0 on w'
A = np.zeros((4, 4)); A[:2, :2] = -B0/m; A[2:, :2] = np.eye(2)/np.sqrt(m)
Phi = expm(A*dt)
Aug = np.zeros((8, 8)); Aug[:4, :4] = A; Aug[:4, 4:] = np.eye(4)
Psi = expm(Aug*dt)[:4, 4:]                                     # int_0^dt e^{A s} ds
Gq = np.zeros((4, 4)); Gq[:2, :2] = 2*Q0/m
VL = np.zeros((8, 8)); VL[:4, :4] = -A; VL[:4, 4:] = Gq; VL[4:, 4:] = A.T
E = expm(VL*dt); Kc = E[4:, 4:].T @ E[:4, 4:]                   # Van Loan: covariance for T = 1
Kc = (Kc + Kc.T)/2; Lc = cholesky(Kc + 1e-300*np.eye(4), lower=True)
J = np.array([[0, -1], [1, 0]])
def fields(x, y):
    T = 1 + a*np.cos(x); th = x + y + lam*np.cos(x); s = np.sin(x)
    A1 = -237/65*s + 2/(3*T**2) - s/(5*T); A2 = -1/(6*T**2) - 4*s/(5*T)
    U = 72/65*s - 9/13*T - 3/26*T*s - 1/(3*T**2); V = -18/65*s + 15/13*T - 6/13*T*s + 2*s/(5*T)
    c2, s2 = np.cos(2*th), np.sin(2*th)
    F = np.stack([A1 + U*c2 + V*s2, A2 + U*s2 - V*c2])
    return T, th, F
def rot(th): return np.cos(th), np.sin(th)
def to_frame(c, s, v): return np.stack([c*v[0] + s*v[1], -s*v[0] + c*v[1]])      # R^T v
def from_frame(c, s, v): return np.stack([c*v[0] - s*v[1], s*v[0] + c*v[1]])     # R v
def APhi(x, y, w, T, th, F):
    c, s = rot(th)
    wp = to_frame(c, s, w)                                   # w' = R^T w
    # Z = R Z0 R^T / T ; dZ/dth = R [J, Z0] R^T / T ; grad th = (1 - lam sin x, 1) ; grad T = (-a sin x, 0)
    JZ = J@Z0 - Z0@J
    q0 = np.einsum('in,ij,jn->n', wp, Z0, wp); qJ = np.einsum('in,ij,jn->n', wp, JZ, wp)
    nx = 1 - lam*np.sin(x); Tx = -a*np.sin(x)
    dPhi_dx = qJ*nx/T - q0*Tx/T**2; dPhi_dy = qJ/T
    Zw = from_frame(c, s, Z0@wp)/T
    return w[0]*dPhi_dx + w[1]*dPhi_dy + 2*(F[0]*Zw[0] + F[1]*Zw[1])
def heat_obs(w, T, th):
    c, s = rot(th); wp = to_frame(c, s, w)
    return (np.einsum('in,ij,jn->n', wp, W0, wp)/T - gam.sum())/m
# initial state: x from rho ~ T^2 (rejection), y uniform, w from the local Gaussian
x = np.empty(0)
while x.size < P:
    cand = rng.uniform(0, 2*np.pi, 2*P); keep = rng.uniform(0, 1.8**2, 2*P) < (1 + a*np.cos(cand))**2
    x = np.concatenate([x, cand[keep]])
x = x[:P]; y = rng.uniform(0, 2*np.pi, P)
T, th, F = fields(x, y)
c, s = rot(th); w = from_frame(c, s, np.sqrt(T)*(np.linalg.cholesky(C0)@rng.standard_normal((2, P))))
burn = int(2.0/dt); nsteps = burn + int(tau/dt)
accA = np.zeros(P); accH = np.zeros(P); t0 = time.time()
for k in range(nsteps):
    T, th, F = fields(x, y)
    if k >= burn:
        accA += APhi(x, y, w, T, th, F); accH += heat_obs(w, T, th)
    c, s = rot(th)
    wp = to_frame(c, s, w); Fp = to_frame(c, s, F)
    z0 = np.concatenate([wp, np.zeros((2, P))])
    f = np.concatenate([Fp/np.sqrt(m), np.zeros((2, P))])
    z = Phi@z0 + Psi@f + np.sqrt(T)*(Lc@rng.standard_normal((4, P)))
    w = from_frame(c, s, z[:2]); dX = from_frame(c, s, z[2:])
    x = (x + dX[0]) % (2*np.pi); y = (y + dX[1]) % (2*np.pi)
n = nsteps - burn
estA = accA/n/np.sqrt(m); estH = accH/n - kappa/m
print(json.dumps(dict(m=m, dt_over_m=r, P=P, tau=tau, seed=seed,
                      E_identity=estA.mean(), se_identity=estA.std()/np.sqrt(P),
                      E_definition=estH.mean(), se_definition=estH.std()/np.sqrt(P), secs=round(time.time() - t0))), flush=True)
