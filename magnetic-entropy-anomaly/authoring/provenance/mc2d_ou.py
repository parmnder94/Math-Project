"""2D step MC with an exact Ornstein-Uhlenbeck velocity step (coefficients frozen over the step).
dV = -(G/m) V dt + sqrt(2 T Gamma)/m dW, G = Gamma(y) - J; stationary velocity covariance (T/m) I.
V' = E V + N, E = exp(-G dt/m), Cov N = (T/m)(I - E E^T); position X' = X + (V + V')/2 dt.
Estimator: crossing sum of (m|V|^2/2)(1/T_new - 1/T_old).  Usage: m P tau seed spm"""
import numpy as np, sys, json, time
m, P, tau, seed, spm = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
rng = np.random.default_rng(seed); dt = m / spm; h = dt / m
Tf = lambda x: np.where(x % (2 * np.pi) < np.pi, 2.0, 1.0)
x = rng.uniform(0, np.pi, P) + np.pi * (rng.uniform(size=P) > 1 / 3); y = rng.uniform(0, 2 * np.pi, P)
v = rng.standard_normal((2, P)) * np.sqrt(Tf(x) / m)
burn = int(2.0 / dt); steps = burn + int(tau / dt); acc = np.zeros(P); t0 = time.time()
for k in range(steps):
    T = Tf(x); c, s = np.cos(y), np.sin(y)
    A11, A22, A12 = -(c*c + 2*s*s) * h, -(s*s + 2*c*c) * h, -(-c*s) * h + 1.0 * h    # A = -G h, G = Gamma - J: G12 = Gamma12 + 1, G21 = Gamma12 - 1
    A21 = -(-c*s) * h - 1.0 * h
    tr2 = (A11 + A22) / 2; det = A11 * A22 - A12 * A21
    D = np.sqrt((tr2**2 - det).astype(complex)); ex = np.exp(tr2)
    sh = np.where(np.abs(D) > 1e-12, np.sinh(D) / np.where(np.abs(D) > 1e-12, D, 1), 1.0)
    ch = np.cosh(D)
    E11 = np.real(ex * (ch + sh * (A11 - tr2))); E22 = np.real(ex * (ch + sh * (A22 - tr2)))
    E12 = np.real(ex * sh * A12); E21 = np.real(ex * sh * A21)
    # Cov = (T/m)(I - E E^T)
    C11 = (T / m) * (1 - (E11**2 + E12**2)); C22 = (T / m) * (1 - (E21**2 + E22**2)); C12 = -(T / m) * (E11 * E21 + E12 * E22)
    l11 = np.sqrt(C11); l21 = C12 / l11; l22 = np.sqrt(np.maximum(C22 - l21**2, 0))
    z1, z2 = rng.standard_normal(P), rng.standard_normal(P)
    v0, v1 = v[0].copy(), v[1].copy()
    v[0] = E11 * v0 + E12 * v1 + l11 * z1
    v[1] = E21 * v0 + E22 * v1 + l21 * z1 + l22 * z2
    x = (x + 0.5 * (v0 + v[0]) * dt) % (2 * np.pi); y = (y + 0.5 * (v1 + v[1]) * dt) % (2 * np.pi)
    if k >= burn:
        Tn = Tf(x); acc += 0.5 * m * (v[0]**2 + v[1]**2) * (1.0 / Tn - 1.0 / T)
S = acc / tau
print(json.dumps(dict(m=m, spm=spm, sqrtm_S=np.sqrt(m) * S.mean(), se=np.sqrt(m) * S.std() / np.sqrt(P), secs=round(time.time() - t0))), flush=True)
