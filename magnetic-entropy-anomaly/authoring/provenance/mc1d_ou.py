"""1D step MC with an exact OU velocity step (T frozen over the step) and exact integrated displacement.
Estimator: crossing sum of (m V^2/2)(1/T_new - 1/T_old).  Usage: m P tau seed spm"""
import numpy as np, sys, json, time
m, P, tau, seed, spm = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
rng = np.random.default_rng(seed); dt = m / spm; a = np.exp(-dt / m)
Tf = lambda x: np.where(x % (2*np.pi) < np.pi, 2.0, 1.0)
x = rng.uniform(0, np.pi, P) + np.pi * (rng.uniform(size=P) > 1/3); v = rng.standard_normal(P) * np.sqrt(Tf(x) / m)
# joint Gaussian of (v_{n+1}, displacement) for OU with gamma/m = 1/m over dt: exact covariance
g = 1 / m; h = dt
var_v = (1 - a**2) / (2 * g)                     # times 2T/m^2 * ... handle by scaling below
burn = int(2.0 / dt); steps = burn + int(tau / dt); acc = np.zeros(P); t0 = time.time()
# For dV = -g V dt + s dW with s^2 = 2 T g / m: V' = a V + N1, X' = X + V (1-a)/g + N2 with
# Var N1 = s^2 (1-a^2)/(2g), Var N2 = s^2/g^2 [h - 2(1-a)/g + (1-a^2)/(2g)], Cov = s^2/(2 g^2) (1-a)^2
c11 = (1 - a**2) / (2 * g); c22 = (h - 2 * (1 - a) / g + (1 - a**2) / (2 * g)) / g**2; c12 = (1 - a)**2 / (2 * g**2)
Lc = np.linalg.cholesky(np.array([[c11, c12], [c12, c22]]))
for k in range(steps):
    T = Tf(x); s = np.sqrt(2 * T * g / m)
    z = rng.standard_normal((2, P)); n1 = s * (Lc[0, 0] * z[0]); n2 = s * (Lc[1, 0] * z[0] + Lc[1, 1] * z[1])
    x = (x + v * (1 - a) / g + n2) % (2 * np.pi); v = a * v + n1
    if k >= burn:
        Tn = Tf(x); acc += 0.5 * m * v * v * (1.0 / Tn - 1.0 / T)
S = acc / tau
print(json.dumps(dict(m=m, spm=spm, sqrtm_S=np.sqrt(m) * S.mean(), se=np.sqrt(m) * S.std() / np.sqrt(P), secs=round(time.time() - t0))), flush=True)
