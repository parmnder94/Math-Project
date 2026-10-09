"""Monte Carlo of the underdamped dynamics: estimates S_dot_m = -<S>/(2 sqrt m) (exact stationary identity).
Usage: python3 mc.py <m> <particles> <time> <seed> [steps per unit m, default 40]"""
import numpy as np, sys, time
m = float(sys.argv[1]); P = int(sys.argv[2]); tau = float(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed); dt = m / float(sys.argv[5]) if len(sys.argv) > 5 else m / 40; eps = np.sqrt(m)
x = rng.uniform(0, 2*np.pi, P); y = rng.uniform(0, 2*np.pi, P)
T = lambda x: 1 + 0.6*np.cos(x)
u = rng.standard_normal((2, P)) * np.sqrt(T(x))       # u = sqrt(m) V
def coeffs(x, y):
    c, s = np.cos(y), np.sin(y)
    # Gamma = R diag(1,2) R^T ; sqrt(2 T Gamma) = R diag(sqrt(2T), sqrt(4T)) R^T
    G11 = c*c + 2*s*s; G22 = s*s + 2*c*c; G12 = -s*c + 2*s*c
    t = T(x); a, b = np.sqrt(2*t), np.sqrt(4*t)
    S11 = c*c*a + s*s*b; S22 = s*s*a + c*c*b; S12 = c*s*(b - a)
    return G11, G12, G22, S11, S12, S22, t
acc = 0.0; n = 0; burn = int(2.0 / dt); steps = burn + int(tau / dt)
t0 = time.time()
for k in range(steps):
    G11, G12, G22, S11, S12, S22, t = coeffs(x, y)
    g1 = -0.6*np.sin(x)
    if k >= burn:
        gu = g1*u[0]; S = gu/(t*t)*(u[0]**2 + u[1]**2 - 4*t); acc += S.mean(); n += 1
    # du = [-(Gamma - J)u/eps^2 + g/eps] dt + sqrt(2 T Gamma)/eps dW ; J u = (-u2, u1)
    du1 = (-(G11*u[0] + G12*u[1]) - u[1]) / m * dt + g1/eps*dt
    du2 = (-(G12*u[0] + G22*u[1]) + u[0]) / m * dt
    dW = rng.standard_normal((2, P)) * np.sqrt(dt)
    du1 += (S11*dW[0] + S12*dW[1]) / eps; du2 += (S12*dW[0] + S22*dW[1]) / eps
    x = (x + u[0]/eps*dt) % (2*np.pi); y = (y + u[1]/eps*dt) % (2*np.pi)
    u = u + np.vstack([du1, du2])
print(f"m={m} P={P} tau={tau} seed={seed}  Sdot_m = {-acc/n/(2*eps):.5f}  ({time.time()-t0:.0f}s)", flush=True)
