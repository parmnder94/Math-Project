"""sim2 with the one-period propagator U_F and Lam = i log U_F computed in 40-digit arithmetic (Taylor integrator)."""
import numpy as np, mpmath as mp, sys, sim2
mp.mp.dps = 40
s3m = mp.sqrt(3)
def mpmat(A): return mp.matrix([[mp.mpc(complex(A[i, j]).real, complex(A[i, j]).imag) for j in range(4)] for i in range(4)])
JxM = mp.matrix([[0, s3m/2, 0, 0], [s3m/2, 0, 1, 0], [0, 1, 0, s3m/2], [0, 0, s3m/2, 0]])
PM = mp.diag([0, 1, 1, 0]); QM = mp.eye(4) - PM
XM = JxM*JxM - (7*PM + 3*QM)/4
def UF_mp(eps, ph, nsteps=80, order=36):
    eps = mp.mpf(eps); ph = mp.mpf(ph)
    a = eps**2 + mp.mpf(4)/3*eps**4 + mp.mpf(148)/45*eps**6; b = eps**2 + mp.mpf(5)/3*eps**4 + mp.mpf(167)/36*eps**6
    H0 = 2*QM + a*(PM - QM)
    h = 2*mp.pi/nsteps; U = mp.eye(4)
    for st in range(nsteps):
        t0 = st*h
        # Taylor coefficients of H(t0 + tau)
        Hk = []
        for k in range(order + 1):
            d1 = mp.cos(t0 + ph + k*mp.pi/2)/mp.factorial(k)                # d^k/dt^k cos(t+ph) / k!
            d2 = 2**k*mp.cos(2*(t0 + ph) + k*mp.pi/2)/mp.factorial(k)
            M = 2*eps*d1*JxM - 2*b*d2*XM
            if k == 0: M = M + H0
            Hk.append(M)
        Uk = [U]
        for k in range(order):
            S = mp.zeros(4)
            for j in range(k + 1): S += Hk[j]*Uk[k - j]
            Uk.append(-1j*S/(k + 1))
        U = mp.zeros(4)
        for k in range(order, -1, -1): U = U*h + Uk[k] if False else U
        # Horner in tau = h
        U = Uk[order]
        for k in range(order - 1, -1, -1): U = U*h + Uk[k]
    return U
_floquet_double = sim2.floquet
def floquet_mp(eps, ph):
    _, L = _floquet_double(eps, ph)                  # averaging map in double precision
    U = UF_mp(eps, ph)
    E, V = mp.eig(U)
    Lam = V*mp.diag([1j*mp.log(e) for e in E])*mp.inverse(V)
    Lam = np.array([[complex(Lam[i, j]) for j in range(4)] for i in range(4)])
    return Lam, L
sim2.floquet = floquet_mp
if __name__ == "__main__":
    for N in [int(float(v)) for v in sys.argv[1:]]:
        p, eps = sim2.prob(N); print(f"N={N:.0e} eps={eps:.5f} p={p:.10f} p-p*={p-1467167323933/5157470703125:+.3e}", flush=True)
