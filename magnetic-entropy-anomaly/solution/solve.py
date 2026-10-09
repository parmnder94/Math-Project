"""Reference solution for magnetic-entropy-anomaly (temperature-step version).

Model: dX = V dt, m dV = [-Gamma(y) V + b J V] dt + sqrt(2 T Gamma) dW on the torus [0, 2pi)^2, with T = Th = 2 on
0 < x < pi and T = Tc = 1 on pi < x < 2pi, Gamma(y) = R(y) diag(g1, g2) R(y)^T, (g1, g2) = (1, 2), b = 1.
Target: K = lim_{m->0} sqrt(m) S_dot_m.

1. Bulk.  In each stripe T is constant, the fast velocity is locally Maxwellian and the limiting position process is
   the overdamped diffusion with mobility (Gamma - bJ)^-1; its stationary density is uniform in each stripe and carries
   no current.  Across a step the normal current vanishes and the exact momentum balance of the kinetic layer,
   d(P_11 gam_22 - P_12 gam_12)/dxi = 0 with gam = Gamma - bJ, makes the pressure rho T continuous (P_12 = 0 in the bulk).
   Hence rho_h Th = rho_c Tc = P and, with mass 1 on the 2pi x 2pi torus, P = 1/(2 pi^2 (1/Th + 1/Tc)) = 1/(3 pi^2).
2. Heat.  No work is done on the particle (the Lorentz force is perpendicular to V), so the heat released into the
   medium at x per unit time is (Gamma:Pi - T tr(Gamma) rho)/m, Pi = m <V V^T> rho.  It vanishes in the Maxwellian
   bulk; it lives in kinetic layers of width sqrt(m) at x = 0 and x = pi.  With xi = x/sqrt(m) and u = sqrt(m) V the
   layer density solves  u1 d_xi p = div_u[(Gamma - bJ) u p + T(xi) Gamma grad_u p]  (T jumps at xi = 0), and
   sqrt(m) S_dot_m -> sum over both lines of  int dy int dxi (Gamma:Pi - T tr(Gamma) rho)/T.
3. Energy balance in the layer: d/dxi (q/2) = -(Gamma:Pi - T tr(Gamma) rho), q = int u1 |u|^2 p du, and q -> 0 in
   the bulk on both sides.  So each line contributes int dy (q(xi=0)/2)(1/T_right - 1/T_left).
4. Symmetry.  Reflecting x (and u1) maps the line x = pi onto x = 0 with (theta, b) -> (-theta, -b), and reflecting
   u2 maps (theta, b) -> (-theta, -b) at fixed line; so q at x = pi is minus q at x = 0, and with q_a(theta) the
   interface flux at x = 0 (cold side xi < 0, hot side xi > 0, normalised to P = 1),
   K = P * 2pi * < -q_a/2 >_theta = -<q_a>_theta / (3 pi).
5. Layer solver (Hermite moment method).  p = phi_Tr(u) g, g = sum c_{n1 n2} h_n1(s1) h_n2(s2), s = u/sqrt(Tr),
   normalised Hermite functions, Tr = 3/2.  The layer equation becomes sqrt(Tr) S1 c' = A_T c with S1 the matrix of
   s1 (Jacobi matrix, invertible for even K1) and A_T the Galerkin matrix of the operator
   (1-t) tr(Gamma) g + (gam s).grad g - 2t (Gamma s).grad g + t Gamma:hess g + (t-1) s^T Gamma s g,  t = T/Tr.
   Bounded solutions on each half line = the Maxwellian at that temperature plus the stable invariant subspace of
   (sqrt(Tr) S1)^-1 A_T (ordered real Schur form, orthonormal and well conditioned).  Continuity of p at xi = 0 fixes
   the amplitudes; rho T = 1 on both sides is consistent (matching residual ~1e-15).  The transverse order K2 = 16 is
   converged to 1e-13; the normal order K1 converges like K1^-2.87 (velocity singularity at u1 = 0 on the step), so
   the results at K1 = 40 and 80 are Richardson-extrapolated (error about 5e-6 relative; the truth uses 120 and 160).  Gamma has period pi in y, and q_a(theta) is smooth, so
   the trapezoid rule with 8 angles is exact to 1e-10.
"""
import json, os
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sps
from math import lgamma

TH, TC, G1, G2, BFIELD = 2.0, 1.0, 1.0, 2.0, 1.0
TREF, K2, NTHETA = 1.5, 16, 8
K1_PAIR, P_EXP = (40, 80), 2.87
OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")


def maxwellian_coeffs(K, t):
    """Coefficients of phi_t / phi_1 in normalised Hermite functions h_n (1D): c_2k = ((t-1)/2)^k sqrt((2k)!)/k!."""
    c = np.zeros(K)
    for n in range(0, K, 2):
        k = n // 2
        if t == 1:
            c[n] = 1.0 if k == 0 else 0.0
        else:
            c[n] = np.sign(t - 1) ** k * np.exp(k * np.log(abs(t - 1) / 2) - lgamma(k + 1) + 0.5 * lgamma(n + 1))
    return c


def ladder(K):
    n = np.arange(K)
    S = np.diag(np.sqrt(n[1:]), -1) + np.diag(np.sqrt(n[1:]), 1)     # s h_n = sqrt(n+1) h_{n+1} + sqrt(n) h_{n-1}
    D = np.diag(np.sqrt(n[1:]), 1)                                     # h_n' = sqrt(n) h_{n-1}
    return S, D


def layer_operator(K1, K2, G, b, t):
    pad = 4
    Ka, Kb = K1 + pad, K2 + pad
    S1, D1 = ladder(Ka)
    S2, D2 = ladder(Kb)
    s = [sps.csr_matrix(np.kron(S1, np.eye(Kb))), sps.csr_matrix(np.kron(np.eye(Ka), S2))]
    d = [sps.csr_matrix(np.kron(D1, np.eye(Kb))), sps.csr_matrix(np.kron(np.eye(Ka), D2))]
    gam = G - b * np.array([[0.0, -1.0], [1.0, 0.0]])
    L = (1 - t) * np.trace(G) * sps.identity(Ka * Kb)
    for i in range(2):
        for j in range(2):
            L = L + (gam[i, j] - 2 * t * G[i, j]) * (s[j] @ d[i]) + t * G[i, j] * (d[i] @ d[j]) \
                + (t - 1) * G[i, j] * (s[i] @ s[j])
    keep = np.array([a * Kb + c for a in range(K1) for c in range(K2)])
    return L.tocsr()[keep][:, keep].toarray()


def interface_flux(theta, K1, K2=K2, b=BFIELD, g=(G1, G2), Tleft=TC, Tright=TH, Tr=TREF):
    """q(xi = 0) = int u1 |u|^2 p du for the layer with Tleft on xi < 0, Tright on xi > 0, rho T = 1 in both bulks."""
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    G = R @ np.diag(g) @ R.T
    S1, _ = ladder(K1)
    Sinv = np.kron(np.linalg.inv(np.sqrt(Tr) * S1), np.eye(K2))
    parts = []
    for T, side in ((Tleft, -1), (Tright, +1)):
        B = Sinv @ layer_operator(K1, K2, G, b, T / Tr)
        keep = (lambda re: re < -1e-6) if side > 0 else (lambda re: re > 1e-6)
        try:
            _, Z, k = sla.schur(B, output="real", sort=lambda re, im: keep(re))
        except np.linalg.LinAlgError:
            _, Z, k = sla.schur(B.astype(complex), output="complex", sort=lambda z: keep(z.real))
        parts.append((Z[:, :k], np.kron(maxwellian_coeffs(K1, T / Tr), maxwellian_coeffs(K2, T / Tr)) / T))
    (Zl, Ml), (Zr, Mr) = parts
    A = np.hstack([Zl, -Zr])
    sol = np.linalg.lstsq(A, (Mr - Ml).astype(A.dtype), rcond=None)[0]
    resid = np.abs(A @ sol - (Mr - Ml)).max()
    assert resid < 1e-10, resid
    c = np.real(Zr @ sol[Zl.shape[1]:] + Mr).reshape(K1, K2)
    # s1 |s|^2 = He3(s1) + 3 He1(s1) + He1(s1)(He2(s2) + 1)  ->  sqrt6 c30 + sqrt2 c12 + 4 c10
    return (np.sqrt(6) * c[3, 0] + np.sqrt(2) * c[1, 2] + 4 * c[1, 0]) * Tr ** 1.5


def K_value(K1, ntheta=NTHETA, **kw):
    thetas = np.arange(ntheta) * np.pi / ntheta
    qa = np.array([interface_flux(th, K1, **kw) for th in thetas])
    P = 1.0 / (2 * np.pi ** 2 * (1 / TH + 1 / TC))
    # line x = 0 (cold left, hot right) gives P (q_a/2)(1/Th - 1/Tc) per unit length; line x = pi gives the same
    # (q_b = -q_a with the sides swapped).  K = 2 pi P (1/Tc - 1/Th) < -q_a >_theta.
    return 2 * np.pi * P * (1 / TC - 1 / TH) * np.mean(-qa)


def main():
    Ka, Kb = K1_PAIR
    va, vb = K_value(Ka), K_value(Kb)
    K = vb + (vb - va) / ((Kb / Ka) ** P_EXP - 1)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump({"K": K}, fh)
    print("K(K1=%d) = %.12f, K(K1=%d) = %.12f, extrapolated K = %.12f" % (Ka, va, Kb, vb, K))


if __name__ == "__main__":
    main()
