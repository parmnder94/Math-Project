"""Reference solution for magnetic-entropy-anomaly: exact small-mass limit of the heat entropy rate minus the
position-path relative entropy rate (fixed magnetic field), computed in exact rational arithmetic with SymPy.

Model: dX = V dt, m dV = [grad T - Gamma(y) V + b J V] dt + sqrt(2 T Gamma) dW on the torus, with
T = 1 + a cos x, Gamma(y) = R(y) diag(g1, g2) R(y)^T.  Shipped parameters: a = 3/5, (g1, g2) = (1, 2), b = 1.

1. Small-mass limit (u = sqrt(m) V): fast Ornstein-Uhlenbeck generator L0 = -(Gamma - bJ)u.grad_u + T Gamma:hess_u
   with invariant Maxwellian of covariance T I; limiting position diffusion with mobility M = (Gamma - bJ)^-1,
   drift M grad T + T div M, diffusion D = T M_s; stationary density uniform; current j = rho0 M_a grad T.
2. I0 = int j^T D^-1 j / rho0 (fixed-field time reversal of the limiting diffusion).
3. Heat rate: S_dot_m = -<S>/(2 sqrt m) with the centred cubic S = (g.u)(|u|^2 - 4T)/T^2; solving L0 chi = S in the
   cubic Hermite sector of the local drag frame gives lim S_dot_m = -(1/4) <chi S>_0.
"""
import json, os
from fractions import Fraction
import sympy as sp

A_AMP, G1, G2, BFIELD = sp.Rational(3, 5), sp.Integer(1), sp.Integer(2), sp.Integer(1)
OUT = os.environ.get("ANSWER_PATH", "/app/output/answer.json")


def delta(a=A_AMP, g1=G1, g2=G2, b=BFIELD, field_in_reversal=True):
    x, y, x1, x2 = sp.symbols('x y xi1 xi2', real=True)
    J = sp.Matrix([[0, -1], [1, 0]])
    R = sp.Matrix([[sp.cos(y), -sp.sin(y)], [sp.sin(y), sp.cos(y)]])
    T = 1 + a * sp.cos(x); Tp = sp.diff(T, x)
    IT = sp.nsimplify(1 - sp.sqrt(1 - a ** 2))               # <T'^2 / T>_x, exact for T = 1 + a cos x
    # ---- heat rate: cubic Hermite sector in the local frame xi = R^T u / sqrt(T) (all coefficients rational)
    G0 = sp.diag(g1, g2) - b * J
    def L0(f):
        grad = sp.Matrix([sp.diff(f, x1), sp.diff(f, x2)]); xi = sp.Matrix([x1, x2])
        return sp.expand(-(G0 * xi).dot(grad) + g1 * sp.diff(f, x1, 2) + g2 * sp.diff(f, x2, 2))
    H = [x1**3 - 3*x1, (x1**2 - 1)*x2, x1*(x2**2 - 1), x2**3 - 3*x2]
    norms = sp.diag(6, 2, 2, 6)
    mons = [x1**3, x1**2*x2, x1*x2**2, x2**3, x1, x2]
    coef = lambda f: sp.Matrix([sp.Poly(f, x1, x2).coeff_monomial(mm) for mm in mons])
    Bm = sp.Matrix.hstack(*[coef(h) for h in H])
    A = sp.Matrix.hstack(*[(Bm.T * Bm).inv() * Bm.T * coef(-L0(h)) for h in H])   # -L0 H_j = sum_i A_ij H_i
    assert all(sp.simplify(e) == 0 for e in (Bm * A - sp.Matrix.hstack(*[coef(-L0(h)) for h in H])))
    p = [sp.Matrix([1, 0, 1, 0]), sp.Matrix([0, 1, 0, 1])]
    C = sp.Matrix(2, 2, lambda i, j: (p[i].T * norms * A.inv() * p[j])[0])          # <p_i A^-1 p_j>
    # S = T^-1/2 h.p, chi = -T^-1/2 sum_j h_j A^-1 p_j, h = R^T grad T  ->  lim S_dot = (1/4) <h^T C h / T>
    h = R.T * sp.Matrix([Tp, 0])
    yavg = lambda f: sp.integrate(sp.expand(sp.simplify(f)), (y, 0, 2 * sp.pi)) / (2 * sp.pi)
    Sdot = sp.nsimplify(yavg((h.T * C * h)[0] / Tp ** 2) / 4) * IT
    # ---- I0: current of the limiting diffusion and its relative entropy rate under fixed-field reversal
    M = sp.simplify((R * sp.diag(g1, g2) * R.T - b * J).inv())
    Ms, Ma = (M + M.T) / 2, (M - M.T) / 2
    jv = Ma * sp.Matrix([Tp, 0])                   # j / rho0 = div(T M_a) = M_a grad T  (M_a is constant)
    I0 = sp.nsimplify(yavg((jv.T * Ms.inv() * jv)[0] / Tp ** 2)) * IT if field_in_reversal else sp.Integer(0)
    return sp.nsimplify(Sdot), sp.nsimplify(I0), sp.nsimplify(Sdot - I0)


def main():
    sdot, i0, d = delta()
    ans = {"delta": str(Fraction(int(sp.numer(d)), int(sp.denom(d))))}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh: json.dump(ans, fh)
    print(f"lim S_dot = {sdot}, I0 = {i0}, Delta = {ans['delta']}")


if __name__ == "__main__":
    main()
