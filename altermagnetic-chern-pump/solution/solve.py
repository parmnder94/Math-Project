"""Reference solution for altermagnetic-chern-pump (exact integer arithmetic, standard library only).

C_2 is the degree of d/|d| : T^4 -> S^4. The derivation (author's golden solution):

1. d_1 = ... = d_4 = 0 forces l_i = pi n_i / 3 with n_i in Z_6. Write b_i = n_i mod 2, a_i = n_i mod 3. There
   p_i = (-1)^{b_i}, u_i = omega^{a_i}, S(p_i, u_i) = u_i^{b_i}, Xi(omega^a, omega^b) = omega^{ab}, so
   Z_1 = omega^A, Z_2 = omega^B with
       A = a1 a2 + a3 a4 + sum_i a_i b_i,   B = a1 a3 + a2 a4 + a1 b2 + a2 b3 + a3 b4 + a4 b1   (mod 3),
   and prod chi = (-1)^{Q_2(b)}, Q_2 = b1 b2 + b2 b3 + b3 b4 + b4 b1. The bracket takes the values 9/5, 9/10, 3/10,
   -3/5, never 0, so d never vanishes (the BdG spectrum is gapped).
2. North-pole regular value: only zeros with d_5 > 0 count, each with chirality c(b) = (-1)^{sum b_i}. Since
   sum_b c(b) = 0, 2 D_l = -81 W_0 + 2 W_A + 2 W_B - 2 W_AB with w(b) = c(b) (-1)^{Q_2(b)} and
   W_0 = sum w, W_A = sum w Z_A, W_B = sum w Z_B, W_AB = sum w Z_AB (zero counts of A, B, A=B=0 over F_3^4).
3. Completing squares over F_3: Z_A = 24 + 9 [rho_A = 0], rho_A = b1 b2 + b3 b4; Z_B = 24 + 9 [rho_B = 0],
   rho_B = b1 b3 + b2 b4. W_AB uses the characters omega^{rA + sB} (Gauss sums), it does not factorize.
4. D_l = -37 on the l-torus, and x -> l has det L = 2 (orientation-preserving double cover): C_2 = 2 D_l = -74.

This script evaluates the character-sum route of steps 2-4 and, as an internal cross-check, the exact
north-pole zero count straight from the definitions in the instruction, in exact arithmetic in Q(omega).
"""
import itertools
import json
import os
from fractions import Fraction as F

OUT = "/app/output"


class Qw:
    """x + y*omega in Q(omega), omega^2 = -1 - omega."""

    def __init__(self, x, y=0):
        self.x, self.y = F(x), F(y)

    def __add__(self, o):
        o = o if isinstance(o, Qw) else Qw(o)
        return Qw(self.x + o.x, self.y + o.y)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, Qw) else Qw(o)
        return Qw(self.x - o.x, self.y - o.y)

    def __rsub__(self, o):
        return Qw(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Qw) else Qw(o)
        xx, xy, yy = self.x * o.x, self.x * o.y + self.y * o.x, self.y * o.y
        return Qw(xx - yy, xy - yy)

    __rmul__ = __mul__

    def __truediv__(self, k):
        return Qw(self.x / k, self.y / k)

    def re(self):
        return self.x - self.y / 2


OMEGA = Qw(0, 1)


def wpow(k):
    r = Qw(1)
    for _ in range(k % 3):
        r = r * OMEGA
    return r


def chi(p, q):
    return (1 + p + q - p * q) / 2


def S(p, u):
    return (1 + p + (1 - p) * u) / 2


def Xi(u, v):
    w, w2 = OMEGA, OMEGA * OMEGA
    return ((1 + v + v * v) + u * (1 + w2 * v + w * v * v) + u * u * (1 + w * v + w2 * v * v)) / 3


def d5_at_zero(n):
    """d_5 at l_i = pi n_i / 3, evaluated from the instruction's definitions: p_i = (-1)^{n_i}, u_i = omega^{n_i}."""
    p = [F((-1) ** k) for k in n]
    u = [wpow(k) for k in n]
    pc = p + [p[0]]
    Z1 = Xi(u[0], u[1]) * Xi(u[2], u[3])
    Z2 = Xi(u[0], u[2]) * Xi(u[1], u[3])
    for i in range(4):
        Z1 = Z1 * S(p[i], u[i])
        Z2 = Z2 * S(pc[i + 1], u[i])
    pref = F(1)
    for i in range(4):
        pref *= chi(pc[i], pc[i + 1])
    return pref * (F(1, 5) + Z1.re() + F(3, 5) * Z2.re())


def det_int(M):
    M = [[F(v) for v in row] for row in M]
    n, det = len(M), F(1)
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] != 0), None)
        if piv is None:
            return F(0)
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [M[r][k] - f * M[c][k] for k in range(n)]
    return det


def character_sum_route():
    def A(a, b):
        return (a[0] * a[1] + a[2] * a[3] + sum(a[i] * b[i] for i in range(4))) % 3

    def B(a, b):
        return (a[0] * a[2] + a[1] * a[3] + a[0] * b[1] + a[1] * b[2] + a[2] * b[3] + a[3] * b[0]) % 3

    def w(b):
        return (-1) ** (sum(b) + b[0] * b[1] + b[1] * b[2] + b[2] * b[3] + b[3] * b[0])

    W0 = WA = WB = WAB = 0
    for b in itertools.product(range(2), repeat=4):
        rhoA, rhoB = b[0] * b[1] + b[2] * b[3], b[0] * b[2] + b[1] * b[3]
        ZA = 24 + 9 * (rhoA % 3 == 0)   # rank-4 hyperbolic form level-set counts over F_3
        ZB = 24 + 9 * (rhoB % 3 == 0)
        # joint count by character orthogonality: (1/9) sum_{r,s} sum_a omega^{rA+sB}
        tot = Qw(0)
        for r in range(3):
            for s in range(3):
                for a in itertools.product(range(3), repeat=4):
                    tot = tot + wpow(r * A(a, b) + s * B(a, b))
        ZAB = tot / 9
        assert ZAB.y == 0 and ZAB.x.denominator == 1
        W0 += w(b); WA += w(b) * ZA; WB += w(b) * ZB; WAB += w(b) * int(ZAB.x)
    two_D = -81 * W0 + 2 * WA + 2 * WB - 2 * WAB
    return {"W0": W0, "WA": WA, "WB": WB, "WAB": WAB, "D_l": F(two_D, 2)}


def direct_zero_count():
    D, gap = 0, None
    for n in itertools.product(range(6), repeat=4):
        v = d5_at_zero(n)
        gap = abs(v) if gap is None else min(gap, abs(v))
        if v > 0:
            D += (-1) ** sum(n)   # sign of the Jacobian 3^4 prod cos(3 l_i)
    return D, gap


def main():
    L = [[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [-1, -1, -1, 1]]
    detL = det_int(L)
    cs = character_sum_route()
    D_direct, gap = direct_zero_count()
    assert gap > 0, "d vanishes somewhere"
    assert cs["D_l"] == D_direct, (cs, D_direct)
    C2 = int(detL * cs["D_l"])
    print(f"W0={cs['W0']} WA={cs['WA']} WB={cs['WB']} WAB={cs['WAB']}  D_l={cs['D_l']}  "
          f"direct D_l={D_direct}  min|d5| at zeros={gap}  det L={detL}  C2={C2}")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "answer.json"), "w") as fh:
        json.dump({"C2": C2}, fh)


if __name__ == "__main__":
    main()
