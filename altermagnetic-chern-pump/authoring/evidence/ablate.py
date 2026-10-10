"""Values produced by plausible shortcuts (exact arithmetic). Usage: python3 ablate.py <out.json>"""
import itertools
import json
import sys
from fractions import Fraction as F

W = {0: F(1), 1: F(-1, 2), 2: F(-1, 2)}   # Re omega^r


def A(a, b):
    return (a[0] * a[1] + a[2] * a[3] + sum(a[i] * b[i] for i in range(4))) % 3


def B(a, b):
    return (a[0] * a[2] + a[1] * a[3] + a[0] * b[1] + a[1] * b[2] + a[2] * b[3] + a[3] * b[0]) % 3


def Q2(b):
    return b[0] * b[1] + b[1] * b[2] + b[2] * b[3] + b[3] * b[0]


def degree(bracket=lambda x, y: F(1, 5) + W[x] + F(3, 5) * W[y], use_h=True, mode="north"):
    """Degree on the l-torus from the zero data. mode='north' counts d5>0 zeros; 'signed' sums c*sgn(d5)."""
    D = 0
    for n in itertools.product(range(6), repeat=4):
        a = [k % 3 for k in n]; b = [k % 2 for k in n]
        v = bracket(A(a, b), B(a, b)) * ((-1) ** Q2(b) if use_h else 1)
        assert v != 0
        c = (-1) ** sum(b)
        D += c * (1 if v > 0 else 0) if mode == "north" else c * (1 if v > 0 else -1)
    return D


def main():
    detL = 2
    ref = detL * degree()
    # independence (factorization) assumption: replace Z_AB by Z_A Z_B / 81 in the character-sum formula
    W0 = WA = WB = WABf = F(0)
    for b in itertools.product(range(2), repeat=4):
        w = (-1) ** (sum(b) + Q2(b))
        ZA = sum(1 for a in itertools.product(range(3), repeat=4) if A(a, b) == 0)
        ZB = sum(1 for a in itertools.product(range(3), repeat=4) if B(a, b) == 0)
        W0 += w; WA += w * ZA; WB += w * ZB; WABf += w * F(ZA * ZB, 81)
    fact = detL * (-81 * W0 + 2 * WA + 2 * WB - 2 * WABf) / 2
    out = {
        "reference C2": ref,
        "covering factor det L = 2 omitted (degree on the l-torus)": degree(),
        "signed north-minus-south sum not halved, times det L": detL * degree(mode="signed"),
        "pairing prefactor prod chi dropped": detL * degree(use_h=False),
        "Z_2 term dropped from d_5": detL * degree(bracket=lambda x, y: F(1, 5) + W[x]),
        "Z_1 term dropped from d_5": detL * degree(bracket=lambda x, y: F(1, 5) + F(3, 5) * W[y]),
        "joint count factorized, Z_AB = Z_A Z_B / 81": str(fact),
        "S^4 orientation reversed": -ref,
    }
    out = {k: (int(v) if isinstance(v, (int, F)) and F(v).denominator == 1 else v) for k, v in out.items()}
    json.dump(out, open(sys.argv[1], "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
