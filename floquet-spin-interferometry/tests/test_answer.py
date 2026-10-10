"""Verifier for floquet-spin-interferometry.

The answer is (p_inf, c1, c2, c3) in p_N = p_inf + c1 eps_N + c2 eps_N^2 + c3 eps_N^3 + O(eps_N^4) for the
interferometer in the instruction, real numbers with no closed form. tests/truth.json holds p_inf and c1 to 30 digits
and c2, c3 to 18. p_inf follows from the eighth-order Floquet effective Hamiltonian (exact Feshbach recurrence) and the
slow-frame connection; c1 from the first-order micromotion kicks at the pulse boundaries; c2 and c3 from a smooth
interpolant F(eps) with p_N = F(eps_N) + O(eps_N^8) built from the exact one-period Floquet propagator and the
micromotion-dressed slow connection. solution/solve.py reproduces all four, and high-precision simulations of the
literal protocol confirm them (authoring/evidence).
The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = {k: float(v) for k, v in json.load(open(os.path.join(HERE, "truth.json"))).items()}
FIELDS = ["p", "c1", "c2", "c3"]
ABS_TOL = 1e-9


def _parse(v):
    if isinstance(v, bool):
        raise ValueError("boolean is not a number")
    x = float(v) if isinstance(v, (int, float)) else float(re.sub(r"\s+", "", str(v)))
    if not math.isfinite(x):
        raise ValueError("non-finite value")
    return x


@pytest.fixture(scope="module")
def answer():
    with open(ANSWER) as fh:
        obj = json.load(fh)
    assert isinstance(obj, dict), "answer.json must hold a JSON object"
    keys = {re.sub(r"[^a-z0-9]", "", str(k).lower()): v for k, v in obj.items()}
    return keys


@pytest.mark.parametrize("field", FIELDS)
def test_field(answer, field):
    assert field in answer, f"missing field {field}"
    value = _parse(answer[field])
    err = abs(value - TRUTH[field])
    assert err <= ABS_TOL, f"{field} = {value!r} differs from the exact value by {err:.2e} (tolerance {ABS_TOL:.0e})"
