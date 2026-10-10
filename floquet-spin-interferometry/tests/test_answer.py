"""Verifier for floquet-spin-interferometry.

The answer is the pair (p_inf, c) in p_N = p_inf + c eps_N + O(eps_N^2) for the interferometer in the instruction,
real numbers with no closed form. tests/truth.json holds both to 30 digits. p_inf follows from the eighth-order Floquet
effective Hamiltonian (exact Feshbach recurrence) and the slow-frame connection; c from the first-order micromotion
kicks at the pulse boundaries. solution/solve.py reproduces both, and high-precision simulations of the literal
protocol confirm them (authoring/evidence).
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
FIELDS = ["p", "c"]
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
