"""Verifier for magnetic-entropy-anomaly.

The answer is the exact rational Delta = lim_{m->0} S_dot_m - I_0 for the model in the instruction.
tests/truth.json holds it as a fraction. It was derived analytically and is recomputed in exact arithmetic
by solution/solve.py; a Monte Carlo of the underdamped dynamics checks the heat-rate limit independently
(authoring/evidence). The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re
from fractions import Fraction

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = Fraction(json.load(open(os.path.join(HERE, "truth.json")))["delta"])
FLOAT_TOL = 1e-12


def _parse(v):
    if isinstance(v, bool):
        raise ValueError("boolean is not a number")
    if isinstance(v, int):
        return Fraction(v), True
    if isinstance(v, float):
        if not math.isfinite(v):
            raise ValueError("non-finite value")
        return v, False
    s = re.sub(r"\s+", "", str(v))
    if re.fullmatch(r"[+-]?\d+(/[+-]?\d+)?", s):
        return Fraction(s), True
    x = float(s)
    if not math.isfinite(x):
        raise ValueError("non-finite value")
    return x, False


@pytest.fixture(scope="module")
def answer():
    with open(ANSWER) as fh:
        obj = json.load(fh)
    assert isinstance(obj, dict), "answer.json must hold a JSON object"
    keys = {re.sub(r"[^a-z0-9]", "", str(k).lower()): v for k, v in obj.items()}
    assert "delta" in keys, "missing field delta"
    return _parse(keys["delta"])


def test_delta(answer):
    value, exact = answer
    if exact:
        assert value == TRUTH, f"delta = {value} is not the exact value"
    else:
        assert abs(value - float(TRUTH)) <= FLOAT_TOL, f"delta = {value} differs from the exact value"
