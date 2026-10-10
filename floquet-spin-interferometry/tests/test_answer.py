"""Verifier for floquet-spin-interferometry.

The answer is the exact N -> infinity probability of measuring |->_C for the interferometer in the instruction.
tests/truth.json holds it as a reduced fraction. It follows from the eighth-order Floquet effective Hamiltonian
(computed by an exact Feshbach recurrence), the slow-frame connection and exact finite matrix algebra; solution/solve.py
reproduces it in exact arithmetic, and direct simulations of the literal protocol approach it (authoring/evidence).
The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re
from fractions import Fraction

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = Fraction(json.load(open(os.path.join(HERE, "truth.json")))["p"])
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
    assert "p" in keys, "missing field p"
    return _parse(keys["p"])


def test_probability(answer):
    value, exact = answer
    if exact:
        assert value == TRUTH, f"p = {value} is not the exact value"
    else:
        assert abs(value - float(TRUTH)) <= FLOAT_TOL, f"p = {value} differs from the exact value"
