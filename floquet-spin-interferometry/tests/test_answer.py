"""Verifier for floquet-spin-interferometry.

The answer is the N -> infinity probability of measuring |->_C for the interferometer in the instruction, a real
number with no closed form. tests/truth.json holds it to 30 digits. It follows from the eighth-order Floquet effective
Hamiltonian (exact Feshbach recurrence), the slow-frame connection and 40-digit matrix exponentials; solution/solve.py
reproduces it, and direct simulations of the literal protocol approach it (authoring/evidence).
The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = float(json.load(open(os.path.join(HERE, "truth.json")))["p"])
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
    assert "p" in keys, "missing field p"
    return _parse(keys["p"])


def test_probability(answer):
    err = abs(answer - TRUTH)
    assert err <= ABS_TOL, f"p = {answer!r} differs from the exact value by {err:.2e} (tolerance {ABS_TOL:.0e})"
