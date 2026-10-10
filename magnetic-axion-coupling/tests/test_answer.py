"""Verifier for magnetic-axion-coupling.

The answer is the axion angle theta (Chern-Simons magnetoelectric coupling) of the occupied bands of the instruction's
model, taken in (-pi, pi]. It has no closed form. tests/truth.json holds it to 12 digits; it was computed from the
Chern-Simons definition in a smooth projected gauge on stretched grids, and an independent second-Chern-form
integration along a gapped path from the time-reversal-symmetric insulator agrees (authoring/evidence).
The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = float(json.load(open(os.path.join(HERE, "truth.json")))["theta"])
ABS_TOL = 1e-3


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
    assert "theta" in keys, "missing field theta"
    return _parse(keys["theta"])


def test_axion_angle(answer):
    err = abs(answer - TRUTH)
    assert err <= ABS_TOL, f"theta = {answer!r} differs from the exact value by {err:.2e} (tolerance {ABS_TOL:.0e})"
