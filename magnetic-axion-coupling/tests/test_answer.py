"""Verifier for magnetic-axion-coupling.

The answer is the orbital magnetoelectric coefficient alpha_zz = dP_z/dB_z of the instruction's model in units of
e^2/h, taken in (-1/2, 1/2]. It has no closed form. tests/truth.json holds it to 10 digits; it was computed from the
converse response dM_z/dE_z of slabs, and an independent Hofstadter-supercell computation of dP_z/dB_z agrees
(authoring/evidence).
The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = float(json.load(open(os.path.join(HERE, "truth.json")))["alpha_zz"])
ABS_TOL = 3e-3


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
    assert "alphazz" in keys, "missing field alpha_zz"
    return _parse(keys["alphazz"])


def test_magnetoelectric_coefficient(answer):
    err = abs(answer - TRUTH)
    assert err <= ABS_TOL, f"alpha_zz = {answer!r} differs from the exact value by {err:.2e} (tolerance {ABS_TOL:.0e})"
