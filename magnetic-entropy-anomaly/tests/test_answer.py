"""Verifier for magnetic-entropy-anomaly.

The answer is K = lim_{m->0} sqrt(m) S_dot_m for the model in the instruction (a temperature step on a torus with a
rotating drag tensor and a magnetic field). It has no closed form: it is fixed by a kinetic boundary-layer problem
at the two temperature jumps. tests/truth.json holds it to about 1e-7 relative. It is computed by the reference
layer solver at high Hermite order and checked by finite-mass solutions and Monte Carlo (authoring/evidence).
The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = float(json.load(open(os.path.join(HERE, "truth.json")))["K"])
REL_TOL = 1e-4


def _parse(v):
    if isinstance(v, bool):
        raise ValueError("boolean is not a number")
    if isinstance(v, (int, float)):
        x = float(v)
    else:
        s = re.sub(r"\s+", "", str(v))
        if re.fullmatch(r"[+-]?\d+/[+-]?\d+", s):
            p, q = s.split("/")
            x = int(p) / int(q)
        else:
            x = float(s)
    if not math.isfinite(x):
        raise ValueError("non-finite value")
    return x


@pytest.fixture(scope="module")
def answer():
    with open(ANSWER) as fh:
        obj = json.load(fh)
    assert isinstance(obj, dict), "answer.json must hold a JSON object"
    keys = {re.sub(r"[^a-z0-9]", "", str(k).lower()): v for k, v in obj.items()}
    assert "k" in keys, "missing field K"
    return _parse(keys["k"])


def test_K(answer):
    rel = abs(answer - TRUTH) / abs(TRUTH)
    assert rel <= REL_TOL, f"K = {answer!r} has relative error {rel:.2e} (tolerance {REL_TOL:.0e})"
