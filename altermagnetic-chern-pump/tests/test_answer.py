"""Verifier for altermagnetic-chern-pump.

The answer is the second Chern number C_2 = deg(d/|d| : T^4 -> S^4), an integer. tests/truth.json holds it. It
follows from the exact north-pole zero count on the l-torus (D_l = -37) and the double covering x -> l (det L = 2);
solution/solve.py reproduces it in exact arithmetic, and preimage counts of eight other regular values agree
(authoring/evidence). The verifier only parses /app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = int(json.load(open(os.path.join(HERE, "truth.json")))["C2"])


def _parse_int(v):
    if isinstance(v, bool):
        raise ValueError("boolean is not an integer")
    if isinstance(v, int):
        return v
    if isinstance(v, float):
        if not math.isfinite(v) or v != int(v):
            raise ValueError(f"{v!r} is not an integer")
        return int(v)
    s = re.sub(r"\s+", "", str(v)).replace("−", "-")
    if not re.fullmatch(r"[+-]?\d+", s):
        raise ValueError(f"{v!r} is not an integer")
    return int(s)


@pytest.fixture(scope="module")
def answer():
    with open(ANSWER) as fh:
        obj = json.load(fh)
    assert isinstance(obj, dict), "answer.json must hold a JSON object"
    keys = {re.sub(r"[^a-z0-9]", "", str(k).lower()): v for k, v in obj.items()}
    assert "c2" in keys, "missing field C2"
    return _parse_int(keys["c2"])


def test_chern_number(answer):
    assert answer == TRUTH, f"C2 = {answer} is not the exact value"
