"""Verifier for torus-entropy-gap.

The answer is Delta = Sdot_0 - sigma_x for the model in the instruction: the finite part of the reservoir entropy
flow in the small-mass limit minus the time-reversal relative-entropy rate of the observed x path of the limiting
diffusion. It is an exact rational number (tests/truth.json). solution/solve.py derives it in exact arithmetic;
Monte Carlo of the literal SDE and an exact stationarity check support it (authoring/evidence).
The answer must be an exact fraction; decimals are not accepted. The verifier never runs agent code.
"""
import json
import os
import re
from fractions import Fraction

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = Fraction(json.load(open(os.path.join(HERE, "truth.json")))["Delta"])


def _parse(v):
    if isinstance(v, bool) or isinstance(v, float):
        raise ValueError("Delta must be an exact fraction string 'p/q'")
    if isinstance(v, int):
        return Fraction(v)
    s = re.sub(r"\s+", "", str(v))
    if not re.fullmatch(r"[+-]?\d+(/[+-]?\d+)?", s):
        raise ValueError(f"not an exact fraction: {v!r}")
    return Fraction(s)


@pytest.fixture(scope="module")
def answer():
    with open(ANSWER) as fh:
        obj = json.load(fh)
    assert isinstance(obj, dict), "answer.json must hold a JSON object"
    keys = {re.sub(r"[^a-z0-9]", "", str(k).lower()): v for k, v in obj.items()}
    assert "delta" in keys, "missing field Delta"
    return _parse(keys["delta"])


def test_delta(answer):
    assert answer == TRUTH, f"Delta = {answer} is not the exact value"
