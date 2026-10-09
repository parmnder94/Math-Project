"""Verifier for heat-trace-corrections.

The five constants are exact mathematical quantities of the random-permutation model in the instruction.
tests/truth.json holds them to about 1e-8 relative or better (authoring/provenance/truth.py: exact
class sums of the heat kernel, Linial-Puder quotient expansion of word-map moments to order 1/N, classes
up to length 13 and pairs of classes up to total length 14; cross-checked against exact finite-N averages
over S_N x S_N and Richardson extrapolation of the unpruned expansion). The verifier only parses
/app/output/answer.json; it never runs agent code.
"""
import json
import math
import os
import re

import pytest

ANSWER = "/app/output/answer.json"
HERE = os.path.dirname(os.path.abspath(__file__))
TRUTH = json.load(open(os.path.join(HERE, "truth.json")))
METRICS = "/logs/verifier/metrics.json"
REL_TOL = 1e-6
FIELDS = ["z_inf", "m", "m1", "v", "v1"]

_metrics = {}


def _dump():
    os.makedirs(os.path.dirname(METRICS), exist_ok=True)
    with open(METRICS, "w") as fh:
        json.dump(_metrics, fh, indent=2)


def _norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).strip().lower())


def _num(v):
    if isinstance(v, bool):
        raise ValueError("boolean is not a number")
    x = float(str(v).strip())
    if not math.isfinite(x):
        raise ValueError("non-finite value")
    return x


@pytest.fixture(scope="module")
def answer():
    with open(ANSWER) as fh:
        obj = json.load(fh)
    assert isinstance(obj, dict), "answer.json must hold a JSON object"
    keys = {_norm(k): v for k, v in obj.items()}
    out = {}
    for f in FIELDS:
        assert _norm(f) in keys, f"missing field {f}"
        out[f] = _num(keys[_norm(f)])
    return out


def test_answer_parses(answer):
    assert set(answer) == set(FIELDS)


@pytest.mark.parametrize("field", FIELDS)
def test_value(answer, field):
    got, want = answer[field], TRUTH[field]
    rel = abs(got - want) / abs(want)
    _metrics[field] = {"reported": got, "relative_error": rel}
    _dump()
    assert rel <= REL_TOL, f"{field}: relative error {rel:.2e} exceeds {REL_TOL:.0e}"
