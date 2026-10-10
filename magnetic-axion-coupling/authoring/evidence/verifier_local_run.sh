#!/bin/bash
# usage: verifier_local_run.sh <dir containing answer.json or empty dir>
# Builds the verifier image from ../../tests and runs test.sh with the given directory mounted at /app/output.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
docker image inspect mac-tests >/dev/null 2>&1 || docker build -q -t mac-tests "$HERE/../../tests" >/dev/null
docker run --rm --network none -v "$1":/app/output:ro mac-tests bash -c \
  "bash /tests/test.sh >/tmp/log 2>&1; tail -1 /tmp/log; echo reward=\$(cat /logs/verifier/reward.txt)"
