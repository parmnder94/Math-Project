#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p /app/output
cd /app
pip install --no-cache-dir --quiet --root-user-action=ignore --disable-pip-version-check python-flint==0.9.0
python3 /solution/solve.py
test -s /app/output/answer.json
