#!/bin/bash
set -euo pipefail
# Match BLAS threads to the 2-CPU task limit (large grading hosts otherwise oversubscribe).
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p /app/output
cd /app
python3 /solution/solve.py
test -s /app/output/answer.json
