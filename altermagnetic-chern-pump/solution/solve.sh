#!/bin/bash
set -euo pipefail
mkdir -p /app/output
cd /app
python3 /solution/solve.py
test -s /app/output/answer.json
