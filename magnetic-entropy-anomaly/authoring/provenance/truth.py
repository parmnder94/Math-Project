"""Ground truth for magnetic-entropy-anomaly: the reference layer solver (solution/solve.py) run at larger normal
Hermite orders K1 = 120 and 160 (K2 = 16, 8 angles), Richardson-extrapolated with the measured exponent 2.87.
Usage: python3 truth.py <solution_dir>   (about 20 min on 2 cores)"""
import importlib.util, json, os, sys
spec = importlib.util.spec_from_file_location("solve", os.path.join(sys.argv[1], "solve.py"))
S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
a, b = S.K_value(120), S.K_value(160)
K = b + (b - a) / ((160 / 120) ** 2.87 - 1)
print(json.dumps({"K_K1_120": a, "K_K1_160": b, "K": K}))
