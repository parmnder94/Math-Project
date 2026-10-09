"""Ground truth for heat-trace-corrections.
Usage: python3 truth.py <out_dir>
Runs solution/solve.py with growing cutoffs (class words up to LMAX, pairs of classes up to total length
LPAIR) and writes <out_dir>/truth.json (LMAX = 13, LPAIR = 14) and <out_dir>/truth_report.json."""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
SOLVE = os.path.join(HERE, "..", "..", "solution", "solve.py")
os.makedirs(sys.argv[1], exist_ok=True)
report = {}
for lmax, lpair in [(10, 10), (11, 12), (12, 13), (13, 14)]:
    out = os.path.join(sys.argv[1], f"run_{lmax}_{lpair}.json"); t0 = time.time()
    subprocess.run([sys.executable, SOLVE], check=True, stdout=subprocess.DEVNULL,
                   env=dict(os.environ, LMAX=str(lmax), LPAIR=str(lpair), ANSWER_PATH=out))
    r = json.load(open(out)); r["seconds"] = round(time.time() - t0); report[f"lmax{lmax}_lpair{lpair}"] = r
    print(lmax, lpair, r, flush=True)
best = report["lmax13_lpair14"]
json.dump({k: best[k] for k in ("z_inf", "m", "m1", "v", "v1")}, open(os.path.join(sys.argv[1], "truth.json"), "w"), indent=2)
json.dump(report, open(os.path.join(sys.argv[1], "truth_report.json"), "w"), indent=2)
