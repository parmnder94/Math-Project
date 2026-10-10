import sys, sim3, sim4
sim3.mp.mp.dps = 64
if __name__ == "__main__":
    for N in [int(float(v)) for v in sys.argv[1:]]:
        p, eps = sim4.prob(N); print(f"N={N:.0e} eps={eps:.8f} p={p:.15f}", flush=True)
