# Solver probes (blind agent attempts during authoring)

Each probe ran as an autonomous agent in its own empty directory, saw only the instruction, had the 9000 s budget, and was forbidden to read authoring files. "High" and "max" are reasoning-effort settings of the same frontier model. These records document calibration and serve as independent derivations of the truth; they are not part of the difficulty argument.

| version | target | high effort | max effort |
|---|---|---|---|
| user's original (smooth T) | Delta = 17/420 | solved (platform difficulty probe) | — |
| + uniform electric field, non-uniform field b(x) | Delta = 0.0134311455014 | solved, 19 and 26 min | solved, 12 and 23 min |
| + first finite-mass coefficient c | Delta and c = 0.00627578091037 | solved, 3 and 8 min | solved, 10 and 22 min |
| temperature step (shipped) | K = 0.0833929606 | solved, 24 min (K = 0.083393) | solved, 27 min (K = 0.08339296) |

Smooth versions are reachable by a generic small-mass (Hilbert) expansion carried to any order. The shipped step version needs the boundary-layer reduction and a half-space kinetic solver; both probes found it independently (same pressure-continuity argument, same energy-flux identity, Hermite moment method with Schur subspaces), and their values agree with `tests/truth.json` to 5e-7 (high) and 4e-9 (max). The max probe also confirmed the asymptotics by a full-SDE Monte Carlo (0.0834 +- 0.0007 at m = 0.0025).
