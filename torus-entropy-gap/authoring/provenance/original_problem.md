# Author-supplied original problem (provenance)

The task is a hardened version of the following problem supplied by the task author (Dr. Prince). The original
statement defined T = 1 + (3/5)cos x, theta = x + y + (2/3)cos x, Gamma = R diag(1,2) R^T, the fast covariance
matrix C = R [[11,-2],[-2,17]] R^T/9, the current h = (2/(5T), (2/3)T'), the force F = (1/T) div(T^2 C) + (Gamma - J)h,
reservoirs a = 1, 2 with friction a and temperature aT, the entropy flow with the subtraction 1/(9m) written out,
and the observed x record "holding every coefficient, including J, fixed" under reversal; it cited Hottovy,
McDaniel, Volpe and Wehr, Commun. Math. Phys. 336 (2015).

Author's derivation (36 steps, summarised): fast covariance Sigma = T C from B0 C0 + C0 B0^T = 2 Q0; Ito form
Sdot_m = (1/m) E[|w|^2/T - 3]; limiting drift and diffusion D = T sym(M C); stationary density T/(4 pi^2) and current
rho[h + (10/9) J grad T]; reconstruction of theta mod pi from the quadratic variation of x and its covariation with
a(theta); lift symmetry y -> y + pi; observed KL rate = full-position rate = 41249/28800; finite heat from a quadratic
Poisson observable and explicit cubic-moment tensors in the reservoir frame (four 4x4 linear systems, rotation
averages, density-weighted integrals) = 8657741/18370800; answer -141231719/146966400.

The reference pipeline of this task (solution/solve.py), which uses a different route (the homogenized generator
-<A L^-1 A> and Sdot_0 = -int f0 A L^-1 A Phi with exact averages), reproduces all three numbers exactly with the
original parameters (authoring/evidence/ablations.json, last entry).
