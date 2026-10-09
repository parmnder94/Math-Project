Use dimensionless units with $k_B=1$. A charged Brownian particle moves on the flat torus $(x,y)\in[0,2\pi)^2$. The surrounding medium is held at two temperatures: $T=2$ on the stripe $0<x<\pi$ and $T=1$ on the stripe $\pi<x<2\pi$, so the temperature jumps across the lines $x=0$ and $x=\pi$. The stationary underdamped dynamics is

$$
d\mathbf X_t=\mathbf V_t\,dt,\qquad
m\,d\mathbf V_t=\big[-\Gamma(y_t)\mathbf V_t+J\mathbf V_t\big]dt+\sqrt{2T(\mathbf X_t)\Gamma(y_t)}\,d\mathbf W_t,
$$

where $m>0$, $\mathbf X_t=(x_t,y_t)$, $\mathbf W_t$ is standard two-dimensional Brownian motion, the square root is the symmetric positive root, and

$$
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
\Gamma(y)=R(y)\begin{pmatrix}1&0\\0&2\end{pmatrix}R(y)^{\mathsf T},\qquad
R(y)=\begin{pmatrix}\cos y&-\sin y\\\sin y&\cos y\end{pmatrix}.
$$

The term $J\mathbf V$ is the Lorentz force from a uniform perpendicular magnetic field; there is no other force. At each $m$, define the stationary entropy flow into the thermal environment by

$$
\dot S_m=\lim_{\tau\to\infty}\frac1\tau\mathbb E\!\left[
\int_0^\tau\frac{\mathbf V_t^{\mathsf T}\Gamma\mathbf V_t}{T}\,dt
-\int_0^\tau\frac{\mathbf V_t^{\mathsf T}\sqrt{2T\Gamma}}{T}\circ d\mathbf W_t
\right],
$$

where coefficients are evaluated at $\mathbf X_t$ and $\circ$ denotes Stratonovich integration.

As $m\to0$ this entropy flow diverges, but the limit $K=\lim_{m\to0}\sqrt m\,\dot S_m$ exists and is finite. Find $K$.

Write `/app/output/answer.json` containing a JSON object with one numeric field, `K`, holding $K$. It is compared with the exact value and must agree to a relative error of at most $10^{-4}$.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
