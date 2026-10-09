Use dimensionless units with $k_B=1$. A charged Brownian particle moves on the flat torus $(x,y)\in[0,2\pi)^2$. Its stationary underdamped dynamics is

$$
d\mathbf X_t=\mathbf V_t\,dt,\qquad
m\,d\mathbf V_t=[\nabla T(\mathbf X_t)-\Gamma(y_t)\mathbf V_t+J\mathbf V_t]dt+\sqrt{2T(\mathbf X_t)\Gamma(y_t)}\,d\mathbf W_t,
$$

where $m>0$, $\mathbf W_t$ is standard two-dimensional Brownian motion, the square root is the symmetric positive root, and

$$
T(x,y)=1+\frac35\cos x,\qquad
J=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
\Gamma(y)=R(y)\begin{pmatrix}1&0\\0&2\end{pmatrix}R(y)^{\mathsf T},\qquad
R(y)=\begin{pmatrix}\cos y&-\sin y\\\sin y&\cos y\end{pmatrix}.
$$

The term $J\mathbf V$ is the Lorentz force from a fixed perpendicular magnetic field; the external force is exactly $\nabla T$. At each $m$, define the stationary entropy flow into the thermal environment by

$$
\dot S_m=\lim_{\tau\to\infty}\frac1\tau\mathbb E\!\left[
\int_0^\tau\frac{\mathbf V_t^{\mathsf T}\Gamma\mathbf V_t}{T}\,dt
-\int_0^\tau\frac{\mathbf V_t^{\mathsf T}\sqrt{2T\Gamma}}{T}\circ d\mathbf W_t
\right],
$$

where coefficients are evaluated at $\mathbf X_t$ and $\circ$ denotes Stratonovich integration.

Let $\mathbb P_\tau^{(0)}$ be the stationary position-path law on $[0,\tau]$ obtained by taking $m\to0$ first. Let $\mathbb P_\tau^{(0),R}$ be its image under $\mathbf X_t\mapsto\mathbf X_{\tau-t}$; the magnetic field is held fixed in this comparison. Define

$$
\mathcal I_0=\lim_{\tau\to\infty}\frac1\tau
D_{\mathrm{KL}}(\mathbb P_\tau^{(0)}\Vert\mathbb P_\tau^{(0),R}),
\qquad
D_{\mathrm{KL}}(P\Vert Q)=\mathbb E_P\ln\frac{dP}{dQ}.
$$

Find the exact value of $\Delta=\lim_{m\to0}\dot S_m-\mathcal I_0$ as a reduced rational number.

Write `/app/output/answer.json` containing a JSON object with one field, `delta`, holding $\Delta$ as a string of the form `"p/q"`. A decimal number is also accepted if it is within $10^{-12}$ of the exact value.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
