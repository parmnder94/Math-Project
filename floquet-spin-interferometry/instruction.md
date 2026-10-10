An isolated spin $ j=3/2 $ has $ \hbar=1 $, identity $ I_4 $, standard $ J_x,J_y,J_z $, and $ J_z|m\rangle=m|m\rangle $. Define

$$
P=|1/2\rangle\langle1/2|+|-1/2\rangle\langle-1/2|,\quad Q=I_4-P,\quad X=J_x^2-\frac{7P+3Q}{4}.
$$

Fix $ \Delta>0 $. For each positive integer $ N $, set

$$
T_N=\frac{2\pi N}{\Delta},\qquad
\epsilon_N=\left(10\,\Delta T_N\right)^{-1/8},
$$

$$
a_N=\epsilon_N^2+\frac43\epsilon_N^4+\frac{148}{45}\epsilon_N^6,\qquad
b_N=\epsilon_N^2+\frac53\epsilon_N^4+\frac{167}{36}\epsilon_N^6.
$$

Three pulses are indexed by $ r=1,2,3 $, with

$$
(\alpha_1,\alpha_2,\alpha_3)=\left(0,\frac{\pi}{3},\frac{2\pi}{3}\right),\qquad
(\phi_1,\phi_2,\phi_3)=\left(\frac{\pi}{6},\frac{\pi}{4},\frac{\pi}{3}\right).
$$

Each pulse uses its own $ 0\le t\le T_N $, $ s=t/T_N $, and

$$
J_r=J_x\cos\alpha_r+J_y\sin\alpha_r,\quad K_r=PJ_rP,\quad
R_r(s)=e^{-i\pi sJ_r},\quad D_r(s)=e^{-i\pi sK_r/2}.
$$

Writing $ \theta_r(t)=\Delta t+\phi_r $, define

$$
H_{r,N}(t)=R_r(s)\left[2\Delta Q+
\Delta D_r(s)\left(2\epsilon_N\cos\theta_r(t)J_x+a_N(P-Q)-2b_N\cos(2\theta_r(t))X\right)D_r^\dagger(s)\right]R_r^\dagger(s).
$$

Let $ U_{r,N} $ and $ \widetilde U_{r,N} $ be the final propagators of $ H_{r,N}(t) $ and $ H_{r,N}(T_N-t) $, with $ i\dot U=HU $, $ U(0)=I_4 $.

A control qubit and the spin start in $ |+\rangle_C\langle+|\otimes P/2 $, where $ |\pm\rangle_C=(|0\rangle_C\pm|1\rangle_C)/\sqrt2 $. Apply

$$
\mathcal V_N=|0\rangle_C\langle0|\otimes U_{3,N}U_{2,N}U_{1,N}
+|1\rangle_C\langle1|\otimes\widetilde U_{1,N}\widetilde U_{2,N}\widetilde U_{3,N}.
$$

With $ \Delta $ fixed, the probability $ p_N $ of measuring $ |-\rangle_C $ has an expansion $ p_N=p_\infty+c_1\epsilon_N+c_2\epsilon_N^2+c_3\epsilon_N^3+O(\epsilon_N^4) $ as $ N\to\infty $. Find $ p_\infty $, $ c_1 $, $ c_2 $ and $ c_3 $.

Write `/app/output/answer.json` containing a JSON object with four numeric fields: `p`, holding $ p_\infty $, and `c1`, `c2`, `c3`, holding $ c_1 $, $ c_2 $, $ c_3 $. Each is compared with the exact value and must agree to an absolute error of at most $ 10^{-9} $.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
