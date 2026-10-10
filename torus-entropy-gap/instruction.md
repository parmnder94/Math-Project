Set $ k_B=1 $. A particle of mass $ m>0 $ moves on the torus $ \mathbf X=(x,y)\in[0,2\pi)^2 $ with velocity $ \mathbf V $. Let

$$
T=1+\frac45\cos x,\qquad \theta=x+y+\frac12\cos x,\qquad
\mathbf r_1=\begin{pmatrix}\cos\theta\\ \sin\theta\end{pmatrix},\qquad
\mathbf r_2=\begin{pmatrix}-\sin\theta\\ \cos\theta\end{pmatrix},\qquad
J=\begin{pmatrix}0&-1\\ 1&0\end{pmatrix},
$$

and let the force $ \mathbf F=(F_1,F_2) $ be

$$
F_1=\phi_1+\psi\cos2\theta+\chi\sin2\theta,\qquad F_2=\phi_2+\psi\sin2\theta-\chi\cos2\theta,
$$

$$
\phi_1=-\frac{237}{65}\sin x+\frac{2}{3T^2}-\frac{\sin x}{5T},\qquad
\phi_2=-\frac{1}{6T^2}-\frac{4\sin x}{5T},
$$

$$
\psi=\frac{72}{65}\sin x-\frac{9}{13}T-\frac{3}{26}T\sin x-\frac{1}{3T^2},\qquad
\chi=-\frac{18}{65}\sin x+\frac{15}{13}T-\frac{6}{13}T\sin x+\frac{2\sin x}{5T}.
$$

The particle is coupled to two heat reservoirs. Reservoir $ a\in\{1,2\} $ acts along $ \mathbf r_a $ with friction coefficient $ \gamma_a $ and temperature $ T_a $, where $ (\gamma_1,\gamma_2)=(1,3) $ and $ (T_1,T_2)=(T,2T) $. All spatial quantities are evaluated at $ \mathbf X_t $, and with $ u_a=\mathbf r_a\cdot\mathbf V_t $ the dynamics is

$$
d\mathbf X_t=\mathbf V_t\,dt,\qquad
m\,d\mathbf V_t=\Big[\mathbf F-\sum_{a=1}^2\gamma_au_a\mathbf r_a+\frac12J\mathbf V_t\Big]dt+\sum_{a=1}^2\sqrt{2\gamma_aT_a}\,\mathbf r_a\,dW_{a,t},
$$

where $ W_1,W_2 $ are independent standard Brownian motions. In the stationary state at mass $ m $, the entropy flow into the reservoirs is

$$
\dot S_m=\lim_{\tau\to\infty}\frac1\tau\,\mathbb E\sum_{a=1}^2\left[\int_0^\tau\frac{\gamma_au_a^2}{T_a}\,dt-\int_0^\tau\sqrt{\frac{2\gamma_a}{T_a}}\,u_a\circ dW_{a,t}\right],
$$

where $ \circ $ denotes the Stratonovich integral. There is a constant $ \kappa $ for which $ \dot S_0=\lim_{m\to0}(\dot S_m-\kappa/m) $ exists.

Let $ (x_t,y_t) $ be the limit in law, as $ m\to0 $, of the stationary position process. Let $ \mathbb Q_\tau $ be the law of the path $ (x_t\bmod2\pi)_{0\le t\le\tau} $, let $ \mathbb Q_\tau^R $ be the law of the reversed path $ (x_{\tau-t}\bmod2\pi)_{0\le t\le\tau} $, and let $ \sigma_x=\lim_{\tau\to\infty}\tau^{-1}D_{\rm KL}(\mathbb Q_\tau\Vert\mathbb Q_\tau^R) $, where $ D_{\rm KL}(P\Vert Q)=\mathbb E_P\ln(dP/dQ) $.

Find $ \Delta=\dot S_0-\sigma_x $. It is a rational number.

Write `/app/output/answer.json` containing a JSON object with one field, `Delta`, holding $ \Delta $ as an exact fraction in a string of the form `"p/q"` with integers $ p $ and $ q $. It is accepted only if it equals the exact value.

You have 9000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
