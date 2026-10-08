# K=4: per-length validation of M_d = sum_{y cyc. red.} F(y^d) (Kronecker power + cyclic closure), d=2,3,4
import numpy as np
from engine4 import Point, M_d, cr_words
B,lam=1.5,0.6;z=5.63+0.3j;A=Point(B,lam,z)
J=20;rad=0.3;ts=rad*np.exp(2j*np.pi*np.arange(J)/J)
class P_: pass
for d in (2,3,4):
    def f(t):
        q=P_();q.K=A.K;q.Phi=A.Phi*t**(1.0/d);q.Om=A.Om;return M_d(q,d)
    vals=np.array([f(t) for t in ts]);c=np.array([np.mean(vals*ts**(-l)) for l in range(J)])
    Lb=7;b=np.zeros(Lb+1,complex)
    for y in cr_words(Lb): b[len(y)]+=A.F(y*d)
    print(f"M_{d} (K=4): per-length max diff {max(abs(c[l]-b[l]) for l in range(1,Lb+1)):.2e} scale {abs(b[1]):.2e}")
