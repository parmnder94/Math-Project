import numpy as np, time
from engine4 import Point, brute_marked
from engine4b import psi_marked
B,lam=1.5,0.6;z=5.63+0.3j
A=Point(B,lam,z);Ab=Point(B,lam,np.conj(z))
J=20;rad=0.3;ts=rad*np.exp(2j*np.pi*np.arange(J)/J)
class P_: pass
def wA(pt,t,p0):
    q=P_();q.K=pt.K;q.Phi=pt.Phi*t**(1.0/p0);q.Om=pt.Om;return q
Lb=6
cases=[[(A,False,2),(A,True,2)],[(A,False,2),(A,False,1),(Ab,True,1)],[(Ab,False,1),(A,True,1),(Ab,True,2)],
       [(A,False,1),(A,False,3)],[(Ab,False,3),(A,True,1)],[(A,False,1),(Ab,True,1),(A,False,1)]]
for cs in cases:
    t0=time.time();P0,i0,p0=cs[0]
    vals=np.array([psi_marked([(wA(P0,t,p0),i0,p0)]+cs[1:]) for t in ts]);c=np.array([np.mean(vals*ts**(-l)) for l in range(J)])
    b=brute_marked(cs,Lb)
    print([("zb" if P is Ab else "z",inv,p) for P,inv,p in cs], f"per-length max diff {max(abs(c[l]-b[l]) for l in range(1,Lb+1)):.2e} scale {abs(b[1]):.2e} t/eval {(time.time()-t0)/J:.2f}s",flush=True)
