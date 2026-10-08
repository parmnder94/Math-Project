import numpy as np, sys
from edge_engine import *
B,lam=1.5,0.6
z=complex(sys.argv[1]);A=Point(B,lam,z);Ab=Point(B,lam,np.conj(z))
J=20;rad=0.3;ts=rad*np.exp(2j*np.pi*np.arange(J)/J)
class P: pass
def wA(t):
    q=P();q.K=A.K;q.Phi=A.Phi*t;q.Om=A.Om;return q
Lb=6
for (BB,CC,iB,iC) in [(A,A,False,False),(A,A,True,False),(A,Ab,False,True),(Ab,A,True,True)]:
    vals=np.array([Psi3(wA(t),BB,CC,iB,iC) for t in ts]);c=np.array([np.mean(vals*ts**(-l)) for l in range(J)])
    b=brute_Psi3(A,BB,CC,Lb,iB,iC,bylen=True)
    print(f"Psi3 invB={iB} invC={iC}: per-length max diff {max(abs(c[l]-b[l]) for l in range(1,Lb+1)):.2e} scale {abs(b[1]):.2e}",flush=True)
