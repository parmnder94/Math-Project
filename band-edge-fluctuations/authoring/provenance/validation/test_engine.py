import numpy as np,time
from edge_engine import *
B,lam=1.5,0.6
import sys; z=complex(sys.argv[1])
A=Point(B,lam,z);Ab=Point(B,lam,np.conj(z))
print("conj-sum spr",A.conj_spr)
for d in [2,3]:
    t0=time.time();e=M_d(A,d);b=brute_M_d(A,d,9);print(f"M_{d}: transfer {e:.12e}  brute(l<=9) {b:.12e}  diff {abs(e-b):.1e}  t={time.time()-t0:.1f}")
for (p1,p2,inv,BB) in [(1,1,False,A),(1,1,True,A),(1,1,False,Ab),(1,1,True,Ab),(1,2,False,Ab),(2,1,True,A)]:
    t0=time.time();e=Psi2(A,BB,p1,p2,inv);b=brute_Psi2(A,BB,p1,p2,8,inv)
    print(f"Psi2 p=({p1},{p2}) inv={inv} B={'zbar' if BB is Ab else 'z'}: transfer {e:.10e} brute(l<=8) {b:.10e} diff {abs(e-b):.1e} t={time.time()-t0:.1f}")
