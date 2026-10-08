import numpy as np, sys
import edge_engine as E
from edge_engine import *
B,lam=1.5,0.6
z=complex(sys.argv[1]);A=Point(B,lam,z);Ab=Point(B,lam,np.conj(z))
J=32;rad=0.35;ts=rad*np.exp(2j*np.pi*np.arange(J)/J)
def weighted(pt,t):
    class P: pass
    q=P();q.K=pt.K;q.Phi=pt.Phi*np.sqrt(t);q.Om=pt.Om;return q
def coeffs(fun):
    vals=np.array([fun(t) for t in ts]);return np.array([np.mean(vals*ts**(-l)) for l in range(J)])
def brute_by_len(fun_word,lmax):
    out=np.zeros(lmax+1,complex)
    for r in cr_words(lmax): out[len(r)]+=fun_word(r)
    return out
Lb=7
# M_2
c=coeffs(lambda t: M_d(weighted(A,t),2))
b=brute_by_len(lambda y:A.F(y*2),Lb)
print("M_2 per-length max diff",max(abs(c[l]-b[l]) for l in range(1,Lb+1)), "scale",abs(b[1]))
for (p1,p2,inv,BB) in [(1,1,False,A),(1,1,True,A),(1,1,False,Ab),(1,1,True,Ab),(1,2,True,Ab)]:
    # weight only slot A's letters? Psi2 uses both slots -> each letter weighted twice; weight via A only
    def f(t):
        class P: pass
        q=P();q.K=A.K;q.Phi=A.Phi*t;q.Om=A.Om
        return Psi2(q,BB,p1,p2,inv)
    c=coeffs(f)
    def fw(r):
        fa=A.F(r*p1);tot=0
        for k in range(len(r)):
            rr=(winv(r[:k])+winv(r[k:])) if inv else r[k:]+r[:k]
            tot+=fa*BB.F(rr*p2)
        return tot
    b=brute_by_len(fw,Lb)
    print(f"Psi2 {p1},{p2} inv={inv}: per-length max diff {max(abs(c[l]-b[l]) for l in range(1,Lb+1)):.2e} scale {abs(b[1]):.2e}")
