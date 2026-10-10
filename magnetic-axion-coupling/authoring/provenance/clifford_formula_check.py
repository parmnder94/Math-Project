"""Literature check at the Neel-only point (J1 = J2 = 0), where H has the Clifford form sum_a d_a Gamma_a: the closed
formula of Li, Wang, Qi and Zhang, Nat. Phys. 6, 284 (2010) against the Chern-Simons code. usage: python3 clifford_formula_check.py"""
import sys, numpy as np, itertools
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from model import PARAMS
from cs_stretched import theta_cs
p=dict(PARAMS,J1=0.0,J2=0.0)
N=128; g=(np.arange(N)+0.5)*2*np.pi/N
KX,KY,KZ=np.meshgrid(g,g,g,indexing='ij')
M0,A1,A2,B1,B2,m=[p[k] for k in ("M0","A1","A2","B1","B2","m")]
d={1:A2*np.sin(KX),2:A2*np.sin(KY),3:A1*np.sin(KZ),4:M0-2*B1*(1-np.cos(KZ))-2*B2*(2-np.cos(KX)-np.cos(KY)),5:m+0*KX}
dd={1:[A2*np.cos(KX),0*KX,0*KX],2:[0*KX,A2*np.cos(KY),0*KX],3:[0*KX,0*KX,A1*np.cos(KZ)],
    4:[-2*B2*np.sin(KX),-2*B2*np.sin(KY),-2*B1*np.sin(KZ)],5:[0*KX,0*KX,0*KX]}
nd=np.sqrt(sum(d[a]**2 for a in d))
def formula(special, others):
    eps=0
    for perm in itertools.permutations(others):
        s=np.linalg.det(np.eye(4)[[others.index(x) for x in perm]])
        i,j,k,l=perm
        eps=eps+s*d[i]*dd[j][0]*dd[k][1]*dd[l][2]
    ds=d[special]
    integrand=(2*nd+ds)/((nd+ds)**2*nd**3)*eps
    return integrand.sum()*(2*np.pi/N)**3/(4*np.pi)
ref=theta_cs(p,128,0.8)[0]
print("CS code at J=0: %.10f"%np.angle(np.exp(1j*ref)))
for special,others in [(4,[1,2,3,5]),(5,[1,2,3,4])]:
    v=formula(special,others)
    print("LWQZ formula, special d%d: %.10f  (wrapped +: %.10f, -: %.10f)"%(special,v,np.angle(np.exp(1j*v)),np.angle(np.exp(-1j*v))))
