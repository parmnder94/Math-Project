# Validate: E fix(w) -> #divisors(d), Cov(fix w, fix w') -> sigma(gcd(d,d')) if same root class (mod inversion), else 0
import numpy as np, sys
from words import INV
from classes import canon
rng=np.random.default_rng(1)
N=int(sys.argv[1]);S=int(sys.argv[2])
W={'a':(0,),'A':(1,),'aa':(0,0),'aaa':(0,0,0),'aaaa':(0,0,0,0),'b':(2,),'ab':(0,2),'ba':(2,0),'BA':(3,1),'abab':(0,2,0,2),'aB':(0,3),'aab':(0,0,2),'bAA':(2,1,1),'aabaab':(0,0,2,0,0,2),'abAB':(0,2,1,3),'bAbaB':(2,1,2,0,3)}
names=list(W)
def apply(w,sig,tau,sigi,taui):
    idx=np.arange(N)
    maps=[sig,sigi,tau,taui]
    for x in reversed(w): idx=maps[x][idx]   # U_w e_v = e_{w(v)}, rightmost letter acts first
    return idx
F=np.zeros((S,len(names)))
for t in range(S):
    sig=rng.permutation(N);tau=rng.permutation(N);sigi=np.argsort(sig);taui=np.argsort(tau)
    for k,nm in enumerate(names):
        F[t,k]=np.sum(apply(W[nm],sig,tau,sigi,taui)==np.arange(N))
mu=F.mean(0);C=np.cov(F.T)
def root(w):
    l=len(w)
    for p in range(1,l+1):
        if l%p==0 and w==w[:p]*(l//p): return w[:p],l//p
def cycred(w):
    w=list(w)
    while len(w)>1 and w[0]==INV[w[-1]]: w=w[1:-1]
    return tuple(w)
from math import gcd
def sigma(n): return sum(d for d in range(1,n+1) if n%d==0)
def tau(n): return sum(1 for d in range(1,n+1) if n%d==0)
info={}
for nm in names:
    u,d=root(cycred(W[nm]));info[nm]=(canon(u),d)
print("word  E[fix] theory")
for k,nm in enumerate(names): print(f"{nm:8s} {mu[k]:.3f} {tau(info[nm][1])}")
print("max |cov - theory| over pairs:")
err=0;worst=None
for i,a in enumerate(names):
    for j,b in enumerate(names):
        th=sigma(gcd(info[a][1],info[b][1])) if info[a][0]==info[b][0] else 0
        e=abs(C[i,j]-th)
        if e>err: err=e;worst=(a,b,C[i,j],th)
print(err,worst, "stat err ~",np.sqrt(2.0/S)*3)
