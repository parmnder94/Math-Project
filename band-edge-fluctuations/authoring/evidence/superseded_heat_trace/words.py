import numpy as np
INV=[1,0,3,2]   # a,A,b,B
NAME='aAbB'
def reduce(w):
    out=[]
    for x in w:
        if out and out[-1]==INV[x]: out.pop()
        else: out.append(x)
    return tuple(out)
def mul(u,v): return reduce(u+v)
def inv(w): return tuple(INV[x] for x in reversed(w))
def s(w): return ''.join(NAME[x] for x in w) or 'e'
# group algebra with 2x2 matrix coefficients
def ga_mul(X,Y):
    Z={}
    for u,cu in X.items():
        for v,cv in Y.items():
            w=mul(u,v);Z[w]=Z.get(w,0)+cu@cv
    return Z
def P_elem(B,lam):
    from cactus import R,sx,sy,sz
    Ta=R(sx,lam);Tb=R(sy,lam);Tc=R(sz,lam)
    return {():B*sz,(0,):Ta,(1,):Ta.conj().T,(2,):Tb,(3,):Tb.conj().T,(0,2):Tc,(3,1):Tc.conj().T}
