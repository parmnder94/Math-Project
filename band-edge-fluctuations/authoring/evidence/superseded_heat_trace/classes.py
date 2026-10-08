import numpy as np, itertools, time
from words import INV, s as wstr
def cyc_reduced_words(l):
    out=[]
    def rec(w):
        if len(w)==l:
            if w[-1]!=INV[w[0]]: out.append(tuple(w))
            return
        for x in range(4):
            if w and x==INV[w[-1]]: continue
            rec(w+[x])
    rec([]);return out
def is_power(w):
    l=len(w)
    for p in range(1,l):
        if l%p==0 and w==w[:p]*(l//p): return True
    return False
def rotations(w): return [w[k:]+w[:k] for k in range(len(w))]
def winv(w): return tuple(INV[x] for x in reversed(w))
def canon(w): return min(rotations(w)+rotations(winv(w)))
def primitive_classes(lmax):
    cls=[]
    for l in range(1,lmax+1):
        for w in cyc_reduced_words(l):
            if not is_power(w) and canon(w)==w: cls.append(w)
    return cls
def ndiv(n): return sum(1 for d in range(1,n+1) if n%d==0)
if __name__=="__main__":
    for l in range(1,9):
        print(l,len([c for c in primitive_classes(l) if len(c)==l]))
