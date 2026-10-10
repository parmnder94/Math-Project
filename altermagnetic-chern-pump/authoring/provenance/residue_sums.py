# Floating-point cross-check of the golden-solution intermediates: north-pole zero count on the l-torus,
# det L, and the binary/ternary character sums W_0, W_A, W_B, W_AB and U_rs (imaginary parts printed in units of sqrt 3).
import itertools, cmath, math
import numpy as np
w=cmath.exp(2j*math.pi/3)
def chi(p,q): return (1+p+q-p*q)/2
def S(p,u): return (1+p+(1-p)*u)/2
def Xi(u,v): return ((1+v+v*v)+u*(1+w*w*v+w*v*v)+u*u*(1+w*v+w*w*v*v))/3
def d5(l):
    p=[math.cos(3*x) for x in l]; u=[cmath.exp(2j*x) for x in l]
    pc=p+[p[0]]
    Z1=Xi(u[0],u[1])*Xi(u[2],u[3])*np.prod([S(p[i],u[i]) for i in range(4)])
    Z2=Xi(u[0],u[2])*Xi(u[1],u[3])*np.prod([S(pc[i+1],u[i]) for i in range(4)])
    return np.prod([chi(pc[i],pc[i+1]) for i in range(4)])*(0.2+Z1.real+0.6*Z2.real)
# degree on l-torus: north-pole preimages
D=0; mn=1e9; sumsigned=0
for n in itertools.product(range(6),repeat=4):
    l=[math.pi*k/3 for k in n]
    v=d5(l); mn=min(mn,abs(v))
    c=np.prod([math.cos(3*x) for x in l])  # sign of Jacobian (3^4 * prod cos)
    s=1 if c>0 else -1
    if v>0: D+=s
    sumsigned+=s*np.sign(v)
print("min|d5| at zeros",mn,"D_l(north)",D,"signed sum/2",sumsigned/2)
L=np.array([[1,1,0,0],[0,1,1,0],[0,0,1,1],[-1,-1,-1,1]]); print("detL",round(np.linalg.det(L)))
# intermediate sums
def A(a,b): return (a[0]*a[1]+a[2]*a[3]+sum(a[i]*b[i] for i in range(4)))%3
def B(a,b): return (a[0]*a[2]+a[1]*a[3]+a[0]*b[1]+a[1]*b[2]+a[2]*b[3]+a[3]*b[0])%3
def wt(b): return (-1)**(sum(b)+b[0]*b[1]+b[1]*b[2]+b[2]*b[3]+b[3]*b[0])
W0=WA=WB=WAB=0; U={}
for b in itertools.product(range(2),repeat=4):
    W0+=wt(b)
    for a in itertools.product(range(3),repeat=4):
        x,y=A(a,b),B(a,b)
        WA+=wt(b)*(x==0); WB+=wt(b)*(y==0); WAB+=wt(b)*(x==0 and y==0)
        for r in range(3):
            for s in range(3): U[(r,s)]=U.get((r,s),0)+wt(b)*w**((r*x+s*y)%3)
print("W0",W0,"WA",WA,"WB",WB,"WAB",WAB)
for k,v in sorted(U.items()): print(k, complex(round(v.real,6),round(v.imag/math.sqrt(3),6)),"(imag in units of sqrt3)")
