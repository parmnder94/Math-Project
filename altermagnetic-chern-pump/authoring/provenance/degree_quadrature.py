# Degree of d/|d| on the x-torus by periodic trapezoid quadrature of det[d, d_1 d, ..., d_4 d]/|d|^5.
# usage: python3 degree_quadrature.py <grid points per axis>   (results in authoring/evidence/quadrature_attempts.json)
import numpy as np, sys, math
w=np.exp(2j*np.pi/3)
def dvec(kx,ky,kz,ph):
    l=[kx+ky,ky+kz,kz+ph,-kx-ky-kz+ph]
    p=[np.cos(3*x) for x in l]; u=[np.exp(2j*x) for x in l]
    chi=lambda a,b:(1+a+b-a*b)/2; S=lambda a,v:(1+a+(1-a)*v)/2
    Xi=lambda a,b:((1+b+b*b)+a*(1+w*w*b+w*b*b)+a*a*(1+w*b+w*w*b*b))/3
    pc=p+[p[0]]
    Z1=Xi(u[0],u[1])*Xi(u[2],u[3])*S(p[0],u[0])*S(p[1],u[1])*S(p[2],u[2])*S(p[3],u[3])
    Z2=Xi(u[0],u[2])*Xi(u[1],u[3])*S(pc[1],u[0])*S(pc[2],u[1])*S(pc[3],u[2])*S(pc[4],u[3])
    d5=chi(pc[0],pc[1])*chi(pc[1],pc[2])*chi(pc[2],pc[3])*chi(pc[3],pc[4])*(0.2+Z1.real+0.6*Z2.real)
    return np.stack([np.sin(3*x) for x in l]+[d5])
N=int(sys.argv[1]); h=1e-5
g=(np.arange(N)+0.37)*2*np.pi/N
tot=0.0; mind=1e9
for i0 in range(N):
    X=np.meshgrid(g[i0:i0+1],g,g,g,indexing='ij'); X=[x.ravel() for x in X]
    d=dvec(*X); mind=min(mind,np.sqrt((d**2).sum(0)).min())
    J=[]
    for k in range(4):
        Xp=[x.copy() for x in X]; Xm=[x.copy() for x in X]; Xp[k]+=h; Xm[k]-=h
        J.append((dvec(*Xp)-dvec(*Xm))/(2*h))
    M=np.stack([d]+J,axis=0).transpose(2,1,0)   # (pts, comp, col)
    det=np.linalg.det(M)
    tot+=(det/np.sqrt((d**2).sum(0))**5).sum()
vol=(2*np.pi/N)**4
print(N,"degree",tot*vol/(8*np.pi**2/3),"min|d| on grid",mind)
