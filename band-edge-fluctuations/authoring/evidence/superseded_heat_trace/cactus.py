import numpy as np
sx=np.array([[0,1],[1,0]],complex); sy=np.array([[0,-1j],[1j,0]]); sz=np.diag([1.,-1]).astype(complex); I2=np.eye(2,dtype=complex)
def R(n,th): return np.cos(th)*I2+1j*np.sin(th)*n
Z2=np.zeros((2,2),complex)
class Cactus:
    def __init__(s,Ta,Tb,Tc,D):
        s.Ta,s.Tb,s.Tc,s.D=Ta,Tb,Tc,D
        d=lambda M:M.conj().T
        s.C=[np.vstack([Tb,Tc]), np.vstack([d(Tb),Ta]), np.vstack([d(Tc),d(Ta)])]
        s.hp=[np.block([[Z2,d(Ta)],[Ta,Z2]]), np.block([[Z2,d(Tc)],[Tc,Z2]]), np.block([[Z2,d(Tb)],[Tb,Z2]])]
        s.others=[(1,2),(0,2),(0,1)]
    def F(s,Sig,z):
        St=Sig[0]+Sig[1]+Sig[2]; base=z*I2-s.D-St
        out=[]
        for r in range(3):
            r1,r2=s.others[r]
            M=np.block([[base+Sig[r1],Z2],[Z2,base+Sig[r2]]])-s.hp[r]
            out.append(s.C[r].conj().T@np.linalg.solve(M,s.C[r]))
        return out
    def G(s,Sig,z): return np.linalg.inv(z*I2-s.D-sum(Sig))
    def solve(s,z,Sig0=None,iters=100000,tol=1e-15,damp=0.5):
        Sig=Sig0 if Sig0 is not None else [ -1j*I2*0.5 for _ in range(3)]
        for it in range(iters):
            N=s.F(Sig,z); err=max(np.abs(N[r]-Sig[r]).max() for r in range(3))
            Sig=[damp*Sig[r]+(1-damp)*N[r] for r in range(3)]
            if err<tol: break
        return Sig,it
    # Newton on complex unknowns (12 complex)
    def newton(s,z,Sig,tol=1e-14,maxit=50):
        def vec(S): return np.concatenate([m.ravel() for m in S])
        def mat(v): return [v[4*r:4*r+4].reshape(2,2) for r in range(3)]
        v=vec(Sig)
        for it in range(maxit):
            r0=vec(s.F(mat(v),z))-v
            if np.abs(r0).max()<tol: break
            J=np.zeros((12,12),complex);h=1e-7
            for j in range(12):
                vv=v.copy();vv[j]+=h;J[:,j]=(vec(s.F(mat(vv),z))-vv-r0)/h   # holomorphic in entries
            v=v-np.linalg.solve(J,r0)
        return mat(v),np.abs(vec(s.F(mat(v),z))-v).max()
def modelA(B,lam,t=1.0):
    return Cactus(t*R(sx,lam),t*R(sy,lam),t*R(sz,lam),B*sz)
def rho(c,x,eta=1e-10,Sig0=None):
    # continuation in eta from 1e-1 down
    Sig=Sig0
    if Sig is None:
        Sig,_=c.solve(x+0.5j,iters=20000,tol=1e-13)
        for e in [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,eta]:
            Sig,res=c.newton(x+1j*e,Sig)
    else:
        Sig,res=c.newton(x+1j*eta,Sig)
    g=c.G(Sig,x+1j*eta);return -np.trace(g).imag/(2*np.pi),Sig
