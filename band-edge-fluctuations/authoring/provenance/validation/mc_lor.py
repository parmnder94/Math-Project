import numpy as np, sys, time
from cactus import R,sx,sy,sz
lam=0.6;B=1.5;x0=float(sys.argv[1]);eta=float(sys.argv[2])
N=int(sys.argv[3]);S=int(sys.argv[4]);seed=int(sys.argv[5])
rng=np.random.default_rng(seed);Ts=[R(sx,lam),R(sy,lam),R(sz,lam)]
P11=B*sz+sum(T+T.conj().T for T in Ts);old=np.linalg.eigvalsh(P11)
h=lambda x: eta/((x-x0)**2+eta**2)
out=[];t0=time.time()
for t in range(S):
    sig=rng.permutation(N);tau=rng.permutation(N)
    Pa=np.zeros((N,N));Pa[sig,np.arange(N)]=1;Pb=np.zeros((N,N));Pb[tau,np.arange(N)]=1;Pc=Pa@Pb
    H=np.kron(B*sz,np.eye(N))
    for T,Pm in zip(Ts,[Pa,Pb,Pc]): H=H+np.kron(T,Pm)+np.kron(T.conj().T,Pm.T)
    ev=np.linalg.eigvalsh(H)
    out.append(h(ev).sum()-h(old).sum())
    if t%2000==0: np.save(f"mclor_{x0}_{eta}_N{N}_s{seed}.npy",np.array(out))
np.save(f"mclor_{x0}_{eta}_N{N}_s{seed}.npy",np.array(out));print("done",time.time()-t0)
