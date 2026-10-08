# Spectral radii of every transfer operator used by the exact resummation at the shipped point (must be < 1)
import numpy as np, itertools
from engine4 import Point, kron_all, slot_mats, INV
z=complex(5.6,0.06);P={'+':Point(1.5,0.6,z),'-':Point(1.5,0.6,np.conj(z))}
def spr(slots):
    Dl=[kron_all([sl[a] for sl in slots]) for a in range(4)];D=Dl[0].shape[0]
    T=np.zeros((4*D,4*D),complex)
    for a in range(4):
        for b in range(4):
            if b!=INV[a]: T[a*D:(a+1)*D,b*D:(b+1)*D]=Dl[b]
    return max(abs(np.linalg.eigvals(T)))
worst=0;rows=[]
print("conjugator superoperator spr:",P['+'].conj_spr,P['-'].conj_spr)
for R in range(1,5):
    for combo in itertools.product([('+',False),('-',False),('+',True),('-',True)],repeat=R):
        if combo[0][1]: continue
        s=spr([slot_mats(P[sg],inv) for sg,inv in combo]);worst=max(worst,s);rows.append((R,combo,s))
for R in range(1,5):
    print(f"tensor power {R}: max spr = {max(s for r,c,s in rows if r==R):.4f}")
print("worst over all transfers:",worst)
