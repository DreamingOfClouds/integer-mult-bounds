import numpy as np, pickle
from setup_perm import *
F=np.load('flats.npy'); T=np.load('ftypes.npy'); Cm=np.load('Cm.npy')
Zs=np.array([[(z>>i)&1 for i in range(9)] for z in range(1,512)])
Par=((Zs@Cm.T)%2).astype(np.float32)           # 511x120
Gram=np.array(U.Gram)
def dvals(Smask, flats):
    """Smask: bool(120) punctured. returns deficiency d for each flat row (bool complement masks)."""
    live=(~Smask).astype(np.float32)
    # kernel active iff no live line in C_f with z=1:  (F*live) @ Par.T == 0
    W=((flats.astype(np.float32)*live)@Par.T)
    act=(W==0).sum(1)           # n_f = number of active z
    return np.floor(np.log2(act+1)+1e-9).astype(int)
if __name__=='__main__':
    a0,out=pickle.load(open('minopts_root.pkl','rb'))
    tr=transporters()
    rootflats=F[T==8]
    ks=out[2]; 
    k25=[m for m in ks if bin(m).count('1')==25]
    for m in k25[:3]:
        S=np.zeros(120,bool); S[mask_to_set(m)]=True
        d=dvals(S,rootflats)
        print('|S|',S.sum(),'root flats d>=1:',(d>=1).sum(),'d>=2:',(d>=2).sum(),'d>=3',(d>=3).sum())
        dall=dvals(S,F)
        for t in range(1,9): print(' type',t,'d>=1',(dall[T==t]>=1).sum(),'d>=2',(dall[T==t]>=2).sum())
