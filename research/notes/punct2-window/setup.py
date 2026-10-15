import sys; sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np, itertools, pickle
from bit import *
BU=BitUnit('E8'); U=BU.U; v=120; L=np.array(U.L,dtype=np.int64)
# reflection perms on lines (doubled coords norm 8 -> reflect: x - 2(x.r)/(r.r) r , r.r=8)
idx={}
for i,l in enumerate(U.L): idx[tuple(l)]=i; idx[tuple(-np.array(l))]=i
def refl_perm(r):
    r=np.array(r); P=[]
    for l in L:
        y=l-(2*(l@r)//8)*r if (2*(l@r))%8==0 else None
        assert y is not None
        P.append(idx[tuple(y)])
    return np.array(P)
perms=[refl_perm(r) for r in U.L]
sr=[(1,-1,-1,-1,-1,-1,-1,1),(2,2,0,0,0,0,0,0),(-2,2,0,0,0,0,0,0),(0,-2,2,0,0,0,0,0),(0,0,-2,2,0,0,0,0),(0,0,0,-2,2,0,0,0),(0,0,0,0,-2,2,0,0),(0,0,0,0,0,-2,2,0)]
A=np.array(sr,dtype=float); Winv=np.linalg.inv(A)
reps=[]
for i in range(8):
    s=Winv[:,i]
    for m in range(1,64):
        if np.allclose(s*m,np.round(s*m)): s=np.round(s*m).astype(np.int64); break
    reps.append((L@s)!=0)   # complement mask
simple=[refl_perm(np.array(r)) for r in sr]
allflats=[]; types=[]
for ti,rep in enumerate(reps):
    cur=np.packbits(rep.astype(np.uint8))[None,:]
    seen={rep.tobytes()}
    frontier=rep[None,:]
    allm=[rep]
    while len(frontier):
        new=[]
        for p in simple:
            nf=frontier[:,np.argsort(p)]  # image set under perm: mask'[p[t]]=mask[t]
            new.append(nf)
        nf=np.concatenate(new); 
        pk=np.ascontiguousarray(np.packbits(nf,axis=1))
        _,ui=np.unique(pk.view(np.dtype((np.void,pk.shape[1]))),return_index=True)
        nf=nf[ui]
        keep=[k for k,row in enumerate(nf) if row.tobytes() not in seen]
        nf=nf[keep]
        for row in nf: seen.add(row.tobytes())
        frontier=nf
        if len(nf): allm.append(nf)
    arr=np.concatenate([a.reshape(-1,120) for a in allm])
    print(ti+1,'compl',rep.sum(),'orbit',len(arr),flush=True)
    allflats.append(arr); types+= [ti+1]*len(arr)
F=np.concatenate(allflats); T=np.array(types)
np.save('flats.npy',F); np.save('ftypes.npy',T); np.save('Cm.npy',BU.Cm); np.save('Om.npy',BU.Om)
print(F.shape)
