import numpy as np, pickle, collections, time, sys
from setup_perm import *
tr=transporters(54)
a0,out=pickle.load(open('minopts_root.pkl','rb'))
base=[(m,1) for m in out[2]]+[(m,2) for m in out[3]]
def pack(ms):
    A=np.zeros((len(ms),2),dtype=np.uint64)
    for i,m in enumerate(ms):
        A[i,0]=m&((1<<64)-1); A[i,1]=m>>64
    return A
# permutation action on packed masks via bit arrays
bits=np.zeros((len(base),120),dtype=bool)
for i,(m,w) in enumerate(base):
    for t in range(120):
        if m>>t&1: bits[i,t]=True
W=np.array([w for m,w in base])
def packbits(b):
    lo=np.zeros(len(b),dtype=np.uint64); hi=np.zeros(len(b),dtype=np.uint64)
    for t in range(64): lo|=b[:,t].astype(np.uint64)<<np.uint64(t)
    for t in range(64,120): hi|=b[:,t].astype(np.uint64)<<np.uint64(t-64)
    return np.stack([lo,hi],1)
PK={}   # flat(root line) -> packed options
for beta in range(120):
    g=tr[beta]
    nb=np.zeros_like(bits); nb[:,g]=bits   # image: line t -> g[t]
    PK[beta]=packbits(nb)
cap=int(sys.argv[1]) if len(sys.argv)>1 else 44
PA=PK[54]; n=len(PA)
# adjacency: for each option i of flat 54, list of (beta,j) with union<=cap
adj=[[] for _ in range(n)]
t0=time.time()
for beta in range(120):
    if beta==54: continue
    PB=PK[beta]
    for a in range(0,n,512):
        u=np.bitwise_count(PA[a:a+512,None,:]|PB[None,:,:]).sum(2)
        ii,jj=np.nonzero(u<=cap)
        for i,j in zip(ii,jj): adj[a+i].append((beta,int(j)))
    if beta%20==0: print('beta',beta,round(time.time()-t0),flush=True)
pickle.dump((adj,),open('adj%d.pkl'%cap,'wb'))
print('total adj',sum(len(x) for x in adj))
