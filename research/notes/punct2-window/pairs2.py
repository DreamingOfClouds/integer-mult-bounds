import numpy as np, pickle, collections, time
from setup_perm import *
tr=transporters()
a0,out=pickle.load(open('minopts_root.pkl','rb')); alpha=54
Gram=np.array(U.Gram)
def pack(ms):
    A=np.zeros((len(ms),2),dtype=np.uint64)
    for i,m in enumerate(ms):
        A[i,0]=m&((1<<64)-1); A[i,1]=m>>64
    return A
def img(ms,g): return [apply_perm(m,g) for m in ms]
print(hasattr(np,'bitwise_count'))
for dA,dB in ((2,2),):
    base=out[2]
    PA=pack(base)
    for name,beta in (('orth',next(t for t in range(120) if Gram[54][t]==0 and t!=54)),('nonorth',next(t for t in range(120) if Gram[54][t]!=0 and t!=54))):
        PB=pack(img(base,tr[beta]))
        hist=collections.Counter()
        for i in range(len(PA)):
            u=PA[i]|PB    # (nB,2)
            c=np.bitwise_count(u).sum(1)
            hist.update(c.tolist())
        print(name,beta,'min union',min(hist),sorted(hist.items())[:12],flush=True)
