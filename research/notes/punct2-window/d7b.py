import numpy as np, pickle, collections
from root_opts import *
exec(open('family.py').read().split("cap=int")[0])
FL=np.load('flats.npy'); TY=np.load('ftypes.npy')
f1=FL[np.flatnonzero(TY==1)[0]]
res,cheap=enum_f(f1,44,maxd=3)
d7m=sorted(set(res.values()))
print(len(d7m))
def toint(row): return int(row[0])|(int(row[1])<<64)
allopts=np.concatenate([PK[b] for b in range(120)])   # 2M x2 uint64
flat_of=np.repeat(np.arange(120),len(PK[54])); wts=np.tile(W,120)
hist=collections.Counter()
for m in d7m:
    lo=np.uint64(m&((1<<64)-1)); hi=np.uint64(m>>64)
    out=np.bitwise_count(allopts[:,0]&~lo)+np.bitwise_count(allopts[:,1]&~hi)
    ok=np.flatnonzero(out<=2)
    hist[len(ok)]+=1
    if len(ok)>0 and hist[len(ok)]==1: print('example',len(ok), collections.Counter(zip(flat_of[ok].tolist(),out[ok].tolist())).most_common(5))
print(hist)
