import numpy as np, pickle, sys
from root_opts import *
from evalS import dvals, F as FF, T as TT
f1=F[np.flatnonzero(T==1)[0]]
res,cheap=enum_f(f1,44,maxd=3)
print('D7 options',[(len(s).bit_length(),pc(m)) for s,m in list(res.items())[:3]], len(res))
rootflats=F[T==8]; d7flats=F[T==1]
best=0
import collections
hist=collections.Counter()
for s,m in res.items():
    S0=np.zeros(120,bool); S0[[t for t in range(120) if m>>t&1]]=True
    # greedy extension: evaluate E on S0 itself and with 2 extra lines chosen among lines that complete kernels
    d=dvals(S0,rootflats); d7=dvals(S0,d7flats)
    E=int(np.sort(np.maximum(d-1,0))[::-1][:8].sum()); n7=int((d7>=2).sum())
    hist[(pc(m),E,n7)]+=1
print(sorted(hist.items()))
