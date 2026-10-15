import numpy as np, pickle, time
from root_opts import *
f0=F[np.flatnonzero(T==8)[0]]
res,cheap=enum_f(f0,44,maxd=4)
# kernels as masks with dims
byd={}
for s_,m in res.items():
    d=len(s_).bit_length(); byd.setdefault(d,set()).add(m)
out={}
for d,ms in byd.items():
    ms=sorted(ms,key=pc)
    mins=[]
    for m in ms:
        if not any((k&m)==k for k in mins): mins.append(m)
    # also dominance by higher d
    out[d]=mins
    print('d',d,'distinct kernels',len(ms),'inclusion-minimal',len(mins),sorted(set(pc(m) for m in mins)))
# dominate d=2 by d>=3 kernels
pickle.dump((int(np.flatnonzero(T==8)[0]),out),open('minopts_root.pkl','wb'))
