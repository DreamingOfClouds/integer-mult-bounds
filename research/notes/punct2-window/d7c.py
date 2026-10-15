import numpy as np, collections
from root_opts import enum_f, F, T, pc
idx=np.flatnonzero(T==1)
masks=set()
per=collections.Counter()
for i in idx:
    r,_=enum_f(F[i],44,maxd=3)
    for s,m in r.items(): masks.add(m)
masks=sorted(masks); print('distinct D7 kernels',len(masks))
A=np.array([[m&((1<<64)-1), m>>64] for m in masks],dtype=np.uint64)
c=collections.Counter()
for k in range(len(masks)):
    out=np.bitwise_count(A[:,0]&~A[k,0])+np.bitwise_count(A[:,1]&~A[k,1])
    c[int((out<=2).sum())]+=1
print(sorted(c.items()))
