import numpy as np
F=np.load('flats.npy').astype(np.float32); T=np.load('ftypes.npy'); Cm=np.load('Cm.npy')
Z=np.array([[(z>>i)&1 for i in range(9)] for z in range(1,512)])
Par=((Z@Cm.T)%2).astype(np.float32)   # 511 x120
print('|O_z| dist',np.unique(Par.sum(1),return_counts=True))
cnt={}
for a in range(0,len(F),50000):
    W=(F[a:a+50000]@Par.T).astype(np.int16)
    for t in range(1,9):
        m=T[a:a+50000]==t
        if m.any():
            w=W[m]
            for thr in (14,27,44):
                cnt[(t,thr)]=cnt.get((t,thr),0)+int((w<=thr).sum())
            cnt[(t,'min')]=min(cnt.get((t,'min'),99),int(w.min()))
for k in sorted(cnt,key=str): print(k,cnt[k])
