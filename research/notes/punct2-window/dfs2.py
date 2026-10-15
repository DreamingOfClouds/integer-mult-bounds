import numpy as np, pickle, sys, collections, time
sys.setrecursionlimit(10000)
from setup_perm import *
exec(open('family.py').read().split("cap=int")[0])
(adj,)=pickle.load(open('adj44.pkl','rb'))
def toint(row): return int(row[0])|(int(row[1])<<64)
pc=lambda x: bin(x).count('1')
thr={25:2,26:3,27:3,28:3,29:3,30:4,31:4,32:4,33:4,34:5,35:5,36:5,37:6,38:6,39:6,40:6,41:7,42:7,43:7,44:7}
viol={}   # (frozenset of (flat,w)) -> union
maxsize=0
n=len(PK[54]); t0=time.time()
for i in range(n):
    m0=toint(PK[54][i]); w0=int(W[i])
    cands=[(b,j,toint(PK[b][j]),int(W[j])) for b,j in adj[i]]
    def rec(start,used,U,wt,fam):
        global maxsize
        u=pc(U)
        maxsize=max(maxsize,len(fam))
        if u>=25 and wt>thr[u]:
            viol[(frozenset(fam),U)]=u
        for k in range(start,len(cands)):
            b,j,m,w=cands[k]
            if b in used: continue
            nU=U|m
            if pc(nU)>44: continue
            rec(k+1,used|{b},nU,wt+w,fam+[(b,w)])
    rec(0,{54},m0,w0,[(54,w0)])
print('maxsize',maxsize,'violating nodes',len(viol),round(time.time()-t0))
pickle.dump(viol,open('viol.pkl','wb'))
c=collections.Counter((u,sum(w for b,w in f)) for (f,U),u in viol.items())
print(sorted(c.items()))
