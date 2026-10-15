import numpy as np, pickle, sys, collections
sys.setrecursionlimit(10000)
from setup_perm import *
exec(open('family.py').read().split("cap=int")[0])   # builds PK, W, base
(adj,)=pickle.load(open('adj44.pkl','rb'))
def toint(row): return int(row[0])|(int(row[1])<<64)
pc=lambda x: bin(x).count('1')
best=collections.defaultdict(int); arg={}
sizes=collections.Counter(len(a) for a in adj)
print('max cands',max(len(a) for a in adj))
n=len(PK[54])
for i in range(n):
    m0=toint(PK[54][i]); w0=int(W[i])
    cands=[(b,j,toint(PK[b][j]),int(W[j])) for b,j in adj[i]]
    def rec(start,used,U,wt,fam):
        u=pc(U)
        if wt>best[u]: best[u]=wt; arg[u]=list(fam)
        if len(fam)>=8: return
        for k in range(start,len(cands)):
            b,j,m,w=cands[k]
            if b in used: continue
            nU=U|m
            if pc(nU)>44: continue
            rec(k+1,used|{b},nU,wt+w,fam+[(b,j,w)])
    rec(0,{54},m0,w0,[(54,i,w0)])
# cumulative
cum=0
for u in range(25,45):
    cum=max(cum,best.get(u,0)); print('union<=',u,'max E (d=2 counts 1, d=3 counts 2) =',cum)
pickle.dump((dict(best),arg),open('dfs_res.pkl','wb'))
