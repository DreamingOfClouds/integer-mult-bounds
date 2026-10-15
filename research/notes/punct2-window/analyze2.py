import numpy as np, pickle, collections, sys, math, time
from analyze import best_D, Dwin, G
from evalS import dvals, F, T
viol=pickle.load(open('viol.pkl','rb'))
unions={}
for (fam,Um),u in viol.items(): unions[Um]=u
print('distinct violating unions',len(unions),collections.Counter(unions.values()))
rootflats=F[T==8]
# line id of root flat k: need mapping; flats order is BFS; find line alpha with orth-pattern == ~flat
rootline=[]
for f in rootflats:
    c=[t for t in range(120) if (G[t]==0).tolist()==(~f).tolist()]
    assert len(c)==1; rootline.append(c[0])
rootline=np.array(rootline)
res=[]
t0=time.time()
for Um,u in unions.items():
    S=np.zeros(120,bool); S[[t for t in range(120) if Um>>t&1]]=True
    d=dvals(S,rootflats)
    R=[int(rootline[k]) for k in range(120) if d[k]>=2]; dval={int(rootline[k]):int(d[k]) for k in range(120) if d[k]>=2}
    bd=best_D(R,dval)
    res.append((u,len(R),sorted(dval.values()),bd,Dwin(u)))
print(round(time.time()-t0))
c=collections.Counter((u,nr,tuple(dv),bd,dw) for u,nr,dv,bd,dw in res)
for k,v in sorted(c.items()): print(k,v)
pickle.dump(res,open('res2.pkl','wb'))
