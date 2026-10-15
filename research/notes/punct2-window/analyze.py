import numpy as np, pickle, collections, sys, math
from fractions import Fraction as Fr
sys.path.insert(0,'.')
from core import Unit
U=Unit('E8'); G=np.array(U.Gram)
viol=pickle.load(open('viol.pkl','rb'))
x=Fr(54,102)
def Dwin(p): return math.floor(72-x*(120-p))+1     # minimal D that would beat P102
def best_D(R,dval):
    """R: list of root-line ids, dval: dict id->d (>=2). maximize sum over 8 terms of dx+dy, free side=1,
    subject to pairwise necessary conditions among R-valued entries."""
    opts=[None]+list(R)
    best=[0]
    # terms chosen in canonical order: sort by (x,y) tuple index to avoid permutations
    entries=[(a,b) for a in opts for b in opts]
    key=lambda e:(-1 if e[0] is None else R.index(e[0]),-1 if e[1] is None else R.index(e[1]))
    entries.sort(key=key)
    val=lambda a:1 if a is None else dval[a]
    def rec(start,chosen,usedx,usedy,tot):
        if len(chosen)==8:
            best[0]=max(best[0],tot); return
        # optimistic bound: remaining terms at most 2*maxd
        rem=8-len(chosen)
        if tot+rem*2*max(dval.values())<=best[0]: return
        for k in range(start,len(entries)):
            a,b=entries[k]
            if a is not None and a in usedx: continue
            if b is not None and b in usedy: continue
            if a is not None and b is not None and G[a][b]==0: continue
            ok=True
            for (c,d) in chosen:
                if a is not None and d is not None and G[a][d]!=0: ok=False;break
                if c is not None and b is not None and G[c][b]!=0: ok=False;break
            if not ok: continue
            rec(k,chosen+[(a,b)],usedx|({a} if a is not None else set()),usedy|({b} if b is not None else set()),tot+val(a)+val(b))
    rec(0,[],frozenset(),frozenset(),0)
    return best[0]
if __name__=='__main__':
    groups={}
    for (fam,Um),u in viol.items():
        key=tuple(sorted(fam))
        groups.setdefault(key,u)
    print('distinct families',len(groups))
    surv=[]
    for key,u in groups.items():
        R=[b for b,w in key]; dval={b:w+1 for b,w in key}
        bd=best_D(R,dval)
        need=min(Dwin(p) for p in range(u,45))
        if bd>=need: surv.append((key,u,bd))
    print('survivors',len(surv))
    pickle.dump(surv,open('surv.pkl','wb'))
    for s in surv[:10]: print(s)
