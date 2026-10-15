import numpy as np, random
from core import *
U=Unit('E8'); G=np.array(U.Gram); Cm=np.load('Cm.npy')
def rank2(vs):
    return rank_f2([int(sum(int(b)<<i for i,b in enumerate(c))) for c in vs])
random.seed(1); viol=0; tight=0
for trial in range(3000):
    S=set(random.sample(range(120),random.choice([0,5,16,25,40])))
    live=[t for t in range(120) if t not in S]
    a,b=random.randrange(120),random.randrange(120)
    R=[t for t in live if G[a][t]!=0]; A=[t for t in live if G[b][t]!=0]
    rows=[sum(1<<s for s in A if G[t][s]==0 or s==t) for t in R]
    r=rank_f2(rows)
    dx=9-rank2([Cm[t] for t in R]); dy=9-rank2([Cm[t] for t in A])
    if r<9-dx-dy: viol+=1
    if r==9-dx-dy: tight+=1
print('violations',viol,'tight',tight)
