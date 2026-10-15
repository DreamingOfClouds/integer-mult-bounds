import numpy as np, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from core import *
U=Unit('E8'); L=np.array(U.L,dtype=np.int64)
idx={}
for i,l in enumerate(U.L): idx[tuple(l)]=i; idx[tuple(-np.array(l))]=i
def refl_perm(r):
    r=np.array(r)
    return np.array([idx[tuple(l-(2*(l@r)//8)*r)] for l in L])
perms=[refl_perm(r) for r in U.L]       # reflection in line i
def mask_to_set(m): return [t for t in range(120) if m>>t&1]
def apply_perm(m,p):
    out=0
    for t in range(120):
        if m>>t&1: out|=1<<int(p[t])
    return out
# group elements mapping line 0 -> every line (BFS)
def transporters(start=54):
    g={start:np.arange(120)}
    fr=[start]
    while fr:
        nf=[]
        for a in fr:
            for r in range(120):
                b=int(perms[r][a])
                if b not in g:
                    g[b]=perms[r][g[a]]; nf.append(b)   # g[b]= refl_r o g[a]
        fr=nf
    return g
