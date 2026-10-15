import numpy as np, sys, pickle
F=np.load('flats.npy'); T=np.load('ftypes.npy'); Cm=np.load('Cm.npy')
Z=np.array([[(z>>i)&1 for i in range(9)] for z in range(1,512)])
Par=((Z@Cm.T)%2).astype(bool)   # row z-1
def kern_masks(f):   # f: complement mask bool(120); per z kernel = Par[z]&f as python int
    return {z+1:int(sum(1<<t for t in np.flatnonzero(Par[z]&f))) for z in range(511)}
def enum(f,cap,maxd=6):
    km=kern_masks(f); cheap={z:m for z,m in km.items() if bin(m).count('1')<=cap}
    # subspaces as frozenset of nonzero ints
    level={frozenset([z]):m for z,m in cheap.items()}
    res={}  # frozenset -> kernel mask
    allres=dict(level)
    cur=level
    for d in range(2,maxd+1):
        nxt={}
        for Zs,m in cur.items():
            for z,mz in cheap.items():
                if z in Zs: continue
                nz=Zs|frozenset(a^z for a in Zs)|frozenset([z])
                if nz in nxt: continue
                if not all(a in cheap for a in nz): continue
                mm=m|mz
                if bin(mm).count('1')<=cap: nxt[nz]=mm
        allres.update(nxt); cur=nxt
        if not nxt: break
    return allres
if __name__=='__main__':
    ti=int(sys.argv[1]); cap=int(sys.argv[2])
    idx=np.flatnonzero(T==ti)
    import collections
    stat=collections.Counter()
    ex=None
    for k,i in enumerate(idx[:int(sys.argv[3])]):
        r=enum(F[i],cap)
        for Zs,m in r.items():
            d=len(Zs).bit_length()  # |Z|=2^d-1
            stat[(d,bin(m).count('1'))]+=1
    print(ti,sorted(stat.items()))
