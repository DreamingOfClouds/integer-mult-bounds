import numpy as np, sys, pickle, collections, time
F=np.load('flats.npy'); T=np.load('ftypes.npy'); Cm=np.load('Cm.npy')
Zs=np.array([[(z>>i)&1 for i in range(9)] for z in range(1,512)])
Par=((Zs@Cm.T)%2).astype(bool)
pc=lambda x: bin(x).count('1')
def enum_f(f,cap,maxd=5):
    km={z+1:sum(1<<int(t) for t in np.flatnonzero(Par[z]&f)) for z in range(511)}
    cheap={z:m for z,m in km.items() if pc(m)<=cap}
    cl=sorted(cheap)
    res={}  # frozenset subspace -> mask
    cur={}
    for i,z1 in enumerate(cl):
        for z2 in cl[i+1:]:
            if z2^z1 in cheap or True:
                z3=z1^z2
                if z3 not in cheap: continue
                m=cheap[z1]|cheap[z2]
                if pc(m)<=cap:
                    cur[frozenset((z1,z2,z3))]=m
    res.update(cur)
    d=2
    while cur and d<maxd:
        nxt={}
        for S_,m in cur.items():
            for z in cl:
                if z in S_: continue
                new=S_|frozenset(a^z for a in S_)|frozenset([z])
                if new in nxt: continue
                if not all(a in cheap for a in new): continue
                mm=m|cheap[z]
                if pc(mm)<=cap: nxt[new]=mm
        res.update(nxt); cur=nxt; d+=1
    return {s:m for s,m in res.items()}, cheap
if __name__=='__main__':
    ti=int(sys.argv[1]); cap=int(sys.argv[2]); n=int(sys.argv[3])
    idx=np.flatnonzero(T==ti)[:n]
    for i in idx:
        t0=time.time()
        res,cheap=enum_f(F[i],cap)
        st=collections.Counter((len(s).bit_length(),pc(m)) for s,m in res.items())
        print(ti,i,'cheap',len(cheap),sorted(st.items()),round(time.time()-t0,1),flush=True)
