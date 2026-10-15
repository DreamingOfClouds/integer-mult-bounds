import sys, random, math, json, time, os
from core import *
def vectors_of_norm(U,n2):
    # span vectors of doubled-norm n2 in the unit's ambient lattice (E8-lattice-like, doubled coords): same parity coords, sum%4==0
    dim=U.dim; out=[]
    def rec(i,cur,rem,par):
        if i==dim:
            if rem==0 and sum(cur)%4==0 and any(cur): out.append(tuple(cur))
            return
        m=int(math.isqrt(rem))
        for c in range(-m,m+1):
            if par is not None and (c%2)!=par: continue
            if (dim-i-1)*0>rem-c*c: continue
            rec(i+1,cur+[c],rem-c*c,c%2 if par is None else par)
    rec(0,[],n2,None)
    # keep one of +-
    S=set();res=[]
    for v in out:
        if tuple(-x for x in v) in S: continue
        S.add(v); res.append(v)
    return res
def in_span(U,vs):
    return vs
def main():
    unit_name=sys.argv[1]; norms=[int(x) for x in sys.argv[2].split(',')]; seed=int(sys.argv[3]); iters=int(sys.argv[4]); start=sys.argv[5]
    tag=sys.argv[6] if len(sys.argv)>6 else 'run'
    rnd=random.Random(seed)
    U=Unit(unit_name); v=U.v; T=U.T
    pool=[]
    for n2 in norms: pool+=vectors_of_norm(U,n2)
    # restrict to span of the unit: for E7/E6/D: project check via rank of Gram: keep vectors whose dot-profile... use orthogonal complement test
    if unit_name!='E8':
        # span basis check: vector in span iff adding to lines does not raise rank; cheap: compute via projection on integer solves skipped -> filter by lines combos: use fractions rank
        base=[list(map(F,l)) for l in U.L]
        def rk(rows):
            rows=[r[:] for r in rows]; rkk=0
            for c in range(len(rows[0])):
                p=next((i for i in range(rkk,len(rows)) if rows[i][c]!=0),None)
                if p is None: continue
                rows[rkk],rows[p]=rows[p],rows[rkk]
                for i in range(len(rows)):
                    if i!=rkk and rows[i][c]!=0:
                        f=rows[i][c]/rows[rkk][c]; rows[i]=[a-f*b for a,b in zip(rows[i],rows[rkk])]
                rkk+=1
            return rkk
        # span basis: first T independent lines
        bas=[];
        for l in U.L:
            if rk([list(map(F,b)) for b in bas]+[list(map(F,l))])>len(bas): bas.append(l)
            if len(bas)==T: break
        bas=[list(map(F,b)) for b in bas]
        pool=[p for p in pool if rk(bas+[list(map(F,p))])==T]
    print('pool',len(pool),flush=True)
    # initial
    if start=='p104':
        d=json.load(open('p104.json')); basis=[U.L[i] for i in d['alpha']]
        live=U.full
        for i in d['F0']+d['F1']: live&=~(1<<i)
    elif start=='e8':
        d=json.load(open('p104.json')); basis=[U.L[i] for i in d['F0']]; live=U.full
    else:
        while True:
            basis=[rnd.choice(pool) for _ in range(T)]
            if U.term_masks(basis) is not None: break
        live=U.full
    def evaluate(basis,live):
        out=U.term_masks(basis)
        if out is None: return None
        Rm,Am=out[0],out[1]
        c,rs=U.cost(Rm,Am,live); vv=bin(live).count('1')
        return (c/vv if vv else 9, c, vv, (Rm,Am))
    cur=evaluate(basis,live); best=cur; bestst=(list(basis),live)
    Tt=float(os.environ.get('TT','0.02')); t0=time.time(); log=open(f'log_{tag}_{unit_name}_{seed}.txt','a')
    for it in range(iters):
        temp=Tt*(1-it/iters)+1e-4
        mv=rnd.random()
        if mv<0.75:
            t=rnd.randrange(v); nl=live^(1<<t)
            if bin(nl).count('1')<T+2: continue
            c,rs=U.cost(cur[3][0],cur[3][1],nl); vv=bin(nl).count('1')
            nf=c/vv
            if nf<=cur[0] or rnd.random()<math.exp(-(nf-cur[0])/temp):
                live=nl; cur=(nf,c,vv,cur[3])
        else:
            j=rnd.randrange(T); nb=list(basis); nb[j]=rnd.choice(pool)
            e=evaluate(nb,live)
            if e is None: continue
            if e[0]<=cur[0] or rnd.random()<math.exp(-(e[0]-cur[0])/temp):
                basis=nb; cur=e
        if cur[0]<best[0]-1e-12:
            best=cur; bestst=(list(basis),live)
            w=U.window(live)
            log.write(json.dumps({'it':it,'x':cur[0],'cost':cur[1],'v':cur[2],'window':w,'basis':[list(b) for b in basis],'live':hex(live)})+'\n'); log.flush()
        if it%2000==0: print(it,'cur %.4f (%d/%d) best %.4f (%d/%d) t=%.0f'%(cur[0],cur[1],cur[2],best[0],best[1],best[2],time.time()-t0),flush=True)
    print('BEST',best[:3])
if __name__=='__main__': main()
