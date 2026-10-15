"""Bit model: B = U Omega U^T over F2 (U: v x h label coords, Omega = Bsub invertible symmetric h x h).
Decoder = Omega-dual bases (f_j,g_j): term j has reader R_j = {t: Omega(u_t,f_j)=1}, source A_j = {s: Omega(u_s,g_j)=1},
charge rank_Q M[R_j & K, A_j & K] = rank(S_R S_A), S_X = sum_{t in X} r_t r_t^T."""
import itertools, random, math, json, sys, time
import numpy as np
from core import *
def popcount(x): return bin(x).count('1')
class BitUnit:
    def __init__(self,name):
        self.U=Unit(name); U=self.U; v=U.v; self.v=v; self.name=name
        # pick basis lines: greedy independent rows of Brow giving invertible principal submatrix
        rows=U.Brow
        h=rank_f2(rows); self.h=h
        # find principal nonsingular h x h : greedy by adding indices keeping principal submatrix full rank
        T=[]; cur=[]
        for t in range(v):
            if rank_f2(cur+[rows[t]])>len(cur):
                cur.append(rows[t]); T.append(t)
            if len(T)==h: break
        assert len(T)==h
        self.T=T
        Bs=[[ (rows[a]>>b)&1 for b in T] for a in T]; self.Bs=Bs
        # coords c_t: solve B[t] restricted rows in T: B[t,s] = sum_i c_ti * (Bs row i col s) ... for s in T: row_t|T = c_t Bs
        # invert Bs over F2
        inv=self.inv2(Bs)
        self.c=[]
        for t in range(v):
            rt=[(rows[t]>>b)&1 for b in T]
            c=[sum(rt[k]*inv[k][i] for k in range(h))%2 for i in range(h)]
            self.c.append(c)
        # verify
        for t in range(v):
            for s in range(0,v,7):
                val=sum(self.c[t][i]*Bs[i][j]*self.c[s][j] for i in range(h) for j in range(h))%2
                assert val==(rows[t]>>s)&1
        self.Om=np.array(Bs,dtype=np.int64)
        self.Cm=np.array(self.c,dtype=np.int64)               # v x h
        # all vectors f in F2^h (nonzero): reader mask
        vecs=np.array([[ (x>>i)&1 for i in range(h)] for x in range(1,1<<h)],dtype=np.int64)  # n x h
        self.vecs=vecs
        self.OmF=(vecs@self.Om)%2                              # n x h : row f -> Omega f (Omega symmetric)
        self.Rbits=((self.Cm@self.OmF.T)%2).T.astype(bool)     # n x v : Omega(u_t,f)=1
        self.OmM=((vecs@self.Om@vecs.T)%2).astype(bool)        # Omega(f,g)
        A=np.array(U.L,dtype=np.int64)
        self.O=np.einsum('ta,tb->tab',A,A)                     # v x dim x dim
        self.dimn=A.shape[1]
    @staticmethod
    def inv2(M):
        n=len(M); A=[row[:]+[int(i==j) for j in range(n)] for i,row in enumerate(M)]
        for c in range(n):
            p=next(r for r in range(c,n) if A[r][c]); A[c],A[p]=A[p],A[c]
            for r in range(n):
                if r!=c and A[r][c]: A[r]=[(x^y) for x,y in zip(A[r],A[c])]
        return [row[n:] for row in A]
    def cost_np(self,MR,MA,live):
        """MR,MA: (m,v) bool arrays of reader/source membership; live: (v,) bool. Returns total, list."""
        lv=live.astype(np.int64)
        SR=np.einsum('jt,tab->jab',MR*lv,self.O); SA=np.einsum('jt,tab->jab',MA*lv,self.O)
        P=SR@SA
        rk=np.linalg.matrix_rank(P.astype(float))
        return int(rk.sum()),rk
    def exact_cost(self,MR,MA,live):
        from fractions import Fraction as F
        tot=0; rs=[]
        for j in range(MR.shape[0]):
            R=[t for t in range(self.v) if MR[j][t] and live[t]]; A=[t for t in range(self.v) if MA[j][t] and live[t]]
            if not R or not A: rs.append(0); continue
            M=[[F(self.U.Gram[t][s]) for s in A] for t in R]
            rk=0
            rows=M; ncol=len(A)
            for c in range(ncol):
                p=next((i for i in range(rk,len(rows)) if rows[i][c]!=0),None)
                if p is None: continue
                rows[rk],rows[p]=rows[p],rows[rk]
                for i in range(rk+1,len(rows)):
                    if rows[i][c]!=0:
                        f=rows[i][c]/rows[rk][c]; rows[i]=[a-f*b for a,b in zip(rows[i],rows[rk])]
                rk+=1
            rs.append(rk); tot+=rk
        return tot,rs
    def window(self,live):
        return rank_f2([self.U.Brow[t]&sum(1<<s for s in range(self.v) if live[s]) for t in range(self.v) if live[t]])
def optS_bit(BU,MR,MA,lam,iters,rnd,live=None,T0=0.6,minv=20):
    v=BU.v
    live=np.ones(v,dtype=bool) if live is None else live.copy()
    c,_=BU.cost_np(MR,MA,live); g=c-lam*live.sum(); best=(g,live.copy(),c)
    for it in range(iters):
        temp=T0*(1-it/iters)+0.02
        t=rnd.randrange(v); live[t]=~live[t]
        if live.sum()<minv: live[t]=~live[t]; continue
        nc,_=BU.cost_np(MR,MA,live); ng=nc-lam*live.sum()
        if ng<=g or rnd.random()<math.exp(-(ng-g)/temp):
            g=ng;c=nc
            if g<best[0]-1e-9: best=(g,live.copy(),c)
        else: live[t]=~live[t]
    return best
