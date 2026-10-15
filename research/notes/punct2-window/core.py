"""PUNCT core: units as integer line vectors, exact dual-basis decoders, exact F2 ranks.
Vectors are in doubled coordinates (E8 roots norm 8; inner products /4 -> G)."""
import itertools, sys
from fractions import Fraction as F
import numpy as np

def e8_roots2():
    R=[]
    for i,j in itertools.combinations(range(8),2):
        for si in(2,-2):
            for sj in(2,-2):
                r=[0]*8; r[i]=si; r[j]=sj; R.append(tuple(r))
    for s in itertools.product((1,-1),repeat=8):
        if s.count(-1)%2==0: R.append(tuple(s))
    return R
def lines_of(R):
    seen=set();L=[]
    for r in R:
        n=tuple(-x for x in r)
        if n in seen or r in seen: continue
        seen.add(r); L.append(r)
    return L
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def dn2(n):
    R=[]
    for i,j in itertools.combinations(range(n),2):
        for si in(2,-2):
            for sj in(2,-2):
                r=[0]*n; r[i]=si; r[j]=sj; R.append(tuple(r))
    return R
def unit(name):
    E8=lines_of(e8_roots2())
    if name=='E8': L=E8
    elif name=='E7':
        r0=E8[0]; L=[r for r in E8 if dot(r,r0)==0]
    elif name=='E6':
        r0=E8[0]; r1=next(r for r in E8 if r!=r0 and dot(r,r0)!=0)
        L=[r for r in E8 if dot(r,r0)==0 and dot(r,r1)==0]
    elif name[0]=='D':
        n=int(name[1:]); L=lines_of(dn2(n))
    return L
def rank_f2(rows):
    basis=[]
    for r in rows:
        for b in basis:
            r=min(r,r^b)
        if r: basis.append(r)
    return len(basis)
class Unit:
    def __init__(self,name):
        self.name=name; self.L=unit(name); self.v=len(self.L)
        L=self.L; v=self.v
        self.dim=len(L[0])
        # rank of span
        self.T=self.rank_span()
        self.Gram=[[dot(a,b)//4 if (dot(a,b)%4==0) else None for b in L] for a in L]
        assert all(g is not None for r in self.Gram for g in r)
        self.Brow=[sum(1<<s for s in range(v) if self.Gram[t][s]==0 or s==t) for t in range(v)]
        self.A=np.array(L,dtype=np.int64)
        self.full=(1<<v)-1
    def rank_span(self):
        rows=[list(map(F,r)) for r in self.L]; rk=0; n=len(rows[0])
        for c in range(n):
            p=next((i for i in range(rk,len(rows)) if rows[i][c]!=0),None)
            if p is None: continue
            rows[rk],rows[p]=rows[p],rows[rk]
            for i in range(len(rows)):
                if i!=rk and rows[i][c]!=0:
                    f=rows[i][c]/rows[rk][c]; rows[i]=[a-f*b for a,b in zip(rows[i],rows[rk])]
            rk+=1
        return rk
    def mask_of(self,dots):
        m=0
        for t,d in enumerate(dots):
            if d!=0: m|=1<<t
        return m
    def term_masks(self,W):
        """W: list of T span vectors (integer). Return (Rmasks, Amasks, exact_ok) for dual-basis decoder."""
        T=self.T; Gm=[[dot(a,b) for b in W] for a in W]
        # inverse via Fractions
        M=[[F(x) for x in row]+[F(int(i==j)) for j in range(T)] for i,row in enumerate(Gm)]
        for c in range(T):
            p=next((r for r in range(c,T) if M[r][c]!=0),None)
            if p is None: return None
            M[c],M[p]=M[p],M[c]; pv=M[c][c]; M[c]=[x/pv for x in M[c]]
            for r in range(T):
                if r!=c and M[r][c]!=0:
                    f=M[r][c]; M[r]=[x-f*y for x,y in zip(M[r],M[c])]
        Gi=[row[T:] for row in M]
        # dual vector w'_j = sum_i w_i Gi[i][j]; dots with lines
        D=[[dot(l,w) for l in self.L] for w in W]   # T x v
        Rm=[];Am=[]
        for j in range(T):
            Rm.append(self.mask_of(D[j]))
            dd=[sum(D[i][t]*Gi[i][j] for i in range(T)) for t in range(self.v)]
            Am.append(self.mask_of(dd))
        return Rm,Am,Gi,D
    def cost(self,Rm,Am,live):
        tot=0; rs=[]
        for R,A in zip(Rm,Am):
            R&=live; A&=live
            rows=[self.Brow[t]&A for t in range(self.v) if R>>t&1]
            r=rank_f2(rows); rs.append(r); tot+=r
        return tot,rs
    def window(self,live):
        return rank_f2([self.Brow[t]&live for t in range(self.v) if live>>t&1])
def verify_exact(U,W,live):
    """exact check sum_j <r_t,w_j><r_s,w'_j> = <r_t,r_s> on live lines (Fractions)."""
    out=U.term_masks(W); Gi=out[2]; D=out[3]; T=U.T
    for t in range(U.v):
        if not live>>t&1: continue
        for s in range(U.v):
            if not live>>s&1: continue
            tot=sum(F(D[j][t])*sum(D[i][s]*Gi[i][j] for i in range(T)) for j in range(T))
            if tot!=U.Gram[t][s]*4: return False
    return True
