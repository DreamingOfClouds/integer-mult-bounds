import numpy as np, sys, time
from ortools.sat.python import cp_model
F=np.load('flats.npy'); T=np.load('ftypes.npy'); Cm=np.load('Cm.npy')
Zs=np.array([[(z>>i)&1 for i in range(9)] for z in range(1,512)])
Par=((Zs@Cm.T)%2).astype(bool)
p=int(sys.argv[1]); tl=float(sys.argv[2]); ftype=int(sys.argv[3]) if len(sys.argv)>3 else 8
cap=p
idx=np.flatnonzero(T==ftype)
m=cp_model.CpModel()
s=[m.NewBoolVar(f's{t}') for t in range(120)]
m.Add(sum(s)<=p)
pick=[];b2=[];b3=[]
nvars=0
for a,i in enumerate(idx):
    f=F[i]
    ys=[]
    for z in range(511):
        K=np.flatnonzero(Par[z]&f)
        if len(K)>cap or len(K)==0: continue
        y=m.NewBoolVar('')
        for t in K: m.AddImplication(y,s[int(t)])
        ys.append(y)
    B2=m.NewBoolVar(''); B3=m.NewBoolVar('')
    m.AddImplication(B3,B2)
    m.Add(sum(ys)>=3).OnlyEnforceIf(B2)
    m.Add(sum(ys)>=7).OnlyEnforceIf(B3)
    b2.append(B2);b3.append(B3)
m.Add(sum(b2)<=8)   # top-8 selection: at most 8 flats counted
m.Maximize(sum(b2)+sum(b3))
sv=cp_model.CpSolver(); sv.parameters.max_time_in_seconds=tl; sv.parameters.num_workers=2
t0=time.time(); r=sv.Solve(m)
print('p',p,'status',sv.StatusName(r),'obj',sv.ObjectiveValue(),'bound',sv.BestObjectiveBound(),'time',round(time.time()-t0,1),flush=True)
