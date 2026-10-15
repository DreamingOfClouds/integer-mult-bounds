import numpy as np, itertools
from ortools.sat.python import cp_model
L=np.load('flats.npy')[:0]
import sys; sys.path.insert(0,'.')
from core import *
U=Unit('E8'); G=np.array(U.Gram)
orth=(G==0)
v=120
frames=[]
def ext(cur,cands):
    if len(cur)==8: frames.append(tuple(cur)); return
    for k,c in enumerate(cands):
        ext(cur+[c],[d for d in cands[k+1:] if orth[c,d]])
ext([],list(range(v)))
print('frames',len(frames))
# also 7-sets maximal? min transversal of all 8-frames
m=cp_model.CpModel(); s=[m.NewBoolVar('') for _ in range(v)]
for fr in frames: m.AddBoolOr([s[t] for t in fr])
m.Minimize(sum(s))
sv=cp_model.CpSolver(); sv.parameters.num_workers=1; sv.parameters.max_time_in_seconds=900
r=sv.Solve(m); print(sv.StatusName(r),sv.ObjectiveValue(),sv.BestObjectiveBound())
