import numpy as np, pickle, collections, time
from setup_perm import *
tr=transporters()
a0,out=pickle.load(open('minopts_root.pkl','rb'))
print('rep flat idx',a0)
# which line is alpha for rep flat? the unique z with weight 1 kernel
ks=out[2]
# root line index: kernel of z0 = {alpha}; find from flats: complement mask of flat contains alpha...
F=np.load('flats.npy'); T=np.load('ftypes.npy')
f0=F[a0]
# alpha = the line in C_f0 whose orthogonal complement is f0's flat
Gram=np.array(U.Gram)
alpha=[t for t in range(120) if set(np.flatnonzero(Gram[t]==0))==set(np.flatnonzero(~f0))-{t} or False]
cand=[t for t in range(120) if (Gram[t]==0).tolist()==(~f0).tolist()]
print('alpha cand',cand)
alpha=cand[0]
print('alpha',alpha)
pickle.dump(alpha,open('alpha.pkl','wb'))
