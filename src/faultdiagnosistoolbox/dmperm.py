import structuralanalysis as sa
import scipy.sparse as sp
import numpy as np

def CSCDict(A):
    A=sp.csc_matrix(A)
    return {'nzmax':A.nnz,
     'm': A.shape[0],
     'n': A.shape[1],
     'p': A.indptr.astype(np.int64),
     'i': A.indices.astype(np.int64),
     'x': A.data.astype(np.float64),
     'nz': -1}  

def dmperm(A):
    return sa.dmperm_internal(CSCDict(A))

def sprank(A):
    return dmperm(A)['rr'][3]

class DMResult:
    def __init__(self):
        self.Mm = []
        self.M0 = []
        self.Mp = []
        self.rowp = []
        self.colp = []
        self.M0eqs = []
        self.M0vars = []
    
class EqBlock:
    row = []
    col = []

    def __init__(self,r,c):
        self.row = r
        self.col = c;

def GetDMParts(X):
    dm = dmperm(X)
    m = X.shape[0]
    n = X.shape[1]
    p = dm['p']
    q = dm['q']
    r = dm['r']
    s = dm['s']
    
    res = DMResult()
    if (p.size==0) or (q.size==0):
        if m>n:
            res.Mm = EqBlock([],[])
            res.Mp = EqBlock(np.arange(0,m),np.arange(0,n))
        else:
            res.Mp = EqBlock([],[])
            res.Mm = EqBlock(np.arange(0,m),np.arange(0,n))

        res.rowp=np.arange(1,m)
        res.colp=np.arange(1,n)
        res.M0=[]
    else:
        idx = 0
        if s[1]>r[1]: # Existance of M-
            res.Mm = EqBlock(np.sort(p[r[0]:r[1]]),np.sort(q[s[0]:s[1]]))
            idx = idx + 1
        else:
            res.Mm = EqBlock([],[])

        while (idx<r.size-1) and (r[idx+1]-r[idx])==(s[idx+1]-s[idx]): # M0 block exists
            res.M0.append(EqBlock(np.sort(p[r[idx]:r[idx+1]]),np.sort(q[s[idx]:s[idx+1]])))
            idx=idx+1

        if idx<r.size-1: # M+ exists
            res.Mp = EqBlock(np.sort(p[r[idx]:r[idx+1]]),np.sort(q[s[idx]:s[idx+1]]))
        else:
            res.Mp = EqBlock([],[])
        res.rowp = p
        res.colp = q

        res.M0eqs = []
        res.M0vars = []

        for hc in res.M0:
            res.M0eqs = np.concatenate((res.M0eqs,hc.row))
            res.M0vars = np.concatenate((res.M0vars,hc.col))
        res.M0eqs = np.sort(res.M0eqs).astype(np.int64) # Ugly!
        res.M0vars = np.sort(res.M0vars).astype(np.int64) # Ugly!
    return res
