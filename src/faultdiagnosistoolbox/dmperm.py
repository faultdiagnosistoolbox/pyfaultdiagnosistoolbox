import structuralanalysis as sa
import scipy.sparse as sp
import numpy as np

def CSCDict(A):
    return {'nzmax':A.nnz,
     'm': A.shape[0],
     'n': A.shape[1],
     'p': A.indptr.astype(np.int64),
     'i': A.indices.astype(np.int64),
     'x': A.data.astype(np.float64),
     'nz': -1}  

def dmperm(A):
    return sa.dmperm_internal(CSCDict(A))

def safedmperm(A):
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
    row = np.array([],dtype=np.int64)
    col = np.array([],dtype=np.int64)

    def __init__(self,r,c):
        self.row = r
        self.col = c;

def MSO(X):
    if sp.issparse(X):
        return sa.findmso_internal(CSCDict(X))
    else:
        return sa.findmso_internal(CSCDict(sp.csc_matrix(X)))

def GetDMParts(X):
    if sp.issparse(X):
        dm = dmperm(X)
    else:
        dm = dmperm(sp.csc_matrix(X))
        
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

        res.M0eqs = np.array([],dtype=np.int64)
        res.M0vars = np.array([],dtype=np.int64)

        for hc in res.M0:
            res.M0eqs = np.concatenate((res.M0eqs,hc.row))
            res.M0vars = np.concatenate((res.M0vars,hc.col))
    return res

def PSODecomposition(X):
    if not IsPSO(X):
        print "PSO Decomposition for PSO structures only, exiting..."
        return;

    n = X.shape[0]
    m = X.shape[1]

    delrows=np.arange(0,n)
    eqclass = []
    trivclass = []
    Mi = np.array([]);

    while len(delrows)>0:
        temprow = np.array([x for x in np.arange(0,n) if not x == delrows[0]])
        dm = GetDMParts(X[temprow,:]);

        if len(dm.M0vars)>0:
            ec_row = np.sort(np.append(temprow[dm.M0eqs], delrows[0]))
            ec = EqBlock(ec_row,np.sort(dm.M0vars))
            eqclass.append(ec)
            Mi = np.concatenate((Mi,ec.col)).astype(np.int64)
            delrows = [x for x in delrows if x not in ec.row]
        else:
            trivclass.append(delrows[0])
            delrows = delrows[1:len(delrows)]
    X0 = np.sort([x for x in np.arange(0,m) if not x in Mi])
    if len(X0)==0:
        X0=[]
    res = {}
    res['eqclass'] = eqclass
    res['trivclass'] = trivclass
    res['X0'] = X0

    p = [];
    q = [];
    
    for ec in eqclass:
        p.append(ec.row)
        q.append(ec.col)
    p.append(trivclass)
    q.append(X0)
    
    res['p'] = p
    res['q'] = q
    return res

def IsPSO( X, *args ):
  if len(args)>0:
    eq=args[0]
  else:
    eq=np.arange(0,X.shape[0])
      
  dm = GetDMParts(X[eq,:])    
  return (len(dm.Mm.row)==0) and (len(dm.M0)==0)
