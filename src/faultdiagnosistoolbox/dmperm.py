import faultdiagnosistoolbox.structuralanalysis as sa
import scipy.sparse as sp
import numpy as np
import copy

def CSCDict(A):
    return {'nzmax':A.nnz,
     'm': A.shape[0],
     'n': A.shape[1],
     'p': A.indptr.astype(np.int64),
     'i': A.indices.astype(np.int64),
     'x': A.data.astype(np.float64),
     'nz': -1}  

def dmperm(A):
    if sp.issparse(A):
        return sa.dmperm_internal(CSCDict(A))
    else:
        return sa.dmperm_internal(CSCDict(sp.csc_matrix(A)))
def srank(A):
    _,_,_,_,rr,_ = dmperm(A)
    return rr[3]

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
        p,q,r,s,_,_ = dmperm(X)
    else:
        p,q,r,s,_,_ = dmperm(sp.csc_matrix(X))
        
    m = X.shape[0]
    n = X.shape[1]
    
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
        print("PSO Decomposition for PSO structures only, exiting...")
        return;

    n = X.shape[0]
    m = X.shape[1]

    delrows=np.arange(0,n)
    eqclass = []
    trivclass = np.array([],dtype=np.int64)
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
            trivclass = np.concatenate((trivclass,[delrows[0]]))
            delrows = delrows[1:len(delrows)]
    X0 = np.sort([x for x in np.arange(0,m) if not x in Mi])
    if len(X0)==0:
        X0=[]
    res = {}
    res['eqclass'] = eqclass
    res['trivclass'] = trivclass
    res['X0'] = X0

    p = np.array([],dtype=np.int64)
    q = np.array([],dtype=np.int64)
    
    for ec in eqclass:
        p = np.concatenate((p,ec.row))
        q = np.concatenate((q,ec.col))
    p = np.concatenate((p,trivclass))
    q = np.concatenate((q,X0))
    
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

#def IsHighIndex(X, **opts):
#    if 'eq' in opts:
#        eq = opts['eq']
#    else:
#        eq = np.arange(0,X.shape[0])
#    X1 = X[eq,:]
#    x1 = np.argwhere(X1==3)
#    row_d = x1[:,0]
#    row_a = np.array([r for r in np.arange(0,X1.shape[0]) if not r in row_d])
#    col_d1 = x1[:,1]
#    col_1 = np.argwhere(X1[row_d,:]==2)[:,1]
#    col_2 = np.array([c for c in range(0,X1.shape[1]) if not (c in col_d1 or c in col_1)])
#    col_2 = [c for c in col_2 if np.any(X1[:,c],axis=0)]
#    Xhod = X1[row_a,:][:,np.concatenate((col_2, col_d1))]
#    
#    return srank(Xhod)<Xhod.shape[1]

def IsHighIndex(X, eq=[]):
    if len(eq)==0:
        eq = np.arange(0,self.X.shape[0])
    X1 = X[eq,:]
    
    col_d1 = np.any(X1==3,axis=0)
    row_a = np.all(X1<=1,axis=1)
    row_d = np.any(X1>1,axis=1)
    col_2 = np.all(X1[row_d,:]==0,axis=0)
    Xhod = X1[row_a,:][:,np.logical_or(col_d1,col_2)]
    nz = sum(np.all(Xhod==0,axis=0))
    return srank(Xhod)<Xhod.shape[1]-nz

def IsLowIndex(X, eq=[]):
    if len(eq)==0:
        eq = np.arange(0,self.X.shape[0])
    return not IsHighIndex(X,eq)

def Mplus( X, causality='mixed' ):
    def Gp(G):
        dm = GetDMParts(G[2])
        if len(dm.Mp.row)==0 or len(dm.Mp.col)==0:
            return (np.array([]),np.array([]),np.array([[]]))
        G1 = copy.deepcopy(G)
        return (G1[0][dm.Mp.row], G1[1][dm.Mp.col], G1[2][dm.Mp.row,:][:,dm.Mp.col])

    def Gm(G):
        dm = GetDMParts(G[2])
        if len(dm.Mm.row)==0 or len(dm.Mm.col)==0:
            return (np.array([]),np.array([]),np.array([[]]))
        G1 = copy.deepcopy(G)
        return (G1[0][dm.Mm.row], G1[1][dm.Mm.col], G1[2][dm.Mm.row,:][:,dm.Mm.col])

    def G0(G):
        dm = GetDMParts(G[2])
        if len(dm.M0eqs)==0 or len(dm.M0vars)==0:
            return (np.array([]),np.array([]),np.array([[]]))
        G1 = copy.deepcopy(G)
        return (G1[0][dm.M0eqs], G1[1][dm.M0vars], G1[2][dm.M0eqs,:][:,dm.M0vars])

    def Gadd(G1,G2):
        c1,x1,A1 = copy.deepcopy(G1)
        c2,x2,A2 = copy.deepcopy(G2)

        c1 = np.unique(np.concatenate((c1,c2)))
        x1 = np.unique(np.concatenate((x1,x2)))
        A1[A2==1]=1; A1[A2==2]=2
        return (c1,x1,A1)
    
    def CGX(G,X):
        if len(X)==0:
            return []
        xIdx = [np.where(G[1]==ii)[0][0] for ii in X]
        idx = np.unique([e[0] for e in np.argwhere(G[2][:,xIdx]>0)])
        if len(idx)==0:
            return []
        
        return G[0][idx]

    def CGE(G,E):
        return G[0][np.where(np.any(Ei,axis=1))[0]]

    def GsubC(G,C):
        c,x,A = copy.deepcopy(G)
        A = np.delete(A,C,axis=0)
        c = np.array([ci for ci in c if not ci in C])
        return (c,x,A)

    def GsubX(G,X):
        c,x,A = copy.deepcopy(G)        
        A = np.delete(A,X,axis=1)
        x = np.array([xi for xi in x if not xi in X])
        return (c,x,A)
    
    Xc = np.array(X.copy())
    # Represent graph as a tuple G=(constraints,variables, adjacency matrix)
    G = (np.arange(0,Xc.shape[0], dtype=np.int64), np.arange(0,Xc.shape[1], dtype=np.int64), Xc)

    if causality is 'mixed':
        return Gp(G)[0]
    
    elif causality is 'int':
        while True:
            # G := G+            
            G = Gp(G)            
                
            # G1 = G - Ed
            G1 = copy.deepcopy(G)
            G1[2][G1[2]==3]=0 # Zero the derivative edges

            # G := G - C(G,X(G1-))
            dm = GetDMParts(G1[2])

            if len(dm.Mm.col)>0 and len(dm.Mm.row)>0:
                G1m = Gm(G1)
                G=GsubC(G,CGX(G,G1m[1]))                
            else:
                break
        return G[0]
    elif causality is 'der':
        Xc = [];
        while True:
            # Gnc := G - Xc
            Gnc = GsubX(G,Xc)
            
            # Gni := Gnc - C(Gnc,Ei)
            Ei = Gnc[2].copy()
            Ei[Ei!=2]=0
            Gni = GsubC(Gnc,CGE(Gnc,Ei))
            
            # Xc:= Xc u X(Gni+ u Gni0)        
            Gnip = Gp(Gni)
            Gni0 = G0(Gni)
            X1 = np.unique(np.concatenate((Gnip[1],Gni0[1]))).astype(np.int64) # X1 = X(Gni+ u Gni0)
            Xc = np.unique(np.concatenate((Xc,X1))).astype(np.int64)

            if len(X1)==0:
                break
            
        # G := (G-C(G,X\Xc))+
        X1 = np.array([xi for xi in G[1] if not xi in Xc]) # X1 = X\Xc
        Gres = Gp(GsubC(G,CGX(G,X1)))
        return Gres[0]
    else:
        return np.array([])

def MixedCausalityMatching(Gin):
    def FindAdmissibleIntEdge(G):
        X = G[2]
        dm=GetDMParts(X)
        for hallComponent in dm.M0:
            Ei = np.argwhere(X[hallComponent.row,:][:,hallComponent.col]==2)
            if len(Ei)>0:
                Ei=Ei[0] # Get the first one
                return np.array([G[0][hallComponent.row[Ei[0]]], G[1][hallComponent.col[Ei[1]]]])
        return []

    G = [Gin[0].copy(), Gin[1].copy(), Gin[2].copy()]
    Gamma = []
    e = FindAdmissibleIntEdge(G)
    while len(e)>0:
        if len(Gamma)==0:
            Gamma = [e]
        else:
            Gamma = np.concatenate((Gamma,[e]),axis=0)
            
        # G := G - C({e}) - X({e})
        rowIdx = np.argwhere(G[0]==e[0])[0,0]
        colIdx = np.argwhere(G[1]==e[1])[0,0]
        G[0] = np.delete(G[0],rowIdx)
        G[1] = np.delete(G[1],colIdx)
        G[2] = np.delete(G[2],rowIdx, axis=0)
        G[2] = np.delete(G[2],colIdx, axis=1)
        e = FindAdmissibleIntEdge(G)
        
    dm = GetDMParts(G[2])
    Gamma_p = np.stack((G[0][dm.rowp][:],G[1][dm.colp][:]), axis=1)
    return np.flipud(np.concatenate((Gamma,Gamma_p),axis=0))
