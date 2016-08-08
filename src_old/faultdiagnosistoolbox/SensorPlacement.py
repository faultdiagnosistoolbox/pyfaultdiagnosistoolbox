import numpy as np
import dmperm
from MHS import MHS

def ParentBlocks(X,b):
        return np.where(X[0:b,b])[0]

def FbDetectability(Xb,fb,dm):
    v = np.array([],dtype=np.int64)
    pb = np.concatenate((ParentBlocks(Xb,fb),[fb]))
    for jj in pb:
        v = np.concatenate((v,dm.M0[jj].col))
    return v

def SensorPlacementDetectability(model,**options):
    if options.has_key('fdet'):
        fdet = options['fdet']
        if isinstance(fdet[0],str):
            fdet = map(lambda fi: model.f.index(fi), fdet)
    else:
        fdet = np.arange(0,len(model.f))
        
    if model.IsUnderdetermined():
        print "Sorry, sensor placement only works for models with no underdetermined parts"
        return
    
    dm = dmperm.GetDMParts(model.X)
    ef = map(lambda fi: np.where(model.F[:,fi])[0][0], fdet)
    nondetIdx = [eIdx for eIdx in np.arange(0,len(ef)) if ef[eIdx] in dm.M0eqs]
    if len(nondetIdx)>0:
        ds = DetectabilitySets( model.X, model.F[:,fdet[nondetIdx]], model.P )
        sensSets = MHS( ds )
    else:
        sensSets = []

    # Compute list of variable names
    s = map(lambda s: list(np.array(model.x)[s]), sensSets)
    
    return (s,sensSets)

def NewSensorEqs(s,nx,nf,Pfault):
    ns = len(s)
    Xs = np.zeros((ns,nx),np.int64)
    Fs = np.zeros((ns,nf),np.int64)
    fs = []
    for sIdx, si in enumerate(s):
        Xs[sIdx,si] = 1
        if si in Pfault:
            nf = nf+1
            Fs = np.hstack((Fs,np.zeros((ns,1),dtype=np.int64))) # Add column for new fault
            Fs[sIdx,-1] = 1 # Add fault
            fs.append(si)
        else:
            Fs[sIdx,:] = np.zeros((1,nf))
    return (Xs, Fs, fs)

def IsolabilitySubProblem(X,F,P,fi):
    ef = np.where(F[:,fi])[0][0]
    n = X.shape[0]
    nf = F.shape[1]

    # Misol = M\{ef}
    Xisol = np.delete(X,ef,axis=0)
    Fisol = np.delete(F,ef,axis=0)

    # Extract just determined part of Xisol
    dm = dmperm.GetDMParts(Xisol)
    X0 = Xisol[dm.M0eqs,:][:,dm.M0vars]
    
    # Find out which faults are included in X0
    feq = np.zeros(Fisol.shape[1],dtype=np.int64)
    for fj in np.arange(0,Fisol.shape[1]):
        e = np.where(Fisol[:,fj])[0]
        if len(e)==0:
            e=-1
        feq[fj] = e

    nondet = np.where(map(lambda f: feq[f] in dm.M0eqs, np.arange(0,Fisol.shape[1])))[0]
    F0 = Fisol[dm.M0eqs,:][:,nondet]

    # Translate P to P0
    P0 = [np.where(dm.M0vars==v)[0][0] for v in P if v in dm.M0vars]

    # Compute detectability sets
    detSets = DetectabilitySets( X0, F0, P0 )

    # Translate back to original variable indices
    detSets = map(lambda ds: dm.M0vars[ds], detSets)
    
    return detSets

def BlockAndFaultOrder(X,F,dm):
    # 1. Construct block adjecency matrix
    #    Xb(i,j) = 1 => bi>bj
    n = len(dm.M0)
    Xb = np.zeros((n,n),dtype=np.int64)
    # 1.1 Determine connected blocks
    for rr in np.arange(0,n):
        for cc in np.arange(0,n):
            Xb[rr,cc] = np.any(X[dm.M0[rr].row,:][:,dm.M0[cc].col])
        
    # 1.2 Traverse block adjacency matrix to determine indirect relationships
    for bb in np.arange(1,n):
        ba = ParentBlocks(Xb,bb)
        iba = np.array([],dtype=np.int64)
        for ll in ba:
            iba = np.concatenate((iba,ParentBlocks(Xb,ll)))
        if len(iba)>0:
            Xb[np.unique(iba),bb]=1
        
    # 2. Construct fault classes and determine maximal elements 
    # 2.1 Determine e_f for each fault (must be 1 equation for each fault) 
    ef = np.zeros(F.shape[1], dtype=np.int64)
    for fi in np.arange(0,F.shape[1]):
        ef[fi] = np.where(F[:,fi])[0]

    # 2.2 Determine block membership for each fault
    efrep = np.unique(ef)
    efb = np.zeros(len(efrep), dtype=np.int64)
    for fi,efi in enumerate(efrep):
        for bjIdx,bj in enumerate(dm.M0):
            if efi in bj.row:
                efb[fi] = bjIdx
            
    # 2.3 Determine blocks corresponding to maximal elements in the
    #     fault class partial order
    maxFaultClasses = np.zeros(len(efb),dtype=np.int64)
    for jj in np.arange(0,len(maxFaultClasses)):
        maxFaultClasses[jj] = len([b for b in ParentBlocks(Xb,efb[jj]) if b in efb])==0

    bFm = efb[maxFaultClasses>0]
    return {'blockorder':Xb, 'bFm':bFm}

def SensPlaceM0(X,F):
    dm0 = dmperm.GetDMParts(X)
    bfOrder = BlockAndFaultOrder( X, F, dm0 )
    return map(lambda b: FbDetectability(bfOrder['blockorder'],b,dm0),bfOrder['bFm'])

def DetectabilitySets(X,F,P):
    dm = dmperm.GetDMParts(X)
    detSets = SensPlaceM0(X[dm.M0eqs,:][:,dm.M0vars],F[dm.M0eqs,:])
    return map(lambda d: [dm.M0vars[x] for x in d if x in P], detSets)

def SensorPlacementIsolability(model):
    Pfault = model.Pfault
    fault = np.arange(0,len(model.x))
    _,spDet = SensorPlacementDetectability(model)
    sIdx = []
    if len(spDet)>0:
        for s in spDet:
            Xs,Fs,fs = NewSensorEqs(s,model.X.shape[1],model.F.shape[1],Pfault)
            X = model.X
            F = np.hstack((model.F,np.zeros((model.F.shape[0],len(fs)),dtype=np.int64)))
            X = np.vstack((X,Xs))
            F = np.vstack((F,Fs))
            P = model.P

            detSets = []
            for fi in np.arange(0,F.shape[1]):
                isolDetSets = IsolabilitySubProblem(X,F,P,fi)
                if len(isolDetSets)>0:
                    detSets = detSets + isolDetSets
            mhs_s = MHS( detSets )
            sIsolIdx = map(lambda m: np.sort(np.concatenate((m,s))), mhs_s)
            sIdx = sIdx + sIsolIdx
    else:
        X = model.X
        F = model.F
        P = model.P

        detSets = []
        for fi in np.arange(0,F.shape[1]):
            isolDetSets = IsolabilitySubProblem(X,F,P,fi)
            if len(isolDetSets)>0:
                detSets = detSets + isolDetSets
        sIdx = MHS( detSets )

    # Remove duplicates
    keepIdx = np.ones(len(sIdx),dtype=np.int64).astype(np.bool)
    for k in np.arange(0,len(sIdx)):
        for l in np.arange(k+1,len(sIdx)):
            if np.array_equal(sIdx[k],sIdx[l]):
                keepIdx[l]=False
    sIdx = [sIdx[k] for k in np.arange(0,len(sIdx)) if keepIdx[k]]

    # Compute list of variable names
    s = map(lambda s: list(np.array(model.x)[s]), sIdx)
    
    return (s,sIdx)
