"""Sensor placement algorithms."""
import numpy as np
# import faultdiagnosistoolbox.dmperm as dmperm
# from faultdiagnosistoolbox import MHS, GetDMParts
from faultdiagnosistoolbox.dmperm import GetDMParts
from faultdiagnosistoolbox.MHS import MHS


def ParentBlocks(X, b):
    """Internal."""
    return np.where(X[0:b, b])[0]


def FbDetectability(Xb, fb, dm):
    """Internal."""
    v = np.array([], dtype=np.int64)
    pb = np.concatenate((ParentBlocks(Xb, fb), [fb]))
    for jj in pb:
        v = np.concatenate((v, dm.M0[jj].col))
    return v


def SensorPlacementDetectability(model, **options):
    """Internal."""
    if 'fdet' in options:
        fdet = options['fdet']
        if isinstance(fdet[0], str):
            # fdet = list(map(lambda fi: model.f.index(fi), fdet))
            fdet = np.array(map(lambda fi: model.f.index(fi), fdet))
    else:
        fdet = np.arange(0, len(model.f))

    if model.IsUnderdetermined():
        print("Sorry, sensor placement only works for models with no underdetermined parts")
        return

    dm = GetDMParts(model.X)
    ef = list(map(lambda fi: np.where(model.F[:, fi])[0][0], fdet))
    nondetIdx = [eIdx for eIdx in np.arange(0, len(ef)) if ef[eIdx] in dm.M0eqs]
    if len(nondetIdx) > 0:
        ds = DetectabilitySets(model.X, model.F[:, fdet[nondetIdx]], model.P)
        sensSets = MHS(ds)
    else:
        sensSets = []

    # Compute list of variable names
    s = list(map(lambda si: list(np.array(model.x)[si]), sensSets))

    return s, sensSets


def NewSensorEqs(s, nx, nf, Pfault):
    """Internal."""
    ns = len(s)
    Xs = np.zeros((ns, nx), np.int64)
    Fs = np.zeros((ns, nf), np.int64)
    fs = []
    for sIdx, si in enumerate(s):
        Xs[sIdx, si] = 1
        if si in Pfault:
            nf = nf + 1
            Fs = np.hstack((Fs, np.zeros((ns, 1), dtype=np.int64)))  # Add column for new fault
            Fs[sIdx, -1] = 1  # Add fault
            fs.append(si)
        else:
            Fs[sIdx, :] = np.zeros((1, nf))
    return Xs, Fs, fs


def IsolabilitySubProblem(X, F, P, fi, Ispec):
    """Internal."""
    ef = np.where(F[:, fi])[0][0]
#    n = X.shape[0]
#    nf = F.shape[1]

    # Misol = M\{ef}
    Xisol = np.delete(X, ef, axis=0)
    Fisol = np.delete(F, ef, axis=0)

    # Extract just determined part of Xisol
    dm = GetDMParts(Xisol)
    X0 = Xisol[dm.M0eqs, :][:, dm.M0vars]

    # Find out which faults are included in X0
    feq = np.zeros(Fisol.shape[1], dtype=np.int64)
    for fj in np.arange(0, Fisol.shape[1]):
        e = np.where(Fisol[:, fj])[0]
        if len(e) == 0:
            e = -1
        else:
            e = e[0]
        feq[fj] = e

    nondet = np.where(list(map(lambda f: feq[f] in dm.M0eqs,
                               np.arange(0, Fisol.shape[1]))))[0]
    # Adapt to isolability specification
    nondetisol = [f for f in nondet if f in np.where(Ispec == 0)[0]]

    # nondetisol = nondet & (Ispec==0);

    F0 = Fisol[dm.M0eqs, :][:, nondetisol]

    # Translate P to P0
    P0 = [np.where(dm.M0vars == v)[0][0] for v in P if v in dm.M0vars]

    # Compute detectability sets
    detSets = DetectabilitySets(X0, F0, P0)

    # Translate back to original variable indices
    detSets = list(map(lambda ds: dm.M0vars[ds], detSets))

    return detSets


def BlockAndFaultOrder(X, F, dm):
    """Internal."""
    # 1. Construct block adjecency matrix
    #    Xb(i,j) = 1 => bi>bj
    n = len(dm.M0)
    Xb = np.zeros((n, n), dtype=np.int64)
    # 1.1 Determine connected blocks
    for rr in np.arange(0, n):
        for cc in np.arange(0, n):
            Xb[rr, cc] = np.any(X[dm.M0[rr].row, :][:, dm.M0[cc].col])

    # 1.2 Traverse block adjacency matrix to determine indirect relationships
    for bb in np.arange(1, n):
        ba = ParentBlocks(Xb, bb)
        iba = np.array([], dtype=np.int64)
        for ll in ba:
            iba = np.concatenate((iba, ParentBlocks(Xb, ll)))
        if len(iba) > 0:
            Xb[np.unique(iba), bb] = 1

    # 2. Construct fault classes and determine maximal elements
    # 2.1 Determine e_f for each fault (must be 1 equation for each fault)
    ef = np.zeros(F.shape[1], dtype=np.int64)
    for fi in np.arange(0, F.shape[1]):
        ef[fi] = np.where(F[:, fi])[0][0]

    # 2.2 Determine block membership for each fault
    efrep = np.unique(ef)
    efb = np.zeros(len(efrep), dtype=np.int64)
    for fi, efi in enumerate(efrep):
        for bjIdx, bj in enumerate(dm.M0):
            if efi in bj.row:
                efb[fi] = bjIdx

    # 2.3 Determine blocks corresponding to maximal elements in the
    #     fault class partial order
    maxFaultClasses = np.zeros(len(efb), dtype=np.int64)
    for jj in np.arange(0, len(maxFaultClasses)):
        maxFaultClasses[jj] = len([b for b in ParentBlocks(Xb, efb[jj]) if b in efb]) == 0

    bFm = efb[maxFaultClasses > 0]
    return {'blockorder': Xb, 'bFm': bFm}


def SensPlaceM0(X, F):
    """Internal."""
    dm0 = GetDMParts(X)
    bfOrder = BlockAndFaultOrder(X, F, dm0)
    return list(map(lambda b: FbDetectability(bfOrder['blockorder'], b, dm0),
                    bfOrder['bFm']))


def DetectabilitySets(X, F, P):
    """Internal."""
    dm = GetDMParts(X)
    detSets = SensPlaceM0(X[dm.M0eqs, :][:, dm.M0vars], F[dm.M0eqs, :])
    detSets = list(map(lambda d: [x for x in dm.M0vars[d] if x in P], detSets))
    return [s for s in detSets if len(s) > 0]


def SensorPlacementIsolability(model, Ispec):
    """Compute sensor placement to achive isolability requirements."""
    if not np.all(Ispec == Ispec.transpose()):
        print("Warning: Isolability not symmetric, making it so...")
        Ispec = np.bitwise_or(Ispec, Ispec.transpose())

    Pfault = model.Pfault
#    fault = np.arange(0,len(model.x))
    _, spDet = SensorPlacementDetectability(model)
    sIdx = []
    if len(spDet) > 0:
        for s in spDet:
            Xs, Fs, fs = NewSensorEqs(s, model.X.shape[1], model.F.shape[1], Pfault)
            X = model.X
            F = np.hstack((model.F, np.zeros((model.F.shape[0], len(fs)), dtype=np.int64)))
            X = np.vstack((X, Xs))
            F = np.vstack((F, Fs))
            P = model.P
            nf = Ispec.shape[0]
            Ispecii = np.vstack((
                np.hstack((Ispec, np.zeros((nf, len(fs))))),
                np.hstack((np.zeros((len(fs), nf)), np.eye(len(fs))))))

            detSets = []
            for fi in np.arange(0, F.shape[1]):
                isolDetSets = IsolabilitySubProblem(X, F, P, fi, Ispecii[fi, :])
                if len(isolDetSets) > 0:
                    detSets = detSets + isolDetSets
            mhs_s = MHS(detSets)
            sIsolIdx = list(map(lambda m: np.sort(np.concatenate((m, s))), mhs_s))
            sIdx = sIdx + sIsolIdx
    else:
        X = model.X
        F = model.F
        P = model.P

        detSets = []
        for fi in np.arange(0, F.shape[1]):
            isolDetSets = IsolabilitySubProblem(X, F, P, fi, Ispec)
            if len(isolDetSets) > 0:
                detSets = detSets + isolDetSets
        sIdx = MHS(detSets)

    # Remove duplicates and supersets from solution

    sIdx = [np.sort(si) for si in sIdx]
    keepIdx = np.full(len(sIdx), True)
    for k in np.arange(0, len(sIdx)):
        for l in np.arange(k + 1, len(sIdx)):
            if np.array_equal(sIdx[k], sIdx[l]):
                keepIdx[l] = False
    sIdx = [sIdx[k] for k in np.arange(0, len(sIdx)) if keepIdx[k]]
    # sIdx = MHS(MHS(sIdx))  # Is there a more efficient way?

    # Compute list of variable names
    s = list(map(lambda si: list(np.array(model.x)[si]), sIdx))

    return s, sIdx
