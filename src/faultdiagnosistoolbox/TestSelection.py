"""Test selection functionality."""

import numpy as np
from faultdiagnosistoolbox.MHS import MHS


def TestSelection(self, arr, method='aminc', isolabilitymatrix=''):
    """Select set of arrs based on structural information."""
    if len(isolabilitymatrix) == 0:
        im = self.IsolabilityAnalysisArrs(arr)
    else:
        im = isolabilitymatrix
    FSM = self.FSM(arr)
    if method == 'aminc':
        return np.sort(aminc(TestSets(FSM, im), FSM.shape[0]))
    elif method == 'full':
        ts = np.array(MHS(TestSets(FSM, im)))
        sort_arg = np.argsort([len(ti) for ti in ts])
        return [np.sort(ts_i) for ts_i in ts[sort_arg]]
#         return np.sort(MHS(TestSets(FSM, im), FSM.shape[0]))
    else:
        print('Unknown test selection method')
        return


def TestSets(FSM, im):
    """internal."""
    isolProperty = np.argwhere(im == 0)
    ts = []
    for iprop in isolProperty:
        ts = ts + list(np.where(np.logical_and(FSM[:, iprop[1]] == 1,
                                               FSM[:, iprop[0]] == 0)))
    return ts


def aminc(pi, n):
    """Find, appriximative, minimal cardinality minimal hitting set."""
    piM = np.zeros((len(pi), n))
    for k, p in enumerate(pi):
        piM[k, p] = 1
    hs = []
    while np.count_nonzero(piM) > 0:
        tIdx = np.argmax(np.sum(piM, axis=0, dtype=np.int64))
        hs.append(tIdx)
        piM[piM[:, tIdx] == 1, :] = 0
    return hs
