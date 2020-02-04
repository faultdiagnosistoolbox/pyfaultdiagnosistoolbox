import csv
import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import sys
import matplotlib.pyplot as plt
import time
import matplotlib.colors as mcolors

def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

class Timer:
    t0 = 0.0
    def tic(self):
        self.t0 = time.time()
    def toc(self):
        return time.time()-self.t0

def RunResgenOnDataSets( resgen, data, params ):
    res = {}
    t = Timer()
    for ds in data.keys():
        sys.stdout.write('  dataset ' + ds)
        state = data[ds]['state_init'].copy()
        Ts = data[ds]['Ts']
        z = data[ds]['z']
        N = len(data[ds]['time'])
        t.tic()        
        res[ds] = resgen(z,state,params,Ts)
        dt=t.toc()
        sys.stdout.write(' (%.2f sec, %.1f x real time @ %d Hz)\n' % (dt, Ts*N/dt, 1/Ts))
        sys.stdout.flush()
    return res

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

summer_cmap = make_colormap([
    (0.000, 0.500, 0.400),
    (0.016, 0.508, 0.400),
    (0.032, 0.516, 0.400),
    (0.048, 0.524, 0.400),
    (0.063, 0.532, 0.400),
    (0.079, 0.540, 0.400),
    (0.095, 0.548, 0.400),
    (0.111, 0.556, 0.400),
    (0.127, 0.563, 0.400),
    (0.143, 0.571, 0.400),
    (0.159, 0.579, 0.400),
    (0.175, 0.587, 0.400),
    (0.190, 0.595, 0.400),
    (0.206, 0.603, 0.400),
    (0.222, 0.611, 0.400),
    (0.238, 0.619, 0.400),
    (0.254, 0.627, 0.400),
    (0.270, 0.635, 0.400),
    (0.286, 0.643, 0.400),
    (0.302, 0.651, 0.400),
    (0.317, 0.659, 0.400),
    (0.333, 0.667, 0.400),
    (0.349, 0.675, 0.400),
    (0.365, 0.683, 0.400),
    (0.381, 0.690, 0.400),
    (0.397, 0.698, 0.400),
    (0.413, 0.706, 0.400),
    (0.429, 0.714, 0.400),
    (0.444, 0.722, 0.400),
    (0.460, 0.730, 0.400),
    (0.476, 0.738, 0.400),
    (0.492, 0.746, 0.400),
    (0.508, 0.754, 0.400),
    (0.524, 0.762, 0.400),
    (0.540, 0.770, 0.400),
    (0.556, 0.778, 0.400),
    (0.571, 0.786, 0.400),
    (0.587, 0.794, 0.400),
    (0.603, 0.802, 0.400),
    (0.619, 0.810, 0.400),
    (0.635, 0.817, 0.400),
    (0.651, 0.825, 0.400),
    (0.667, 0.833, 0.400),
    (0.683, 0.841, 0.400),
    (0.698, 0.849, 0.400),
    (0.714, 0.857, 0.400),
    (0.730, 0.865, 0.400),
    (0.746, 0.873, 0.400),
    (0.762, 0.881, 0.400),
    (0.778, 0.889, 0.400),
    (0.794, 0.897, 0.400),
    (0.810, 0.905, 0.400),
    (0.825, 0.913, 0.400),
    (0.841, 0.921, 0.400),
    (0.857, 0.929, 0.400),
    (0.873, 0.937, 0.400),
    (0.889, 0.944, 0.400),
    (0.905, 0.952, 0.400),
    (0.921, 0.960, 0.400),
    (0.937, 0.968, 0.400),
    (0.952, 0.976, 0.400),
    (0.968, 0.984, 0.400),
    (0.984, 0.992, 0.400),
    (1.000, 1.000, 0.400)])

def PlotConfusionMatrix(C):
    """Plot a confusion matrix in a suitable colormap."""
    nf = C.shape[0]
    plt.imshow(C, cmap=summer_cmap)
    for fi in range(nf):
        for fj in range(nf):
            plt.text(fi, fj, '%.1f' % (C[fj, fi]*100), ha='center', va='center', color='k')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xticks(np.arange(C.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(C.shape[0]+1)-.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)