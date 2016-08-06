import csv
import numpy as np
import scipy.sparse as sp
import scipy.io as sio
import sys
import matplotlib.pyplot as plt
import time

def ReadMatrixCSV(fileName):
    with open(fileName, 'rt') as csvfile:
        matReader = csv.reader(csvfile, delimiter=' ', quotechar='|')
        size = matReader.next() # Assume first line contains size information
        X=np.zeros((int(size[0]),int(size[1])),dtype='float64');
        for edge in matReader:
            X[int(edge[0]),int(edge[1])] = float(edge[2])
        return sp.csc_matrix(X.astype(np.int64))

def BoxOff(*argin):
    if len(argin)>0:
        ax=argin[0]
    else:
        ax=plt.gca();
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


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
