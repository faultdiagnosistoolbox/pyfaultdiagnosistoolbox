import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os
import sys
from importlib import reload
new_paths = ['../Misc/', '../src', '../src/faultdiagnosistoolbox/']
[sys.path.append(d) for d in new_paths if not d in sys.path];

import numpy as np
import scipy
from scipy.stats import gaussian_kde
from misc import *
import GetMeasurementData as gm


import VEP4Engine
model = VEP4Engine.model;

dataDir = '/Users/frisk/Diagnos/Work/EngineDiagnosis/Work/Data/'
dataSets = {'NF' : 'driving_cycle_FTP75_highway_no_fault_dataset1_16-01-20.mat',
            'fyp_im' : 'driving_cycle_FTP75Highway_fault_y_pim_dataset_7_26-01-2016.mat',
            'fyw_af' : 'driving_cycle_FTP75Highway_fault_y_waf_dataset_5_26-01-2016.mat',
            'fyp_ic' : 'driving_cycle_FTP75Highway_fault_y_pic_dataset_8_26-01-2016.mat',
            'fyT_ic' : 'driving_cycle_FTP75Highway_fault_y_Tic_dataset_6_26-01-2016.mat'}
data = {}
sys.stdout.write('Loading data ')
for dd in dataSets:
    sys.stdout.write('.')
    data[dd] = gm.GetMeasurementData(dataDir + dataSets[dd])
sys.stdout.write(' Finished loading ' + str(len(data.keys())) + ' datasets\n')

import ResGen_1649_1
r1 = RunResgenOnDataSets( ResGen_1649_1.ResGen_1649_1, data, VEP4Engine.diag_par )


