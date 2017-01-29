import numpy as np
import matplotlib.pyplot as plt
import sys
new_paths = ['../Misc/', '../src/','../src/faultdiagnosistoolbox/']

[sys.path.append(d) for d in new_paths if not d in sys.path]

from misc import *
import faultdiagnosistoolbox as fdt


# %% Define the model structure
modelDef={}
modelDef['type'] = 'MatrixStruc'
modelDef['X'] = [
    [0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 1, 0, 1],
    [0, 1, 0, 0, 0, 1, 0],
    [1, 1, 1, 0, 0, 0, 0],
    [1, 1, 0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0],
    [0, 0, 0, 1, 0, 0, 0]]

modelDef['F'] = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 0, 1],
    [0, 0, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [0, 0, 0, 0]]

modelDef['Z'] = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]]
model = fdt.DiagnosisModel(modelDef, name='Example from Commault et.al')

model.Lint()


# %% Plot model structure
plt.figure(10)
model.PlotModel()


# %% Perform isolability analysis

plt.figure(20)
model.IsolabilityAnalysis(plot=True);


sensSets,_ = model.SensorPlacementIsolability()
print("Found " + str(len(sensSets)) + " sensor sets")


# %% Add first sensor set and redo isolability analysis

model2 = model.copy()
model2.AddSensors(sensSets[0])
plt.figure(30)
model2.IsolabilityAnalysis(plot=True);

