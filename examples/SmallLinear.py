# %% Initial imports
import numpy as np
import matplotlib.pyplot as plt
import sys
new_paths = ['../Misc/', '../src/','../src/faultdiagnosistoolbox/']
[sys.path.append(d) for d in new_paths if not d in sys.path]
from misc import *
import faultdiagnosistoolbox as fdt
import scipy.sparse as sp

# %% Define structural model

modelDef={}
modelDef['type'] = 'MatrixStruc'
modelDef['X'] = [
    [1, 1, 0, 0, 1],
    [0, 1, 1, 1, 0],
    [0, 0, 1, 0, 1],
    [0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1]]
modelDef['F'] = [
    [0, 0, 0, 0],
    [0, 0, 0, 0],
    [1, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]]
modelDef['Z'] = []
model = fdt.DiagnosisModel(modelDef, name='Small linear model')
model.Lint()

# %% Plot structural model
plt.figure(10)
model.PlotModel()

# %% Determine sensors need to make faults detectable
sDet,_ = model.SensorPlacementDetectability()
print("Sensor sets: " + str(sDet))
model2 = model.copy()
model2.AddSensors(sDet[0])
df,ndf = model2.DetectabilityAnalysis()
print("Detectable faults: " + str(df))


# %% Determine sensors needed to make faults isolable when new sensors can not fail
model.SensorLocationsWithFaults() # No new sensors may become faulty
sIsol,_ = model.SensorPlacementIsolability()
print("Sensor sets: ")
print(sIsol)

model3 = model.copy()
model3.AddSensors(sIsol[0])
plt.figure(20)
model3.IsolabilityAnalysis(plot=True);

plt.figure(21)
model3.PlotDM(eqclass=True,fault=True)


# %% Determine sensors needed to make faults isolable when all new sensors can fail
model.SensorLocationsWithFaults(model.x) # Set all new sensors may become faulty
sIsol,_ = model.SensorPlacementIsolability()
print("Sensor sets: ")
print(sIsol)

model4 = model.copy()
model4.AddSensors(sIsol[0])
plt.figure(30)
model4.IsolabilityAnalysis(plot=True);

plt.figure(31)
model4.PlotDM(eqclass=True,fault=True)

