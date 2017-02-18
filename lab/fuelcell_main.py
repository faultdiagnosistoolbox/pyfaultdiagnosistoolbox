# %matplotlib
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

new_paths = ['../src', '../Misc']
[sys.path.append(d) for d in new_paths if not d in sys.path];
 
from misc import BoxOff

# %% Define the model and basic validity check
import PEMFCModel
model = PEMFCModel.model;
possSensLocations = ['Pcmp', 'Tcain', 'Tcmp', 'Wcmp', 'Vfc', 'Tfc', 'pca',
                     'pan', 'psmca', 'psman']
model.PossibleSensorLocations( possSensLocations )

# %% Basic validity check
model.Lint()

# %% Plot the structural model
plt.figure(10)
model.PlotModel(verbose=False)

# %% Sensor placement analysis - detectability
detSens = model.SensorPlacementDetectability()

# %% Sensor placement analysis - isolability
sens = model.SensorPlacementIsolability()
