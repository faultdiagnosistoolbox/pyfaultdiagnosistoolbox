# %matplotlib
import matplotlib.pyplot as plt
from seaborn import despine

# %% Define the model and basic validity check
from PEMFCModel import model  # noqa
possSensLocations = ['Pcmp', 'Tcain', 'Tcmp', 'Wcmp', 'Vfc', 'Tfc', 'pca',
                     'pan', 'psmca', 'psman']
model.PossibleSensorLocations(possSensLocations)

# %% Basic validity check
model.Lint()

# %% Plot the structural model
_, ax = plt.subplots(num=10, clear=True)
model.PlotModel(ax=ax, verbose=False)
despine(ax=ax)


# %% Sensor placement analysis - detectability
detSens, _ = model.SensorPlacementDetectability()

# %% Sensor placement analysis - isolability
print("Performing sensor placement analysis...")
sens, _ = model.SensorPlacementIsolability()
print(f'Done! Found {len(sens)} sensor sets.')

# %% Add 1 sensor, detectability for all faults
mDet = model.copy()
mDet.AddSensors(['pan'])

mDet.DetectabilityAnalysis()

mDet.IsLowIndex()

_, ax = plt.subplots(1, 2, num=20, clear=True)
mDet.IsolabilityAnalysis(ax=ax[0])
mDet.IsolabilityAnalysis(ax=ax[1], causality='int')
