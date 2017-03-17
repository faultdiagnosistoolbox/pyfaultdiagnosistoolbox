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
detSens,_ = model.SensorPlacementDetectability()

# %% Sensor placement analysis - isolability
print("Performing sensor placement analysis...")
sens,_ = model.SensorPlacementIsolability()
print('Done! Found %d sensor sets.' % len(sens))

# %% Add 1 sensor, detectability for all faults
mDet = model.copy()
mDet.AddSensors(['pan'])

mDet.DetectabilityAnalysis()

mDet.IsLowIndex()

plt.figure(20)
plt.subplot(1,2,1)
mDet.IsolabilityAnalysis(plot=True)
plt.subplot(1,2,2)
mDet.IsolabilityAnalysis(causality='int',plot=True)

# %%
def TestDetectabilityAnalysis(model, causality='mixed'):
    """ Performs a structural detectability analysis
  
    Outputs
    -------
      df  : Detectable faults
      ndf : Non-detectable faults"""
    MplusCausal = lambda X: dmperm.Mplus(X,causality=causality)
      
    dm = MplusCausal(model.X)

    df = [model.f[fidx] for fidx in np.arange(0, model.F.shape[1]) 
        if np.argwhere(model.F[:,fidx])[:,0] in dm]
    ndf = [model.f[fidx] for fidx in np.arange(0, model.F.shape[1])
        if np.argwhere(model.F[:,fidx])[:,0] not in dm]
    return (df,ndf)
