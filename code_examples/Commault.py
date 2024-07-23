import matplotlib.pyplot as plt
import faultdiagnosistoolbox as fdt

# Sensor placement analysis of the model in
#  "Structural analysis for the sensor location problem
#  in fault detection and isolation" by C. Commault,
#  J. Dion and S.Y. Agha, Proceedings of Safeprocess'06
#  Beijing, China.
#
#  To compare results, remember the difference that in this
#  approach, in contrast to the Commault et.al. paper,
#  the fault signals are not considered measurable.

# %% Define the model structure
model_def = {'type': 'MatrixStruc',
             'X': [[0, 0, 0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 0, 1],
                   [0, 0, 0, 0, 0, 1, 0],
                   [1, 0, 0, 0, 1, 0, 1],
                   [0, 1, 0, 0, 0, 1, 0],
                   [1, 1, 1, 0, 0, 0, 0],
                   [1, 1, 0, 1, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0],
                   [0, 0, 1, 1, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0]],
             'F': [[1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 1],
                   [0, 0, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 0]],
             'Z': [[0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [0, 0, 0],
                   [1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1]]}

model = fdt.DiagnosisModel(model_def, name='Example from Commault et.al')

model.Lint()

# Plot model structure
_, ax = plt.subplots(num=10)
model.PlotModel(ax=ax, verbose=True)

# Perform isolability analysis
_, ax = plt.subplots(num=20)
model.IsolabilityAnalysis(ax=ax)

# Perform sensor placement analysis
sensSets, _ = model.SensorPlacementIsolability()
print(f"Found {len(sensSets)} sensor sets")


# Add first sensor set and redo isolability analysis
model2 = model.copy()
model2.AddSensors(sensSets[0])
_, ax = plt.subplots(num=30)
_ = model2.IsolabilityAnalysis(ax=ax)

plt.show()
