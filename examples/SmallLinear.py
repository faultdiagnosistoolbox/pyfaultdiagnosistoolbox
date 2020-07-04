import matplotlib.pyplot as plt
import faultdiagnosistoolbox as fdt

# %% Define structural model
model_def= {'type': 'MatrixStruc',
            'X': [[1, 1, 0, 0, 1],
                  [0, 1, 1, 1, 0],
                  [0, 0, 1, 0, 1],
                  [0, 0, 0, 1, 1],
                  [0, 0, 0, 0, 1]],
            'F': [[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [1, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]],
            'Z': []}
model = fdt.DiagnosisModel(model_def, name='Small linear model')
model.Lint()

# Plot structural model
plt.figure(10)
model.PlotModel()

# %% Determine sensors need to make faults detectable
sDet, _ = model.SensorPlacementDetectability()
print(f"Sensor sets: {sDet}")
model2 = model.copy()
model2.AddSensors(sDet[0])
df, ndf = model2.DetectabilityAnalysis()
print(f"Detectable faults: {df}")


# Determine sensors needed to make faults isolable when new sensors can not fail
model.SensorLocationsWithFaults() # No new sensors may become faulty
sIsol, _ = model.SensorPlacementIsolability()
print(f"Sensor sets: {sIsol}")

model3 = model.copy()
model3.AddSensors(sIsol[0])
plt.figure(20)
_ = model3.IsolabilityAnalysis(plot=True)

plt.figure(21)
model3.PlotDM(eqclass=True, fault=True)


# Determine sensors needed to make faults isolable when all new sensors can fail
model.SensorLocationsWithFaults(model.x) # Set all new sensors may become faulty
sIsol, _ = model.SensorPlacementIsolability()
print(f"Sensor sets: {sIsol}")

model4 = model.copy()
model4.AddSensors(sIsol[0])
plt.figure(30)
model4.IsolabilityAnalysis(plot=True)

plt.figure(31)
model4.PlotDM(eqclass=True, fault=True)

plt.show()
