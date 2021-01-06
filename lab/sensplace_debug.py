# %load_ext autoreload
# %autoreload 2

import faultdiagnosistoolbox as fdt

modelDef = {'type': 'VarStruc',
            'x': ['x1', 'x2', 'x3', 'x4', 'x5'],
            'f': ['f1', 'f2', 'f3', 'f4'],
            'z': ['u'],
            'rels': [
                ['x1', 'x2', 'x5'],
                ['x2', 'x3', 'x4'],
                ['x3', 'x5', 'f1', 'f2'],
                ['x4', 'x5', 'f3'],
                ['x5', 'u', 'f4']]}
model = fdt.DiagnosisModel(modelDef, name='Simple model')
model.PossibleSensorLocations(model.x)

sens_isol, _ = model.SensorPlacementIsolability()
print(sens_isol)

model2 = model.copy()
#model2.SensorLocationsWithFaults() # No faulty sensors
model2.SensorLocationsWithFaults(model.x) # All sensors may fail
sens2_isol, _ = model2.SensorPlacementIsolability()
print(sens2_isol)
