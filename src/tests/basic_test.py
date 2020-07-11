import numpy as np
import faultdiagnosistoolbox as fdt
import sympy as sym

# Define the model structure
modelDef = {
    'type': 'MatrixStruc',
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

model = fdt.DiagnosisModel(modelDef, name='Example from Commault et.al')

pendulum_def = {'type': 'Symbolic',
                'x': ['dx', 'dy', 'dw', 'dz', 'T', 'x', 'y', 'w', 'z'],
                'f': [], 'z': [],
                'parameters': ['g', 'L']}

sym.var(pendulum_def['x'])
sym.var(pendulum_def['parameters'])

pendulum_def['rels'] = [
    -dx + w,
    -dy + z,
    -dw + T * x,
    -dz + T * y - g,
    x**2 + y**2 - L**2,
    fdt.DiffConstraint('dx', 'x'),
    fdt.DiffConstraint('dy', 'y'),
    fdt.DiffConstraint('dw', 'w'),
    fdt.DiffConstraint('dz', 'z')]

pendulum = fdt.DiagnosisModel(pendulum_def, name='Pendulum')


def test_isolabilityanalysis():
    assert np.all(model.IsolabilityAnalysis() == 
                  np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [0, 0, 1, 1]]))


sensSets, idx = model.SensorPlacementIsolability()


def test_sensorplacement():
    assert len(sensSets) == len(idx)
    assert sensSets == [['x6', 'x5'], ['x6', 'x7']]
    assert np.all([list(np.array(model.x)[id_i]) == sens_i for (sens_i, id_i) in zip(sensSets, idx)])


model2 = model.copy()
model2.AddSensors(sensSets[0])


def test_isolabilityanalysis_new_sensors():
    assert np.all(model2.IsolabilityAnalysis() == np.array([[1, 0, 0, 0],
                                                            [0, 1, 0, 0],
                                                            [0, 0, 1, 0],
                                                            [0, 0, 0, 1]]))


def test_pantelides_algorithm():
    s_idx, nu = pendulum.Pantelides()
    assert s_idx == 3
    assert np.all(nu == np.array([1, 1, 0, 0, 2, 1, 1, 0, 0]))
