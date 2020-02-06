import numpy as np
import faultdiagnosistoolbox as fdt
import pytest

# Define the model structure
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
