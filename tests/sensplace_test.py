import faultdiagnosistoolbox as fdt
import sympy as sym
import numpy as np

modelDef = {
    "type": "Symbolic",
    "x": ["x1", "x2", "x3", "x4", "x5", "dx1", "dx2", "dx3", "dx4", "dx5"],
    "f": ["f1", "f2", "f3", "f4"],
    "z": ["u"],
}

x1, x2, x3, x4, x5, dx1, dx2, dx3, dx4, dx5 = [sym.Symbol(vi) for vi in modelDef["x"]]
f1, f2, f3, f4 = [sym.Symbol(vi) for vi in modelDef["f"]]
u = sym.Symbol(modelDef["z"][0])

modelDef["rels"] = [
    -dx1 - x1 + x2 + x5,
    -dx2 - 2 * x2 + x3 + x4,
    -dx3 - 3 * x3 + x5 + f1 + f2,
    -dx4 - 4 * x4 + x5 + f3,
    -dx5 - 5 * x5 + u + f4,
    fdt.DiffConstraint("dx1", "x1"),
    fdt.DiffConstraint("dx2", "x2"),
    fdt.DiffConstraint("dx3", "x3"),
    fdt.DiffConstraint("dx4", "x4"),
    fdt.DiffConstraint("dx5", "x5"),
]

model = fdt.DiagnosisModel(modelDef, name="Linear model")
model.PossibleSensorLocations(["x1", "x2", "x3", "x4", "x5"])

model_lumped = model.copy()
model_lumped.Structural()
model_lumped.LumpDynamics()


def test_model_lumping():
    assert model.Redundancy() == 0
    assert model_lumped.modelType == "Structural"
    assert model_lumped.IsStatic()
    assert model_lumped.ne() == 5 and model_lumped.nx() == 5
    assert model_lumped.Redundancy() == 0


def test_sensplace():
    sens_isol, sens_isol_idx = model.SensorPlacementIsolability()
    sens_lumped_isol, sens_lumped_isol_idx = model_lumped.SensorPlacementIsolability()

    sens_isol_idx = [np.sort(si) for si in sens_isol_idx]
    sens_lumped_isol_idx = [np.sort(si) for si in sens_lumped_isol_idx]

    assert len(sens_isol_idx) == 5

    assert len(sens_isol_idx) == len(sens_isol)
    assert len(sens_lumped_isol_idx) == len(sens_lumped_isol)
    assert len(sens_isol_idx) == len(sens_lumped_isol_idx)

    for si1 in sens_isol_idx:
        assert np.any(
            [(len(si1) == len(si) and np.all(si1 == si)) for si in sens_lumped_isol_idx]
        )


def test_sensplace_faults():
    model.SensorLocationsWithFaults(["x1", "x2", "x3", "x4", "x5"])
    model_lumped.SensorLocationsWithFaults(["x1", "x2", "x3", "x4", "x5"])

    sens_isol, sens_isol_idx = model.SensorPlacementIsolability()
    sens_lumped_isol, sens_lumped_isol_idx = model_lumped.SensorPlacementIsolability()

    assert len(sens_isol_idx) == 9

    sens_isol_idx = [np.sort(si) for si in sens_isol_idx]
    sens_lumped_isol_idx = [np.sort(si) for si in sens_lumped_isol_idx]

    assert len(sens_isol_idx) == len(sens_isol)
    assert len(sens_lumped_isol_idx) == len(sens_lumped_isol)
    assert len(sens_isol_idx) == len(sens_lumped_isol_idx)

    for si1 in sens_isol_idx:
        assert np.any(
            [(len(si1) == len(si) and np.all(si1 == si)) for si in sens_lumped_isol_idx]
        )
