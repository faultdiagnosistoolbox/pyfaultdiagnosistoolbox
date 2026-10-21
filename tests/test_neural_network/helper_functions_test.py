import pytest
import numpy as np
from test_engine_model import LiU_ICE_model
from faultdiagnosistoolbox.NeuralNetworkIntegration.NeuralNetworkIntegration import (
    NeuralNetworkIntegration,
)
from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import (
    SsmodelElement,
)
import helpers_test_data as data


@pytest.fixture
def nni():
    """Initializes and returns an instance of NeuralNetworkIntegration"""
    model = LiU_ICE_model
    return NeuralNetworkIntegration(model)


def test_get_icm(subtests, nni):
    matching, x, z, expected_out = data.TEST_ALL

    class match:
        def __init__(self, c, r):
            self.col = c
            self.row = r

    X = nni.model.X
    Z = nni.model.Z
    nni.model.X = x
    nni.model.Z = z
    test_match = [match(m["col"], m["row"]) for m in matching]

    test_ssmodel = nni.get_int_causality_model(test_match, 85)
    expected_xout = [o["xout"] for o in expected_out]

    assert test_ssmodel is not None, "SSModel should not be None"

    with subtests.test("Expected output"):
        for e in range(len(test_ssmodel)):
            test = test_ssmodel[e]

            with subtests.test("Expected output", e=e):
                assert (
                    test.xout in expected_xout
                ), f"XOUT: Expected xout to be one of {expected_xout} but got {test.xout}"
                assert isinstance(
                    test, SsmodelElement
                ), f"Each element in SSModel should be an instance of SsmodelElement but got {type(test)} at index {e}"

            ret = expected_out[expected_xout.index(test.xout)]

            with subtests.test("Expected x input", e=e):
                assert len(test.xin) == len(
                    ret["xin"]
                ), f"XIN: Expected SSmodel nr {e} to have {len(ret["xin"])} elements but have {len(test.xin)}"
                assert len(test.xin) == len(
                    [x for x in test.xin if x in ret["xin"]]
                ), f"XIN: Expected xin to contain {ret['xin']} but got {test.xin}"

            with subtests.test("Expected z input", e=e):
                assert len(test.zin) == len(
                    ret["zin"]
                ), f"ZIN: Expected SSmodel nr {e} to have {len(ret["zin"])} elements but have {len(test.zin)}"

                assert len(test.zin) == len(
                    [z for z in test.zin if z in ret["zin"]]
                ), f"ZIN: Expected zin to contain {ret['zin']} but got {test.zin}"

            with subtests.test("Expected element type", e=e):
                assert (
                    test.elem_type == ret["type"]
                ), f"Expected elem_type to be {ret['type']} but got {test.elem_type}"

            with subtests.test("Expected xout_int", e=e):
                assert len(test.xout_int) == len(
                    ret["xout_int"]
                ), f"XOUT_INT: Expected SSmodel nr {e} to have {len(ret["xout_int"])} elements but have {len(test.xout_int)}"
                assert len(test.xout_int) == len(
                    [x for x in test.xout_int if x in ret["xout_int"]]
                ), f"XOUT_INT: Expected xin to contain {ret['xout_int']} but got {test.xout_int}"

    # Clean up
    nni.model.X = X
    nni.model.Z = Z


def test_find_state_dependencies_for_state_variable(subtests, nni):
    # Test data
    X, Z, expected_state_vars, expected_sensor_vars, idx, causality = data.FSD_INT_DATA
    nni.x0 = np.asarray(X)
    nni.z0 = np.asarray(Z)

    state_vars, sensor_vars, xx = nni.find_state_dependencies(idx)
    print("RECURSION: ", xx)
    assert len(state_vars) == len(expected_state_vars)
    assert len(sensor_vars) == len(expected_sensor_vars)

    with subtests.test("Expected state variables for variable with int causality"):
        for idx, var in enumerate(expected_state_vars):
            assert (
                var == state_vars[idx]
            ), f"Expected state variable at index {idx} to be {var} got {state_vars[idx]}"

    with subtests.test("Expected sensor variables for variable with int causality"):
        for idx, var in enumerate(expected_sensor_vars):
            assert (
                var == sensor_vars[idx]
            ), f"Expected state variable at index {idx} to be {var} got {sensor_vars[idx]}"


def test_find_state_dependencies_for_output(subtests, nni):
    X, Z, expected_state_vars, expected_sensor_vars, idx, causality = data.FSD_OUT_DATA
    nni.x0 = np.asarray(X)
    nni.z0 = np.asarray(Z)

    state_vars, sensor_vars, xx = nni.find_state_dependencies(idx)
    print("RECURSION: ", xx)
    assert len(state_vars) == len(expected_state_vars)
    assert len(sensor_vars) == len(expected_sensor_vars)

    with subtests.test("Expected state variables for variable with int causality"):
        for idx, var in enumerate(expected_state_vars):
            assert (
                var == state_vars[idx]
            ), f"Expected state variable at index {idx} to be {var} got {state_vars[idx]}"

    with subtests.test("Expected sensor variables for variable with int causality"):
        for idx, var in enumerate(expected_sensor_vars):
            assert (
                var == sensor_vars[idx]
            ), f"Expected state variable at index {idx} to be {var} got {sensor_vars[idx]}"


def test_get_sub_matricies(subtests, nni):
    x, z, expected_x0, expected_z0, eqs, var = data.SUB_MATRIX_DATA

    class sub_matrix:
        X = x
        Z = z

    nni.equations = eqs
    nni.variables = var

    nni.generate_sub_matrices()

    with subtests.test(msg="Expected X0 and Z0 dimensions"):
        assert nni.x0.shape[0] == len(expected_x0)
        assert nni.z0.shape[0] == len(expected_z0)
        assert nni.x0.shape[1] == len(expected_x0[0])
        assert nni.z0.shape[1] == len(expected_z0[0])

    with subtests.test(msg="Expected X0 data"):
        for r in range(nni.x0.shape[0]):
            for c in range(nni.x0.shape[1]):
                assert nni.x0[r][c] == expected_x0[r][c]

    with subtests.test(msg="Expected Z0 data"):
        for r in range(nni.z0.shape[0]):
            for c in range(nni.z0.shape[1]):
                assert nni.z0[r][c] == expected_z0[r][c]


def test_analyze_matching(subtests, nni):
    matching, expected_eq, expected_var = data.TEST_MATHCING_DATA

    class match:
        def __init__(self, c, r):
            self.col = c
            self.row = r

    test_match = [match(m["col"], m["row"]) for m in matching]

    nni.analyze_matching(test_match)

    assert len(nni.equations) == len(expected_eq)
    assert len(nni.variables) == len(expected_var)

    with subtests.test(msg="Expected equations"):
        for e in range(len(nni.equations)):
            assert nni.equations[e] == expected_eq[e]

    with subtests.test(msg="Expected variables"):
        for v in range(len(nni.variables)):
            assert nni.variables[v] == expected_var[v]
