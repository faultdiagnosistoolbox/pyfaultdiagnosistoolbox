import pytest
import numpy as np
import os
from test_engine_model import LiU_ICE_model
from faultdiagnosistoolbox.NeuralNetworkIntegration.NeuralNetworkIntegration import (
    NeuralNetworkIntegration,
)
from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import (
    GenConfigParams,
)
import helpers_test_data as data


@pytest.fixture
def nni():
    """Initializes and returns an instance of NeuralNetworkIntegration"""
    model = LiU_ICE_model
    return NeuralNetworkIntegration(model)


@pytest.fixture
def config_setup(tmp_path):
    """Creates a temporary directory, and returns the directory, filename and path"""
    test_dir = str(tmp_path / "configs")
    os.mkdir(test_dir)
    file_name = "model_1_config"
    model_type = "latent"

    path = os.path.join(test_dir, file_name)
    full_path = path + f"_{model_type}.yml"

    yield test_dir, file_name, full_path, model_type

    # Cleanup
    if os.path.exists(full_path):
        os.remove(full_path)


def _get_complete_mso():
    """
    Return mso 585 (index 584) from a set of msos
    """
    return LiU_ICE_model.MSO()[584]


def _get_mso_585_inputs() -> tuple[list[str], list[str], list[str]]:
    """
    Return inputs for mso 584 in the LiU_ICE_model
    """
    # Inputs
    signals = [
        "intercooler_pressure",
        "intercooler_temperature",
        "intake_manifold_pressure",
        "air_mass_flow",
        "engine_speed",
        "throttle_position",
        "wastegate_position",
        "injected_fuel_mass",
        "ambient_temperature",
        "ambient_pressure",
    ]

    y_variables = [
        "y_p_ic",
        "y_T_ic",
        "y_p_im",
        "y_W_af",
        "y_omega_e",
        "y_alpha_th",
        "y_u_wg",
        "y_wfc",
        "y_T_amb",
        "y_p_amb",
    ]

    x_variables = [
        "p_ic",
        "T_ic",
        "p_im",
        "W_af",
        "omega_e",
        "alpha_th",
        "u_wg",
        "wfc",
        "T_amb",
        "p_amb",
    ]
    return signals, y_variables, x_variables


class FakeModel:
    def __init__(self):
        self.syme = ["eq1", "eq2"]
        self.X = np.array([[1, 0], [0, 1]])
        self.Z = np.array([[1, 0], [0, 1]])
        self.x = ["x1", "x2"]
        self.z = ["z1", "z2"]
        self.name = "fake_model"


class FakeMatching:
    def __init__(self):
        self.matching = [self]

    @property
    def row(self):
        return [0, 1]

    @property
    def col(self):
        return [0, 1]

    def __iter__(self):
        return iter([self])


class FakeParams:
    def __init__(self):
        self.signals = ["s1", "s2"]
        self.y_vars = ["y1", "y2"]
        self.x_vars = ["x1", "x2"]


def test_bad_init_input(subtests):
    """Test input validation for __init__ function"""
    with subtests.test("test init with input None"):
        with pytest.raises(ValueError):
            NeuralNetworkIntegration(None)


def test_config_file_input_validation(subtests):
    """Test input validation in generate_config_file"""
    model = FakeModel()
    nni = NeuralNetworkIntegration(model)

    with subtests.test("Bad res_eq input"):
        with pytest.raises(TypeError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq="bad_idx",
                params=FakeParams(),
                file_name="test",
                type="latent",
            )

        with pytest.raises(IndexError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq=10,
                params=FakeParams(),
                file_name="test",
                type="latent",
            )

    with subtests.test("Bad gamma input"):
        with pytest.raises(ValueError):
            nni.generate_config_file(
                gamma=None,
                path=".",
                res_eq=0,
                params=FakeParams(),
                file_name="test",
                type="latent",
            )
    with subtests.test("Invalid model type"):
        with pytest.raises(ValueError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq=0,
                params=FakeParams(),
                file_name="test",
                type="invalid_type",
            )

    with subtests.test("Bad params"):
        with pytest.raises(AttributeError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq=0,
                params={},
                file_name="test",
                type="greybox",
            )

        bad_params = FakeParams()
        bad_params.x_vars = []
        with pytest.raises(ValueError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq=0,
                params=bad_params,
                file_name="test",
                type="greybox",
            )


def test_config_file_error_handling(subtests, monkeypatch, config_setup):
    model = FakeModel()
    nni = NeuralNetworkIntegration(model)
    with subtests.test("Bad ssmodel returns"):
        monkeypatch.setattr(nni, "get_int_causality_model", lambda *a, **k: [])

        with pytest.raises(RuntimeError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq=0,
                params=FakeParams(),
                file_name="test",
                type="greybox",
            )

        class FakeElem:
            def __init__(self):
                self.elem_type = "out"
                self.xout = ["x"]
                self.zin = []
                self.xin = []
                self.xout_int = ["x"]

        monkeypatch.setattr(
            nni,
            "get_int_causality_model",
            lambda *a, **k: [FakeElem(), FakeElem(), FakeElem()],
        )

        with pytest.raises(RuntimeError):
            nni.generate_config_file(
                gamma=FakeMatching(),
                path=".",
                res_eq=0,
                params=FakeParams(),
                file_name="test",
                type="greybox",
            )

    with subtests.test("Try creating file with taken name"):
        test_dir, file_name, _, model_type = config_setup
        params = GenConfigParams(
            _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
        )
        m = LiU_ICE_model
        real_nni = NeuralNetworkIntegration(m)
        mso = _get_complete_mso()
        gamma = m.Matching(np.setdiff1d(mso, mso[3]))

        real_nni.generate_config_file(
            gamma, test_dir, mso[3], params, file_name, model_type
        )

        with pytest.raises(ValueError):
            real_nni.generate_config_file(
                gamma, test_dir, mso[3], params, file_name, model_type
            )
