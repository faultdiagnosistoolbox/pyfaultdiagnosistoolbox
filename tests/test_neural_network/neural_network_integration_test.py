import pytest
import os
import yaml
import time
from faultdiagnosistoolbox.NeuralNetworkIntegration import NeuralNetworkIntegration
from test_engine_model import LiU_ICE_model


@pytest.fixture
def nni():
    """Initializes and returns an instance of NeuralNetworkIntegration"""
    return NeuralNetworkIntegration()


@pytest.fixture
def config_setup(tmp_path):
    """Creates a temporary directory, and returns the directory, filename and path"""
    test_dir = str(tmp_path / _get_config_path())
    os.mkdir(test_dir)
    file_name = _get_config_name()
    full_path = os.path.join(test_dir, file_name)

    yield test_dir, file_name, full_path

    # Cleanup
    if os.path.exists(full_path):
        os.remove(full_path)


def _get_config_path():
    """Returns a test config folder path name"""
    config_path = "configs"
    return config_path


def _get_config_name():
    """Returns a test config file name"""
    name = "model_1_config.yaml"
    return name


def _get_complete_mso():
    """
    Return a mso from a set of msos
    """
    return LiU_ICE_model.MSO()[0]


def test_config_file_created_at_destination(config_setup, nni):
    """
    Tests that the file is generated at the correct given destination and with the correct filename.

    Solves: Issue #31
    """
    # Setup for the test path and file
    test_dir, file_name, expected_full_path = config_setup

    mso = _get_complete_mso()

    # Call for the config file generator
    nni.generate_config_file(mso=mso, name=file_name, path=test_dir)

    # Check that the file and path was created at correct locations and with correct names
    assert os.path.exists(
        expected_full_path
    ), f"The path was not created correctly at: {expected_full_path}!"
    assert os.path.isfile(
        expected_full_path
    ), f"Expected the file to be created at {expected_full_path}, but the file was not found!"


def test_analyze_MSO(nni):
    """
    Tests should verify correct inputs (parameters, variables, derivatives and integrals)
    for the config-file.

    Solves: Issue #26 - test analyze mso
    """
    msos = LiU_ICE_model.MSO()
    assert len(msos) > 0
    expected_keys = {"variables", "parameters", "derivatives", "integrals"}

    # Iterates through every mso
    for mso in msos:
        result = nni.analyze_mso(mso)

        assert isinstance(result, dict)

        for key in expected_keys:
            assert key in result, f"{key} missing from analyze_mso result"
            assert isinstance(result[key], list), f"{key} should be a list"

        assert len(result["variables"]) > 0, "No variables detected in MSO"


def test_generated_config_file_structure(config_setup, nni):
    """
    Test that calls generate_config_file(mso : list[np.ndarray], name: str = "", path: str = "")
    and creates the config-file.

    The config-file should be compatible with neural_residual toolbox and used
    directly to initialize a nn automatically.

    The config file data should also be correctly set up according to the config creators input data.

    Solves: Issue #27
    """
    # Setup for the test path and file
    test_dir, file_name, full_path = config_setup
    mso = _get_complete_mso()

    # Call for the config file generator
    nni.generate_config_file(mso=mso, name=file_name, path=test_dir)

    # Check that the file was generated correctly and with the correct content
    accepted_file_formats = [".yaml", ".yml"]
    assert os.path.isfile(
        full_path
    ), "The config file was not created and thus couldn't be found!"

    # Split the suffix from the file and format to lower
    assert (
        os.path.splitext(full_path)[-1].lower() in accepted_file_formats
    ), "The file format is not correct! Should be '.yml' or '.yaml'"

    # Check that the file contains the correct information
    with open(full_path, "r") as config_file:
        config_data = yaml.safe_load(config_file)

        # First see that the file is not empty
        assert config_data is not None, "Config file was empty!"

        # Then check that the file contains the rest of the needed information
        required_keys = [
            "description",
            "dataset",
            "signals",
            "zeroed_signals",
            "dynamic",
            "predictors",
        ]
        for key in required_keys:
            # NOTE TO SELF - Maybe zeroed signals is not required here and could be treated as default null if missing?
            assert key in config_data, f"The created config file missed the {key} tag"

        # Check that there is some data for the required ones
        needs_data = ["description", "dataset", "signals", "dynamic", "predictors"]
        for key in needs_data:
            assert config_data[key], f"Key '{key}' exists but lacks data!"

        # Confirm correct data type
        string_types = ["description", "dataset"]
        dictionary_types = ["signals", "dynamic", "predictors"]
        for key in string_types:
            assert isinstance(
                config_data[key], dict
            ), f"The {key} was not in correct format (should have been of type string)"
        for key in dictionary_types:
            assert isinstance(
                config_data[key], dict
            ), f"The {key} was not in correct format (should have been of type dictionary)"

        # Check that the signals variables are all in the form y_...
        for signal in config_data["signals"]:
            assert signal.startswith(
                "y"
            ), f"Signal {signal} is not a known sensor or input"

        # Confirm that all predictors are defined in the signals model aswell
        for predictor in config_data["predictors"]:
            assert predictor in [
                signal.key() for signal in config_data["signals"]
            ], f"Predictor: {predictor} is an unknown variable!"


def test_generated_config_file_data(config_setup, nni):
    """
    Test that the data of the VEP4Engine model is correct for the greybox.yml config file.

    Solves: Issue #27
    """

    # Setup for the test path and file
    test_dir, file_name, full_path = config_setup
    mso = _get_complete_mso()

    # Call for the config file generator
    nni.generate_config_file(mso=mso, name=file_name, path=test_dir)

    # Check that the file contains the correct data
    with open(full_path, "r") as config_file:
        config_data = yaml.safe_load(config_file)

        # Test that signals have correct data
        correct_signals = {
            "intercooler_pressure": "y_p_ic",
            "intercooler_temperature": "y_T_ic",
            "intake_manifold_pressure": "y_p_im",
            "air_mass_flow": "y_W_af",
            "engine_speed": "y_omega_e",
            "throttle_position": "y_alpha_th",
            "wastegate_position": "y_u_wg",
            "injected_fuel_mass": "y_wfc",
            "ambient_temperature": "y_T_amb",
            "ambient_pressure": "y_p_amb",
        }
        correct_dynamic = {
            "m_t": {
                "states": ["m_t"],
                "inputs": ["y_W_af", "y_alpha_th", "y_omega_e"],
            },
            "T_t": {
                "states": ["T_t", "m_t"],
                "inputs": ["y_T_amb", "y_p_amb"],
            },
        }
        correct_predictors = {
            "y_p_ic": {
                "states": ["m_t"],
                "inputs": [],
            },
            "y_T_ic": {
                "states": ["T_t", "m_t"],
                "inputs": ["y_p_amb", "y_u_wg"],
            },
        }

        # Test that the data is correct
        signals = config_data.get("signals")
        dynamic = config_data.get("dynamic")
        predictors = config_data.get("predictors")

        for key, value in correct_signals.items():
            assert key in signals, "Expected signals key is wrong"
            assert signals[key] == value, "Expected signals value is wrong"

        for key, value in correct_dynamic.items():
            assert key in dynamic, "Expected dynamic key is wrong"
            assert dynamic[key] == value, "Expected dynamic value is wrong"

        for key, value in correct_predictors.items():
            assert key in predictors, "Expected predictors key is wrong"
            assert predictors[key] == value, "Expected predictors value is wrong"

        assert config_data["zeroed_signals"] is None, "zeroed_signals should be null"


def test_time(config_setup, nni):
    """
    Tests that the system generates config-file within 10s
    At most 50 eq for MSO. Can fail in CI pipeline.

    Solves: Issue #31
    """
    mso = _get_complete_mso()

    assert len(mso) <= 50, "Mso has more than 50 equations"

    # Setup for the test path and file
    test_dir, file_name, _ = config_setup

    runs = 5
    times = []

    for _ in range(runs):
        start = time.perf_counter()
        nni.generate_config_file(mso, file_name, test_dir)
        end = time.perf_counter()
        times.append(end - start)

    avg_time = sum(times) / len(times)

    assert avg_time <= 10, f"Average time was above 10 seconds: {avg_time}"
