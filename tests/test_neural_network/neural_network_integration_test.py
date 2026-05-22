import pytest
import numpy as np
import os
import yaml
import time
from test_engine_model import LiU_ICE_model
from faultdiagnosistoolbox.NeuralNetworkIntegration.NeuralNetworkIntegration import (
    NeuralNetworkIntegration,
)
from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import (
    GenConfigParams,
)


@pytest.fixture
def nni():
    """Initializes and returns an instance of NeuralNetworkIntegration"""
    model = LiU_ICE_model
    return NeuralNetworkIntegration(model)


@pytest.fixture
def model():
    model = LiU_ICE_model
    return model


@pytest.fixture
def config_setup(tmp_path):
    """Creates a temporary directory, and returns the directory, filename and path"""
    test_dir = str(tmp_path / _get_config_path())
    os.mkdir(test_dir)
    file_name = _get_config_name()
    model_type = _get_model_type()
    model_name = _get_model_name()
    mso_index = _get_mso_index()

    path = os.path.join(test_dir, file_name)
    full_path = path + f"_{model_type}.yml"

    yield test_dir, file_name, full_path, model_type, model_name, mso_index

    # Cleanup
    if os.path.exists(full_path):
        os.remove(full_path)


def _get_config_path():
    """Returns a test config folder path name"""
    config_path = "configs"
    return config_path


def _get_config_name():
    """Returns a test config file name"""
    name = "model_1_config"
    return name


def _get_model_type():
    """Returns a test model type"""
    model_type = "latent"
    return model_type


def _get_model_name():
    """Returns a test model name"""
    model_name = LiU_ICE_model.name
    return model_name


def _get_mso_index():
    """Returns a test mso index"""
    mso_index = 585
    return mso_index


def _get_complete_mso():
    """
    Return mso 585 (index 584) from a set of msos
    """
    return LiU_ICE_model.MSO()[584]


def _get_mso_585_inputs() -> tuple[list[str], list[str], list[str]]:
    """
    Return inputs for mso 585 in the LiU_ICE_model
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


def test_config_file_created_at_destination(config_setup, nni, model):
    """
    Tests that the file is generated at the correct given destination and with the correct filename.
    """
    # Setup for the test path and file
    test_dir, file_name, expected_full_path, model_type, model_name, mso_index = (
        config_setup
    )

    mso = _get_complete_mso()
    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    nni.generate_config_file(gamma, test_dir, mso[3], params, file_name, model_type)

    # Check that the file and path was created at correct locations and with correct names
    assert os.path.exists(
        expected_full_path
    ), f"The path was not created correctly at: {expected_full_path}!"
    assert os.path.isfile(
        expected_full_path
    ), f"Expected the file to be created at {expected_full_path}, but the file was not found!"


def test_generated_config_file_structure(config_setup, nni, model):
    """
    Test that calls generate_config_file(mso : list[np.ndarray], name: str = "", path: str = "")
    and creates the config-file.

    The config-file should be compatible with neural_residual toolbox and used
    directly to initialize a nn automatically.

    The config file data should also be correctly set up according to the config creators input data.

    """
    # Setup for the test path and file
    test_dir, file_name, full_path, model_type, model_name, mso_index = config_setup
    mso = _get_complete_mso()
    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    nni.generate_config_file(gamma, test_dir, mso[3], params, file_name, model_type)

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
                config_data[key], str
            ), f"The {key} was not in correct format (should have been of type string)"
        for key in dictionary_types:
            assert isinstance(
                config_data[key], dict
            ), f"The {key} was not in correct format (should have been of type dictionary)"

        # Check that the signals variables are all in the form y_...
        for signal in config_data["signals"].values():
            assert signal.startswith(
                "y"
            ), f"Signal {signal} is not a known sensor or input"

        # Confirm that all predictors are defined in the signals model aswell
        for predictor in config_data["predictors"]:
            assert (
                predictor in config_data["signals"].values()
            ), f"Predictor: {predictor} is an unknown variable!"


def test_generated_latent_data(config_setup, nni, model):
    """
    Verified that the generated Latent YAML matches the expected structure
    from the MATLAB reference for MSO 585.
    """
    test_dir, file_name, _, _, model_name, mso_number = config_setup
    model_type = "latent"
    full_path = os.path.join(test_dir, f"{file_name}_{model_type}.yml")

    mso = _get_complete_mso()

    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    nni.generate_config_file(gamma, test_dir, mso[3], params, file_name, model_type)

    with open(full_path, "r") as f:
        config_data = yaml.safe_load(f)

        # Verify Description
        assert "description" in config_data, "Description tag is missing"
        assert (
            config_data["description"]
            == f"Generated {model_type} configuration file: {file_name} for model: {model_name}"
        )
        assert config_data["dataset"] == model_name, "Dataset name mismatch"
        assert config_data["zeroed_signals"] is None, "zeroed_signals should be null"

        # Verify Dynamic
        assert "latent" in config_data["dynamic"]
        latent_cfg = config_data["dynamic"]["latent"]
        assert latent_cfg["states"] == ["latent"]
        assert latent_cfg["num_latents"] == 13

        # Check that all the inputs are correct
        expected_inputs = [
            "y_T_amb",
            "y_W_af",
            "y_alpha_th",
            "y_omega_e",
            "y_p_amb",
            "y_p_ic",
            "y_p_im",
            "y_u_wg",
            "y_wfc",
        ]
        assert sorted(latent_cfg["inputs"]) == sorted(expected_inputs)

        # Verify Predictor
        assert "y_T_ic" in config_data["predictors"]
        pred_cfg = config_data["predictors"]["y_T_ic"]
        assert pred_cfg["states"] == ["latent"]
        assert sorted(pred_cfg["inputs"]) == sorted(
            ["y_p_ic", "y_p_im", "y_alpha_th", "y_T_amb"]
        )


def test_generated_greybox_data(config_setup, nni, model):
    """
    Test that the generated Greybox YAML matches the expected structure
    and data from the VEP4Engine model for MSO 585.
    """
    # Setup for the test path and file
    test_dir, file_name, _, _, model_name, mso_number = config_setup
    model_type = "greybox"

    # Ensure we use the correct filename suffix for the greybox assertion
    full_path = os.path.join(test_dir, f"{file_name}_{model_type}.yml")

    mso = _get_complete_mso()
    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    nni.generate_config_file(gamma, test_dir, mso[3], params, file_name, model_type)

    # Check that the file exists before attempting to open it
    assert os.path.exists(full_path), f"Greybox file was not created at {full_path}"

    # Load and verify the content of the generated YAML file
    with open(full_path, "r") as config_file:
        config_data = yaml.safe_load(config_file)

        # Expected data mapping for signals
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

        signals = config_data.get("signals")
        dynamic = config_data.get("dynamic")
        predictors = config_data.get("predictors")

        # Verify Description
        assert "description" in config_data, "Description tag is missing"
        assert (
            config_data["description"]
            == f"Generated {model_type} configuration file: {file_name} for model: {model_name}"
        )
        assert config_data["dataset"] == model_name, "Dataset name mismatch"
        assert config_data["zeroed_signals"] is None, "zeroed_signals should be null"

        # Verify Signals Mapping
        for key, value in correct_signals.items():
            assert key in signals, f"Missing signal key: {key}"
            assert signals[key] == value, f"Incorrect signal mapping for {key}"

        # Verify Dynamic Equations (States and Inputs)
        # Check m_t
        assert "m_t" in dynamic, "m_t missing from dynamic equations"
        assert set(dynamic["m_t"]["states"]) == {"T_em", "T_t", "m_em", "m_t", "wg_pos"}
        assert dynamic["m_t"]["inputs"] == ["y_p_amb"]

        # Check omega_tc (Edge case: no inputs)
        assert "omega_tc" in dynamic, "omega_tc missing from dynamic equations"
        assert set(dynamic["omega_tc"]["states"]) == {
            "T_af",
            "T_c",
            "T_em",
            "T_t",
            "omega_tc",
            "m_em",
            "m_t",
            "m_af",
            "m_c",
        }
        assert dynamic["omega_tc"]["inputs"] == []

        # Check wg_pos
        assert "wg_pos" in dynamic
        assert dynamic["wg_pos"]["states"] == ["wg_pos"]
        assert dynamic["wg_pos"]["inputs"] == ["y_u_wg"]

        # Verify Predictors
        assert "y_T_ic" in predictors, "y_T_ic missing from predictors"
        pred_y_T_ic = predictors["y_T_ic"]
        assert set(pred_y_T_ic["states"]) == {"T_ic", "m_ic", "T_c", "m_c"}
        assert set(pred_y_T_ic["inputs"]) == {
            "y_p_ic",
            "y_p_im",
            "y_alpha_th",
            "y_T_amb",
        }


def test_time_greybox(config_setup, nni, model):
    """
    Tests that the system generates config-file within 10s for greybox model type,
    With an mso with at most 50 eq.
    """
    mso = _get_complete_mso()
    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    # Setup for the test path and file - for any MSO less than 50 equations
    if len(mso) <= 50:
        test_dir, file_name, _, _, model_name, mso_index = config_setup
        model_type = "greybox"

        runs = 5
        times = []

        for _ in range(runs):
            start = time.perf_counter()
            nni.generate_config_file(
                gamma, test_dir, mso[3], params, file_name, model_type
            )
            end = time.perf_counter()
            times.append(end - start)

        avg_time = sum(times) / len(times)

        assert avg_time <= 10, f"Average time was above 10 seconds: {avg_time}"


def test_time_latent(config_setup, nni, model):
    """
    Tests that the system generates config-file within 10s for latent model type,
    With an mso with at most 50 eq.
    """
    mso = _get_complete_mso()
    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    # Setup for the test path and file - for any MSO less than 50 equations
    if len(mso) <= 50:
        test_dir, file_name, _, model_type, model_name, mso_index = config_setup

        runs = 5
        times = []

        for _ in range(runs):
            start = time.perf_counter()
            nni.generate_config_file(
                gamma, test_dir, mso[3], params, file_name, model_type
            )
            end = time.perf_counter()
            times.append(end - start)

        avg_time = sum(times) / len(times)

        assert avg_time <= 10, f"Average time was above 10 seconds: {avg_time}"


def test_config_file_created_through_model(config_setup, model):
    """
    Tests that the file is generated at the correct given destination and with the correct filename.

    Solves: Issue #31
    """
    # Setup for the test path and file
    test_dir, file_name, expected_full_path, model_type, _, _ = config_setup

    mso = _get_complete_mso()
    params = GenConfigParams(
        _get_mso_585_inputs()[1], _get_mso_585_inputs()[2], _get_mso_585_inputs()[0]
    )
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))
    model.GenerateConfigFile(test_dir, gamma, mso[3], params, file_name, model_type)

    # Check that the file and path was created at correct locations and with correct names
    assert os.path.exists(
        expected_full_path
    ), f"The path was not created correctly at: {expected_full_path}!"
    assert os.path.isfile(
        expected_full_path
    ), f"Expected the file to be created at {expected_full_path}, but the file was not found!"
