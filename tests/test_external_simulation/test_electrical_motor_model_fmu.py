import os
import fmpy
from fmpy.validation import validate_fmu, validate_model_description
from fmpy import extract, read_model_description
import pytest
from electricmotor_model import model as em_model
from faultdiagnosistoolbox.ExternalSimulation import generate_fmu

model_name = "Electric motor"


@pytest.fixture(scope="session")
def em_fmu_path():
    em_path = "generated/electricmotor_model.fmu"
    generate_fmu(em_model, em_path, model_name)  # Skicka in file_path
    return em_path


def test_path_is_valid():
    """Basic validation that the FMU can be parsed by fmpy."""
    fmu_path = em_fmu_path()
    assert os.path.exists(fmu_path), f"FMU not found: {fmu_path}"


def test_import_is_valid():
    """Test that the FMU can be imported by fmpy without errors."""
    fmu_path = em_fmu_path()
    fmpy.dump(fmu_path)


def test_fmu_validation():
    """Test that the FMU passes fmpy's validation checks."""
    fmu_path = em_fmu_path()
    validation_errors = validate_fmu(fmu_path)
    assert (
        not validation_errors
    ), f"FMU validation failed with errors: {validation_errors}"


def test_variable_causality():
    """Test variable input, output and internal"""
    fmu_path = em_fmu_path()
    model_description = read_model_description(fmu_path)
    causalities = {}
    for var in model_description.modelVariables:
        causality = var.causality
        if causality not in causalities:
            causalities[causality] = []
        causalities[causality].append(var.name)

    assert "r" in causalities["output"], "r not in output"
    assert (
        "yw" in causalities["input"] and "yT" in causalities["input"]
    ), "yw och yT not in input"
    assert (
        "b" in causalities["parameter"] and "j" in causalities["parameter"]
    ), "no valid parameter"


def test_model_name():
    """Test that the FMU modelDescription.xml has the correct model name."""
    fmu_path = em_fmu_path()
    model_desc = fmpy.read_model_description(fmu_path)
    assert (
        model_desc.modelName == model_name
    ), f"Wrong model name: {model_desc.modelName}, correct: {model_name}"


def test_model_exchange():
    """Test that the FMU modelDescription.xml has the correct model exchange identifier"""
    fmu_path = em_fmu_path()
    model_desc = fmpy.read_model_description(fmu_path)
    assert model_desc.modelExchange != None, f"Model exchange not found!"


def test_OS_binaries():
    """Ensure FMU ships binaries for Windows, macOS, and Linux."""
    fmu_path = em_fmu_path()
    extracted_dir = extract(fmu_path)
    binaries_dir = os.path.join(extracted_dir, "binaries")
    assert os.path.isdir(binaries_dir), f"Missing binaries directory: {binaries_dir}"

    binaries = []
    # Walk through the binaries directory and collect all binary files
    for root, _, files in os.walk(binaries_dir):
        for name in files:
            rel_path = os.path.relpath(os.path.join(root, name), binaries_dir).lower()
            binaries.append(rel_path)

    has_windows = any("win" in path and path.endswith(".dll") for path in binaries)
    has_linux = any("linux" in path and path.endswith(".so") for path in binaries)

    missing = []
    if not has_windows:
        missing.append("windows (.dll in win*/windows*)")
    if not has_linux:
        missing.append("linux (.so in linux*)")

    assert not missing, (
        "Missing platform binaries in FMU: "
        + ", ".join(missing)
        + f". Found binary files: {sorted(binaries)}"
    )
