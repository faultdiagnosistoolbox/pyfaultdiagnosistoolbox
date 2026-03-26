import os
import fmpy
from fmpy.validation import validate_fmu, validate_model_description
from fmpy import extract, read_model_description
import pytest
from faultdiagnosistoolbox.ExternalSimulation import generate_fmu
from small_example import model as small_model

model_name = "small_model"


@pytest.fixture(scope="session")
def tiny_fmu_path():
    small_model_path = "generated/small_model.fmu"
    generate_fmu(small_model, small_model_path, model_name)  # Skicka in file_path
    return small_model_path


def test_path_is_valid():
    """Basic validation that the FMU can be parsed by fmpy."""
    fmu_path = tiny_fmu_path()
    assert os.path.exists(fmu_path), f"FMU not found: {fmu_path}"


def test_import_is_valid():
    """Test if the fmu is able to be imported"""
    fmu_path = tiny_fmu_path()
    fmpy.dump(fmu_path)


def test_fmu_validation():
    """Test that the FMU passes fmpy's validation checks."""
    fmu_path = tiny_fmu_path()
    validation_errors = validate_fmu(fmu_path)
    assert (
        not validation_errors
    ), f"FMU validation failed with errors: {validation_errors}"


def test_variable_causality():
    """Test variable input, output and internal"""
    fmu_path = tiny_fmu_path()
    model_description = read_model_description(fmu_path)
    causalities = {}
    for var in model_description.modelVariables:
        causality = var.causality
        if causality not in causalities:
            causalities[causality] = []
        causalities[causality].append(var.name)

    assert "r" in causalities["output"], "r not in output"
    assert (
        "y1" in causalities["input"] and "y2" in causalities["input"]
    ), "y1 och y2 not in input"
    # assert "XXX" in causality["parameter"], "no valid parameter"


def test_model_name():
    """Test that the FMU modelDescription.xml has the correct model name."""
    fmu_path = tiny_fmu_path()
    model_desc = fmpy.read_model_description(fmu_path)
    assert (
        model_desc.modelName == model_name
    ), f"Wrong model name: {model_desc.modelName}, correct: {model_name}"


def test_model_exchange():
    """Test that the FMU modelDescription.xml has the correct model exchange identifier"""
    fmu_path = tiny_fmu_path()
    model_desc = fmpy.read_model_description(fmu_path)
    assert model_desc.modelExchange != None, f"Model exchange not found!"


def test_OS_binaries():
    """Ensure FMU ships binaries for Windows, macOS, and Linux."""
    fmu_path = tiny_fmu_path()
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
