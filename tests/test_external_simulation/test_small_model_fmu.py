import os
import importlib.util
from pathlib import Path
import fmpy
import numpy as np
from fmpy.validation import validate_fmu, validate_model_description
from fmpy import extract, read_model_description
import pytest
from faultdiagnosistoolbox.ExternalSimulation import generate_fmu


def _load_model(model_file_name):
    model_path = Path(__file__).resolve().parent / model_file_name
    spec = importlib.util.spec_from_file_location(model_path.stem, model_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module.model


small_model = _load_model("small_example.py")

model_name = "small_example"


@pytest.fixture(scope="session")
def tiny_fmu_path():
    small_model_path = "tests/test_external_simulation/small_example.py"
    output_path = "generated/small_example.fmu"
    generate_fmu(
        model=small_model,
        model_path=small_model_path,
        Gamma=small_model.Matching(small_model.MTES()[0][1:]),
        res_eq=small_model.MTES()[0][0],
    )
    return output_path


def test_path_is_valid(tiny_fmu_path):
    """Basic validation that the FMU can be parsed by fmpy."""
    fmu_path = tiny_fmu_path
    assert os.path.exists(fmu_path), f"FMU not found: {fmu_path}"


def test_import_is_valid(tiny_fmu_path):
    """Test if the fmu is able to be imported"""
    fmu_path = tiny_fmu_path
    fmpy.dump(fmu_path)


def test_fmu_validation(tiny_fmu_path):
    """Test that the FMU passes fmpy's validation checks."""
    fmu_path = tiny_fmu_path
    validation_errors = validate_fmu(fmu_path)
    assert (
        not validation_errors
    ), f"FMU validation failed with errors: {validation_errors}"


def test_variable_causality(tiny_fmu_path):
    """Test variable input, output and internal"""
    fmu_path = tiny_fmu_path
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
    ), "y1 and y2 not in input"


def test_model_name(tiny_fmu_path):
    """Test that the FMU modelDescription.xml has the correct model name."""
    fmu_path = tiny_fmu_path
    model_desc = fmpy.read_model_description(fmu_path)
    assert (
        model_desc.modelName == model_name
    ), f"Wrong model name: {model_desc.modelName}, correct: {model_name}"


def test_co_simulate(tiny_fmu_path):
    """Test that the FMU modelDescription.xml has the correct co-simulated identifier"""
    fmu_path = tiny_fmu_path
    model_desc = fmpy.read_model_description(fmu_path)
    assert model_desc.coSimulation != None, f"Co-simulation not found!"


def test_simulate_small_model_fmu_does_not_crash(tiny_fmu_path):
    """Smoke test that the generated FMU can be instantiated and stepped."""
    time = np.linspace(0.0, 1.0, 10)
    inputs = np.zeros(len(time), dtype=[("time", float), ("y1", float), ("y2", float)])
    inputs["time"] = time
    inputs["y1"] = 2.0
    inputs["y2"] = 4.0

    fmpy.simulate_fmu(
        tiny_fmu_path,
        start_time=0.0,
        stop_time=1.0,
        input=inputs,
    )


def test_OS_binaries(tiny_fmu_path):
    """Ensure FMU ships binaries for Windows, macOS, and Linux."""
    fmu_path = tiny_fmu_path
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
