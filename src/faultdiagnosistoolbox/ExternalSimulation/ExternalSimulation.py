"""Module for generating FMUs from diagnosis models and MSO selections."""

from pathlib import Path
import importlib.util
import shutil
import subprocess
import sys
import tempfile

from .generatePython import generate_fmu_python_file
from .metadata_helpers import extract_fmu_variable_metadata


def _stage_toolbox_package(temp_dir):
    """Stage the toolbox package with compiled extensions for FMU packaging."""
    source_package = Path("src/faultdiagnosistoolbox")
    staged_package = temp_dir / source_package.name
    shutil.copytree(source_package, staged_package)

    dmpermlib_spec = importlib.util.find_spec("faultdiagnosistoolbox.dmpermlib")
    if dmpermlib_spec is None or dmpermlib_spec.origin is None:
        raise ModuleNotFoundError(
            "Could not find compiled extension 'faultdiagnosistoolbox.dmpermlib'. "
            "Build or install faultdiagnosistoolbox before generating the FMU."
        )

    dmpermlib_path = Path(dmpermlib_spec.origin)
    shutil.copy2(dmpermlib_path, staged_package / dmpermlib_path.name)
    return staged_package


def generate_fmu(model, model_path, Gamma, res_eq, fmu_name=None):
    """Generate an FMU from a diagnosis model and MSO selection.

    Parameters
    ----------
    model : DiagnosisModel
        Symbolic diagnosis model.
    model_path : str
        Path to the model file.
    Gamma : Gamma
        MSO selection for the model.
    res_eq : int
        Index of the equation to be used as residual in the generated FMU.
    fmu_name : str, optional
        Output filename for the generated FMU. The `.fmu` extension is added
        automatically when omitted.
    """

    metadata = extract_fmu_variable_metadata(model, Gamma, res_eq)
    model_name = Path(model_path).stem

    # Generate temporary Python file that is used for creating fmu
    with tempfile.TemporaryDirectory(prefix="fdt_fmu_") as temp_dir:
        temp_path = Path(temp_dir)
        fmu_script = generate_fmu_python_file(
            metadata=metadata,
            model_path=model_path,
            model_name=model_name,
            res_eq=res_eq,
            output_path=temp_path / f"{model_name}_fmu.py",
        )
        print("Temporary Python file for FMU generation created successfully.")

        toolbox_package = _stage_toolbox_package(temp_path)

        subprocess.run(
            [
                sys.executable,
                "-m",
                "pythonfmu",
                "build",
                "-f",
                str(fmu_script),
                str(toolbox_package),
                model_path,
            ],
            check=True,
        )
    print("Temporary Python file removed after FMU generation.")

    # Create generated directory
    generated_dir = Path("generated")
    generated_dir.mkdir(exist_ok=True)

    # Move fmu-file into generated directory
    src = Path(f"{model_name}.fmu")
    output_name = Path(fmu_name).name if fmu_name else src.name
    if not output_name.endswith(".fmu"):
        output_name = f"{output_name}.fmu"
    output_fmu = generated_dir / output_name
    shutil.move(src, output_fmu)

    # Final check that the fmu-file has been generated correctly
    if not output_fmu.exists():
        raise FileNotFoundError("FMU build failed, .fmu file not found")

    print(f"FMU generated successfully: {output_fmu}")
    print(
        "Manual for running the FMU in MATLAB is available in "
        "docs/run_fmu_in_matlab.rst"
    )
    return output_fmu
