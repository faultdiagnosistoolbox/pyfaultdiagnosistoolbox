"""Module for generating FMUs from diagnosis models and MSO selections."""

from .metadata_helpers import (
    extract_fmu_variable_metadata,
)
from pathlib import Path
import importlib.util
import subprocess
import shutil
import sys
import tempfile
from .generatePython import generate_fmu_python_file


class FMUGenerator:
    """Class to generate an FMU from a diagnosis model and MSO selection."""

    def __init__(self, model, model_path, Gamma, res_eq):
        """Initialize the FMU generator with model and MSO selection.

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
            Defaults to 0 (first equation).
        """
        self.model = model
        self.model_path = model_path
        self.Gamma = Gamma
        self.res_eq = res_eq

        self.metadata = extract_fmu_variable_metadata(model, Gamma, res_eq)

    def generate_fmu(self):
        """Generate the FMU based on the initialized model and MSO selection."""

        # Generate Python file that is used for creating fmu
        model_name = Path(self.model_path).stem
        fmu_script = generate_fmu_python_file(
            metadata=self.metadata,
            model_path=self.model_path,
            model_name=model_name,
            res_eq=self.res_eq,
            output_path=f"{model_name}_fmu.py",
        )
        print("Temporary Python file for FMU generation created successfully.")
        with tempfile.TemporaryDirectory(prefix="fdt_fmu_") as temp_dir:
            toolbox_package = self._stage_toolbox_package(Path(temp_dir))

            # Build FMU using pythonFMU
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "pythonfmu",
                    "build",
                    "-f",
                    str(fmu_script),
                    str(toolbox_package),
                    self.model_path,
                ],
                check=True,
            )

        # Remove temporary python file
        fmu_script.unlink()

        # Create generated directory
        generated_dir = Path("generated")
        generated_dir.mkdir(exist_ok=True)

        # Move fmu-file into generated directory
        src = Path(f"{Path(self.model_path).stem}.fmu")
        output_fmu = generated_dir / src.name
        shutil.move(src, output_fmu)

        # Final check that the fmu-file has been generated correctly
        if not output_fmu.exists():
            raise FileNotFoundError("FMU build failed, .fmu file not found")

    def _stage_toolbox_package(self, temp_dir):
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
