"""Module for generating FMUs from diagnosis models and MSO selections."""

from .metadata_helpers import (
    extract_fmu_variable_metadata,
)
from pathlib import Path
import subprocess
import shutil
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
        generate_fmu_python_file(
            metadata=self.metadata,
            model_path=self.model_path,
            model_name=Path(self.model_path).stem,
            res_eq=self.res_eq,
        )

        # Build FMU using pythonFMU
        subprocess.run(
            ["pythonfmu", "build", "-f", "temp.py", "src/faultdiagnosistoolbox"],
            check=True,
        )

        # Remove temporary python file
        Path("temp.py").unlink()

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
