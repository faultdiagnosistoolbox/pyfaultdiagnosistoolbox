from .metadata_helpers import (
    extract_fmu_variable_metadata,
)


class FMUGenerator:
    """Class to generate an FMU from a diagnosis model and MSO selection."""

    def __init__(self, model, output_path, model_path, Gamma, res_eq):
        """Initialize the FMU generator with model and MSO selection.

        Parameters
        ----------
        model : DiagnosisModel
            Symbolic diagnosis model.
        output_path : str
            Directory where generated FMU files are written.
        model_path : str
            Path to the model file.
        Gamma : Gamma
            MSO selection for the model.
        res_eq : int
            Index of the equation to be used as residual in the generated FMU.
            Defaults to 0 (first equation).
        """
        self.model = model
        self.output_path = output_path
        self.model_path = model_path
        self.Gamma = Gamma
        self.res_eq = res_eq

        self.metadata = extract_fmu_variable_metadata(model, Gamma, res_eq)

    def generate_fmu(self):
        """Generate the FMU based on the initialized model and MSO selection.
        Note: This is a placeholder implementation. Actual FMU generation logic
        would go here.
        """

        raise NotImplementedError(
            "FMU packaging is not implemented. Use extract_fmu_variable_metadata "
            "for variable extraction."
        )
