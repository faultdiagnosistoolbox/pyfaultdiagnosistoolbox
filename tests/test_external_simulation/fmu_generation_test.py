from faultdiagnosistoolbox.ExternalSimulation.metadata_helpers import (
    extract_fmu_variable_metadata,
)
from electricmotor_model import model
from faultdiagnosistoolbox.ExternalSimulation import FMUGenerator


def _get_Gamma_and_res_eq():
    mtes = model.MTES()
    red_eq = mtes[0][0]
    M0 = [e for e in mtes[0] if e != red_eq]
    Gamma = model.Matching(M0)
    return Gamma, red_eq


def test_generated_fmu():
    Gamma, red_eq = _get_Gamma_and_res_eq()
    model_path = "tests/test_external_simulation/electricmotor_model.py"
    generator = FMUGenerator(
        model=model,
        model_path=model_path,
        Gamma=Gamma,
        res_eq=red_eq,
    )
    generator.generate_fmu()


if __name__ == "__main__":
    test_generated_fmu()
    print("FMU generation test passed.")
