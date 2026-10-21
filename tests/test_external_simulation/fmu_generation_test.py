from pathlib import Path

from faultdiagnosistoolbox.ExternalSimulation.metadata_helpers import (
    extract_fmu_variable_metadata,
)
from electricmotor_model import model
from faultdiagnosistoolbox.ExternalSimulation import generate_fmu


def _get_Gamma_and_res_eq():
    mtes = model.MTES()
    red_eq = mtes[0][0]
    M0 = [e for e in mtes[0] if e != red_eq]
    Gamma = model.Matching(M0)
    return Gamma, red_eq


def test_generated_fmu():
    Gamma, red_eq = _get_Gamma_and_res_eq()
    model_path = "tests/test_external_simulation/electricmotor_model.py"
    fmu_path = generate_fmu(
        model=model,
        model_path=model_path,
        Gamma=Gamma,
        res_eq=red_eq,
        fmu_name="custom_electricmotor_model",
    )
    assert fmu_path == Path("generated/custom_electricmotor_model.fmu")
    assert fmu_path.exists()


if __name__ == "__main__":
    test_generated_fmu()
    print("FMU generation test passed.")
