from pathlib import Path

from faultdiagnosistoolbox.ExternalSimulation.generatePython import (
    generate_fmu_python_file,
)
from faultdiagnosistoolbox.ExternalSimulation.metadata_helpers import (
    extract_fmu_variable_metadata,
)
from electricmotor_model import model


def _get_Gamma_and_res_eq():
    mtes = model.MTES()
    red_eq = mtes[0][0]
    M0 = [e for e in mtes[0] if e != red_eq]
    Gamma = model.Matching(M0)
    return Gamma, red_eq


def test_generate_python(tmp_path):
    Gamma, red_eq = _get_Gamma_and_res_eq()
    metadata = extract_fmu_variable_metadata(model, Gamma=Gamma, res_eq=red_eq)
    model_path = "tests/test_external_simulation/electricmotor_model.py"
    model_name = "test_model"
    output_path = tmp_path / "test_model_fmu.py"

    generated_path = generate_fmu_python_file(
        metadata,
        model_path,
        model_name,
        res_eq=red_eq,
        output_path=output_path,
        residual_name="r1",
    )

    assert generated_path == output_path
    assert output_path.exists()


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as temp_dir:
        test_generate_python(Path(temp_dir))
    print("Test passed.")
