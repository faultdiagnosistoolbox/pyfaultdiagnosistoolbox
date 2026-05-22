from pathlib import Path

from small_example import model


def _get_Gamma_and_res_eq():
    mtes = model.MTES()
    red_eq = mtes[0][0]
    M0 = [e for e in mtes[0] if e != red_eq]
    Gamma = model.Matching(M0)
    return Gamma, red_eq


def test_fmu_model_integration():
    model_path = "tests/test_external_simulation/small_example.py"
    Gamma, red_eq = _get_Gamma_and_res_eq()
    fmu_path = model.GenerateFMU(
        model_path, Gamma, red_eq, fmu_name="small_example_custom"
    )
    assert fmu_path == Path("generated/small_example_custom.fmu")
    assert fmu_path.exists()


if __name__ == "__main__":
    test_fmu_model_integration()
