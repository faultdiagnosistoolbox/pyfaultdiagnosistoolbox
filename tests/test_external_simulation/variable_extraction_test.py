from small_example import model


from faultdiagnosistoolbox.ExternalSimulation.metadata_helpers import (
    extract_fmu_variable_metadata,
)


def _get_Gamma_and_res_eq():
    mtes = model.MTES()
    res_eq = mtes[0][0]
    M0 = [e for e in mtes[0] if e != res_eq]
    Gamma = model.Matching(M0)
    return Gamma, res_eq


def test_extract_fmu_variable_metadata_for_remove_later_selection():
    Gamma, res_eq = _get_Gamma_and_res_eq()
    metadata = extract_fmu_variable_metadata(model, Gamma=Gamma, res_eq=res_eq)

    assert "parameters" in metadata, "Metadata missing 'parameters' key"
    assert "signals" in metadata, "Metadata missing 'signals' key"
    assert set(metadata["parameters"]).issubset(
        set(model.parameters)
    ), "Metadata parameters not subset of model parameters"
    assert set(metadata["signals"]).issubset(
        set(model.z)
    ), "Metadata signals not subset of model signals"


def test_extract_fmu_variable_metadata_has_deterministic_variable_ordering():
    Gamma, res_eq = _get_Gamma_and_res_eq()
    metadata1 = extract_fmu_variable_metadata(model, Gamma=Gamma, res_eq=res_eq)
    metadata2 = extract_fmu_variable_metadata(model, Gamma=Gamma, res_eq=res_eq)

    assert metadata1 == metadata2, "Variable ordering is not deterministic"


def test_residual_variables():
    Gamma, res_eq = _get_Gamma_and_res_eq()
    metadata = extract_fmu_variable_metadata(model, Gamma=Gamma, res_eq=res_eq)

    # Check that the residual variable is included in the signals
    assert "y1" in metadata["signals"], "Residual variable 'y1' not included in signals"
    assert "y2" in metadata["signals"], "Residual variable 'y2' not included in signals"


if __name__ == "__main__":
    test_extract_fmu_variable_metadata_for_remove_later_selection()
    test_extract_fmu_variable_metadata_has_deterministic_variable_ordering()
    test_residual_variables()
    print("All tests passed.")
