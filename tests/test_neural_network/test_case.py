import numpy as np

from faultdiagnosistoolbox.NeuralNetworkIntegration import NeuralNetworkIntegration
from faultdiagnosistoolbox.NeuralNetworkIntegration.helpers import GenConfigParams
from test_engine_model import LiU_ICE_model

if __name__ == "__main__":
    # Inputs
    signals = [
        "intercooler_pressure",
        "intercooler_temperature",
        "intake_manifold_pressure",
        "air_mass_flow",
        "engine_speed",
        "throttle_position",
        "wastegate_position",
        "injected_fuel_mass",
        "ambient_temperature",
        "ambient_pressure",
    ]

    y_variables = [
        "y_p_ic",
        "y_T_ic",
        "y_p_im",
        "y_W_af",
        "y_omega_e",
        "y_alpha_th",
        "y_u_wg",
        "y_wfc",
        "y_T_amb",
        "y_p_amb",
    ]

    x_variables = [
        "p_ic",
        "T_ic",
        "p_im",
        "W_af",
        "omega_e",
        "alpha_th",
        "u_wg",
        "wfc",
        "T_amb",
        "p_amb",
    ]
    model = LiU_ICE_model
    nni = NeuralNetworkIntegration(model)
    print("Mso type: ", type(model.MSO()))
    m = model.MSO()

    mso = model.MSO()[584]
    print("Mso 585: ", mso)

    params = GenConfigParams(y_variables, x_variables, signals)
    gamma = model.Matching(np.setdiff1d(mso, mso[3]))

    generated_file = nni.generate_config_file(
        gamma=gamma,
        path="",
        res_eq=mso[3],
        params=params,
        file_name="config_585",
        type="greybox",
    )

    print(f"YAML-fil skapad: {generated_file}")
