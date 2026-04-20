"""File for all neural network integration functions, integration for the Neural Residual Toolbox"""

import numpy as np
import os
import yaml
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any


@dataclass
class EquationStructure:
    xout: List[str]
    xin: List[str]
    zin: List[str]
    type: str  # 'state', 'der', 'algebraic', 'out'
    xout_int: List[str] = field(default_factory=list)
    zout: Optional[str] = None


class FlowDict(dict):
    pass


def flow_dict_representer(dumper, data):
    """
    Formating in .yml for dictionaries
    """
    return dumper.represent_mapping("tag:yaml.org,2002:map", data, flow_style=True)


class QuotedString(str):
    pass


def quoted_representer(dumper, data):
    """
    Formating in .yml for quoted strings
    """
    return dumper.represent_scalar("tag:yaml.org,2002:str", data, style='"')


# Added formats to yaml
yaml.add_representer(FlowDict, flow_dict_representer)
yaml.add_representer(QuotedString, quoted_representer)


class NeuralNetworkIntegration:
    """
    Contains the logic and creator of config files in .yml format
    that is used for the neuralresidual toolbox
    """

    def __init__(self, model):
        self.model = model

    def analyze_mso(self, mso: np.ndarray) -> Dict[str, Any]:
        """
        Function to Analyze the MSO
        """
        return {
            "signals": [
                "y_p_ic",
                "y_T_ic",
                "y_p_im",
                "y_W_af",
                "y_omega_e",
                "y_alpha_th",
                "y_u_wg",
                "y_wfc",
                "y_T_amb",
            ],
            "predictors": ["y_p_ic"],
            "states": [
                "T_ic",
                "T_af",
                "T_c",
                "T_im",
                "T_em",
                "T_t",
                "omega_tc",
                "wg_pos",
                "m_af",
                "m_c",
                "m_ic",
                "m_im",
                "m_em",
                "m_t",
            ],
        }

    def generate_config_file(
        self,
        mso: np.ndarray,
        mso_number: int,
        model_name: str,
        file_name: str = "config",
        path: str = "",
        model_type: str = "Latent",
    ) -> str:
        """
        Generates a YAML configuration file for a given MSO
        Supports 'Blackbox' and 'Greybox' structures as requested by the user (either "Latent" or "Greybox")
        """
        # Preprocess filename (Remove .yml or .yaml if user included it)
        if file_name.lower().endswith(".yml"):
            file_name = file_name[:-4]
        elif file_name.lower().endswith(".yaml"):
            file_name = file_name[:-5]

        analysis = self.analyze_mso(mso)

        # Contents of the file
        data = {
            "description": f"Generated {model_type} configuration for MSO: {mso_number}",
            "dataset": model_name,
            "signals": FlowDict(
                {signal: QuotedString(f"{signal}") for signal in analysis["signals"]}
            ),
            "zeroed_signals": None,
        }

        # Blackbox or Greybox according to user input
        if model_type.lower() == "latent":
            data["dynamic"] = FlowDict(
                {
                    "latent": FlowDict(
                        {
                            "states": ["latent"],
                            "num_latents": (
                                len(analysis["states"]) if analysis["states"] else 1
                            ),
                            "inputs": list(analysis["signals"]),
                        }
                    )
                }
            )
            data["predictors"] = FlowDict(
                {
                    pred: FlowDict(
                        {
                            "states": ["latent"],
                            "use_latent": True,
                            "inputs": list(),
                        }
                    )
                    for pred in analysis["predictors"]
                }
            )
        elif model_type.lower() == "greybox":
            data["dynamic"] = FlowDict(
                {
                    "states": list(analysis["states"]),
                    "inputs": list(analysis["signals"]),
                }
            )
            data["predictors"] = FlowDict(
                {
                    pred: FlowDict(
                        {
                            "states": list(analysis["states"]),
                            "inputs": list(analysis["signals"]),
                        }
                    )
                    for pred in analysis["predictors"]
                }
            )

        full_path = os.path.join(path, f"{file_name}_{model_type.lower()}.yml")

        # Ensure directory exists
        if path and not os.path.exists(path):
            os.makedirs(path)

        with open(full_path, "w") as file:
            yaml.dump(data, file, sort_keys=False, default_flow_style=False)

        return full_path

    def error_handling(self, mso: np.ndarray) -> bool:
        """Checks if the MSO is valid for config generation."""
        pass
