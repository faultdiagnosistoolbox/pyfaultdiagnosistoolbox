"""File for all neural network integration functions, integration for the Neural Residual Toolbox"""

import numpy as np
import os
import yaml
from typing import Dict, Any
from faultdiagnosistoolbox import DiagnosisModel
import faultdiagnosistoolbox.NeuralNetworkIntegration.helpers as h

# Added formats to yaml
yaml.add_representer(h.FlowDict, h.flow_dict_representer)
yaml.add_representer(h.QuotedString, h.quoted_representer)


class NeuralNetworkIntegration:
    """
    Contains the logic and creator of config files in .yml format
    that is used for the neuralresidual toolbox
    """

    def __init__(self, model=None):
        self.model: DiagnosisModel = model
        self.equations = []
        self.variables = []

    def analyze_mso(self, mso: np.ndarray) -> Dict[str, Any]:
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

    def get_int_causality_model(self, model, matching, red_eq):
        """
        Build a state-space–like model under integer causality assumptions.
            The function analyzes the structural X and Z matrices together with
            a matching to identify state, derivative, and output equations.
            Variable dependencies are extracted and stored in a structured
            representation.

        Args:
            model (DiagnosisModel): model containing X, Z, and var names.
            matching (Matching): matching info defining equation–variable pairs.
            red_eq (int): index of the reduced residual equation.

        Returns:
            ssmodel (list[SsmodelElement]: structured description of states,
                derivatives, and outputs with their variable dependencies.
        """

        self.model = model
        self.analyze_matching(matching)  # Get variables and equations matching
        self.generate_sub_matrices()  # Generate submatrices x0 and z0 for the reduced equation

        # Skip residuals that are differential constraints not equations
        red_eq_text = self.model.syme[red_eq]
        if np.size(red_eq_text) != 1:
            return None

        ssmodel = []
        for v_idx in range(len(self.variables)):
            xo_int = []
            eqs_w_var = self.x0[:, v_idx]
            var_name = self.model.x[self.variables[v_idx]]

            if np.any(eqs_w_var == h.STATE_VALUE):
                elem_type = "state"
                row_mask = eqs_w_var == h.STATE_VALUE
                cols = np.where(self.x0[row_mask, :] == h.RELATED_VALUE)[1]
                xo_int = [self.model.x[self.variables[c]] for c in cols]

            elif var_name in str(red_eq_text):
                elem_type = "out"

            else:
                continue  # Skip variables that are not part of the model under the given causality

            # Get state and sensor dependencies for the current variable
            state_vars, sensor_vars, checked_idx = self.find_state_dependencies(v_idx)

            xin, xout, zin = self.get_io_variables(state_vars, sensor_vars, v_idx)

            ssmodel.append(
                h.SsmodelElement(
                    xout=xout,
                    xin=xin,
                    zin=zin,
                    elem_type=elem_type,
                    xout_int=xo_int,
                )
            )

        return ssmodel

    def find_state_dependencies(self, idx, u_idx=None):
        """
        Recursively determine state and sensor dependencies for a given variable.
            Starting from a specified index, this function traverses the
            structural X and Z matrices to find all state variables and sensor
            dependencies that influence the selected variable.

        Args:
            idx (int): index of the variable/equation to analyze.
            u_idx (list[int] or None): indices of all variables visited during
                the recursive search.

        Returns:
            state_idx (list[int]): indices of state variables affecting the
                selected variable.
            sensor_idx (list[int]): indices of sensor/input variables affecting
                the selected variable.
            used_idx (list[int]): indices of all variables visited during the
                recursive search.
        """

        sensor_idx = set()
        state_idx = set()
        used_idx = set() if u_idx is None else u_idx

        var_in_eq = [x for x, val in enumerate(self.x0[idx]) if val != 0]
        for var_idx in var_in_eq:
            eqs_w_var = self.x0[:, var_idx]

            if idx == var_idx:
                if self.x0[idx][var_idx] == h.RELATED_VALUE:
                    state_idx.add(var_idx)
            elif var_idx in used_idx:
                continue  # Skip already visited variables to prevent infinite loop
            elif np.any(eqs_w_var == h.RELATED_VALUE):
                state_idx.add(var_idx)
            else:
                used_idx.add(var_idx)
                new_sta__idx, new_sen_idx, new_uidx = self.find_state_dependencies(
                    var_idx, used_idx
                )
                state_idx.update(new_sta__idx)
                sensor_idx.update(new_sen_idx)
                used_idx.update(new_uidx)

            used_idx.add(var_idx)
            sensor_idx.update([z for z, val in enumerate(self.z0[idx]) if val != 0])

        return sorted(state_idx), sorted(sensor_idx), used_idx

    def analyze_matching(self, matching):
        """
        Analyzes the given matchings and saves the idexes of the equations
            (match rows) in self.equations and variables (match columns) in
            self.variables

        Args:
            matching (Matching): matching info defining equation–variable pairs.

        Returns:
            void
        """
        self.equations = []
        self.variables = []

        for match in matching:
            self.equations.extend([r for r in match.row])
            self.variables.extend([c for c in match.col])

        # Reverse orders to make plot upper triangular
        self.equations.reverse()
        self.variables.reverse()

    def generate_sub_matrices(self):
        """
        Generates the submatrices of model.X (the relation between unknown
            variables and the models equations) in self.x0 and model.Z (the
            relation between known variables and the models equations) in self.z0

        Returns:
            void
        """
        self.x0 = []
        self.z0 = []

        for r in self.equations:
            self.z0 = (
                np.vstack((self.z0, self.model.Z[r]))
                if len(self.z0) > 0
                else self.model.Z[r]
            )
            x = []
            for c in self.variables:
                x.append(self.model.X[r][c])

            self.x0 = np.vstack((self.x0, x)) if len(self.x0) > 0 else x

        return

    def get_io_variables(self, state, sensor, var):
        """
        Get the names of the variables identified as the output and intput
        signals of the model

        Args:
            state (list[int]): variables representing the state of the model
            sensor (list[int]): variables of the signals that affect the model
            var (list[int]): the variable(s) of the output(s)

        Returns:
            tuple (list[str], list[str], list[str]): containing the x_input
                (names of the x input variables), x_output (names of the output
                variables), z_input (names of the z input variables)
        """
        x_input = [self.model.x[self.variables[v]] for v in state]
        x_output = [self.model.x[self.variables[var]]]
        z_input = [self.model.z[v] for v in sensor]

        return x_input, x_output, z_input

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
            "signals": h.FlowDict(
                {signal: h.QuotedString(f"{signal}") for signal in analysis["signals"]}
            ),
            "zeroed_signals": None,
        }

        # Blackbox or Greybox according to user input
        if model_type.lower() == "latent":
            data["dynamic"] = h.FlowDict(
                {
                    "latent": h.FlowDict(
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
            data["predictors"] = h.FlowDict(
                {
                    pred: h.FlowDict(
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
            data["dynamic"] = h.FlowDict(
                {
                    "states": list(analysis["states"]),
                    "inputs": list(analysis["signals"]),
                }
            )
            data["predictors"] = h.FlowDict(
                {
                    pred: h.FlowDict(
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
