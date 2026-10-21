"""File for all neural network integration functions, integration for the Neural Residual Toolbox"""

import numpy as np
import os
from faultdiagnosistoolbox import DiagnosisModel, Matching
import faultdiagnosistoolbox.NeuralNetworkIntegration.helpers as h


class NeuralNetworkIntegration:
    """
    Contains the logic and creator of config files in .yml format
    that is used for the neuralresidual toolbox.
    """

    def __init__(self, model):
        """
        Init function for a NeuralNetworkIntegration object

        Parameters:
        ----------
        model (DiagnosisModel) : Symbolic diagnosis model

        """

        if model is None:
            raise ValueError("Model cannot be None")

        self.model: DiagnosisModel = model
        self.equations = []
        self.variables = []

    def get_int_causality_model(self, matching, res_eq):
        """
        Build a state-space–like model under integer causality assumptions. The
         function analyzes the structural X and Z matrices together with a
         matching to identify state, derivative, and output equations. Variable
         dependencies are extracted and stored in a structured representation.

        Parameters:
        ----------
        model (DiagnosisModel)  : model containing X, Z, and var names.
        matching (Matching)     : matching info defining equation–variable pairs.
        res_eq (int)            : index of the residual equation.

        Returns:
        list : a structured description of states, derivatives, and outputs
         with their variable dependencies.
        """

        self.analyze_matching(matching)  # Get variables and equations matching
        self.generate_sub_matrices()  # Generate submatrices x0 and z0 for the reduced equation

        # Skip residuals that are differential constraints not equations
        red_eq_text = self.model.syme[res_eq]
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
        Starting from a specified index, this function traverses the structural
        X and Z matrices to find all state variables and sensor dependencies that
        influence the selected variable.

        Parameters:
        ----------
        idx (int)           : index of the variable/equation to analyze.
        u_idx (list | None) : indices of all variables visited during the
         recursive search, is by default set to None

        Returns:
        tuple (list, list, list) : containg state_idx (indices of state variables
         affecting the selected variable), sensor_idx (indices of sensor/input
         variables affecting the selected variable) and used_idx (indices of all
         variables visited during the recursive search)
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
        Analyzes the given matchings and saves the idexes of the equations (match
         rows) in self.equations and variables (match columns) in self.variables

        Parameters:
        ----------
        matching (Matching) : matching info defining equation–variable pairs.

        Returns:
        ----------
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
         variables and the models equations) in self.x0 and model.Z (the relation
         between known variables and the models equations) in self.z0

        Returns:
        ----------
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

        Parameters:
        ----------
        state (list[int])  : variables representing the state of the model
        sensor (list[int]) : variables of the signals that affect the model
        var (list[int])    : the variable(s) of the output(s)

        Returns:
        ----------
        tuple (list, list, list) : containing x_input (names
         of the x input variables), x_output (names of the outputvariables) and
         z_input (names of the z input variables)
        """

        x_input = [self.model.x[self.variables[v]] for v in state]
        x_output = [self.model.x[self.variables[var]]]
        z_input = [self.model.z[v] for v in sensor]

        return x_input, x_output, z_input

    def generate_config_file(
        self,
        gamma: Matching,
        path: str,
        res_eq: int,
        params: h.GenConfigParams,
        file_name: str,
        type: str,
    ):
        """
        Generates a configuration file in YAML format

        Parameters:
        ----------
        path (str)               : directory path where the file should be saved.
        gamma (Matching)         : Matching for the residual generator
        params (GenConfigParams) : en GenConfigParams instans
        res_eq (int)             : index to equation to use as residual equation
        file_name (str)          : base name for the output file
        type (str)               : structure type, either 'Latent' or 'Greybox'

        Returns:
        ----------
        str : the path to the generated YAML file
        """
        self.validate_config_input(path, gamma, params, res_eq, file_name, type)

        ssmodel = self.get_int_causality_model(gamma.matching, res_eq)
        if not ssmodel:
            raise RuntimeError(
                "Could not generate config file for given data due to empty causility model"
            )

        # Check if more than one out
        num_out = sum(1 for elem in ssmodel if elem.elem_type.lower() == "out")
        if num_out > 2:
            raise RuntimeError(
                "Could not generate config file for given data, to many out states in causility model"
            )

        # Preprocess filename
        if file_name.lower().endswith(".yml"):
            file_name = file_name[:-4]
        elif file_name.lower().endswith(".yaml"):
            file_name = file_name[:-5]

        # Contents of the file
        data = {
            "description": h.QuotedString(
                f"Generated {type.lower()} configuration file: {file_name} for model: {self.model.name}"
            ),
            "dataset": self.model.name,
            "signals": h.FlowDict(
                {
                    name: h.QuotedString(f"{var}")
                    for name, var in zip(params.signals, params.y_vars)
                }
            ),
            "zeroed_signals": "null",
        }
        # Map x variables with y variables
        x_to_y = {x.lower(): y for x, y in zip(params.x_vars, params.y_vars)}

        # Latent model type
        if type.lower() == "latent":
            all_zin = sorted(list(set(z for eq in ssmodel for z in eq.zin)))
            num_states = sum(1 for eq in ssmodel if eq.elem_type.lower() == "state")

            data["dynamic"] = h.FlowDict(
                {
                    "latent": h.FlowDict(
                        {
                            "states": ["latent"],
                            "num_latents": num_states if num_states > 0 else 1,
                            "inputs": all_zin,
                        }
                    )
                }
            )

            data["predictors"] = h.FlowDict(
                {
                    y_var: h.FlowDict(
                        {
                            "states": ["latent"],
                            "use_latent": True,
                            "inputs": list(z for z in eq.zin),
                        }
                    )
                    for eq in ssmodel
                    if eq.elem_type.lower() == "out"
                    for y_var in [x_to_y.get(eq.xout[0].lower())]
                    if y_var is not None
                }
            )

        # Greybox model type
        elif type.lower() == "greybox":
            data["dynamic"] = h.FlowDict(
                {
                    eq.xout_int[0]: h.FlowDict(
                        {
                            "states": list(x for x in eq.xin),
                            "inputs": list(z for z in eq.zin),
                        }
                    )
                    for eq in ssmodel
                    if eq.elem_type.lower() == "state"
                }
            )

            data["predictors"] = h.FlowDict(
                {
                    y_var: h.FlowDict(
                        {
                            "states": list(x for x in eq.xin),
                            "inputs": list(z for z in eq.zin),
                        }
                    )
                    for eq in ssmodel
                    if eq.elem_type.lower() == "out"
                    for y_var in [x_to_y.get(eq.xout[0].lower())]
                    if y_var is not None
                }
            )

        full_path = os.path.join(path, f"{file_name}_{type.lower()}.yml")
        if path and not os.path.exists(path):
            os.makedirs(path)

        if os.path.isfile(full_path):
            raise ValueError("File already exists!")

        # Create the path and file with the content
        with open(full_path, "w") as file:
            # Write to file formating and line breaks
            file.write("# General Description\n")
            file.write(f"description: {data['description']}\n\n")
            file.write(f"dataset: {data['dataset']}\n\n")
            file.write(f"signals: {data['signals']}\n\n")
            file.write("# Signals to be zeroed\n")
            file.write(f"zeroed_signals: {data['zeroed_signals']}\n\n")
            file.write("# Dynamic Equations\n")
            file.write(f"dynamic: {data['dynamic']}\n\n")
            file.write("# Predictors\n")
            file.write(f"predictors: {data['predictors']}")

        return full_path

    def validate_config_input(self, path, gamma, params, res_eq, file_name, type):
        """Validates the input given to generate_config_file"""

        required = ["signals", "y_vars", "x_vars"]
        for attr in required:
            if not hasattr(params, attr):
                raise AttributeError(f"params missing attribute '{attr}'")

        # If params.signals is empty set
        params.signals = params.signals if params.signals else params.y_vars

        if len(params.x_vars) != len(params.y_vars):
            raise ValueError("params must have the same amount of x_vars and y_vars")
        elif len(params.signals) != len(params.y_vars):
            raise ValueError("params must have the same amount of signals and y_vars")
        elif not isinstance(type, str):
            raise ValueError("type must be a string")
        elif type.lower() not in ["latent", "greybox"]:
            raise ValueError("model_type must be 'latent' or 'greybox'")
        elif gamma is None:
            raise ValueError("gamma cannot be None")
        elif not isinstance(path, str):
            raise ValueError("path must be a string")
        elif not isinstance(file_name, str):
            raise ValueError("file_name must be a string")
        elif not isinstance(res_eq, int):
            raise TypeError("res_eq must be an integer")
        elif res_eq < 0 or res_eq >= len(self.model.syme):
            raise IndexError("res_eq out of range")
