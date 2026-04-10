"""runtime helpers for residual evaluation in an FMU core loop.

This module shows how to:
- Generate Python residual evaluation code for a selected MSO set and residual equation.
- Load the generated code as a Python module.
- Provide a simple runtime interface to evaluate the residual for given signal and parameter values.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import os
import numpy as np

from .metadata_helpers import (
    extract_fmu_variable_metadata,
)


@contextmanager
def _working_directory(path: Path):
    """Temporarily switch working directory."""
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def _load_module(module_path: Path, unique_module_name: str):
    """Load a generated Python module from a file path."""
    spec = spec_from_file_location(unique_module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module from '{module_path}'")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def import_module_from_path(path):
    """Import a Python module from a file path."""
    module_name = Path(path).stem
    spec = spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load module from '{path}'")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@dataclass
class FMURuntime:
    """Container that evaluates generated residual functions.

    Parameters
        ----------
    model : DiagnosisModel
    Symbolic diagnosis model.

    residual_name: str
    Name of the residual function to evaluate, e.g. "r1".

    module: object
    Loaded Python module containing the generated residual function.

    metadata: object
    Metadata about the model variables used for residual evaluation, extracted by `extract_fmu_variable_metadata`.

    states: dict[str, float]
    Dictionary of state variable values for each residual.
    Keys are residual names, values are dicts of state variable names to current values.
    """

    model: object
    residual_name: str
    module: object
    metadata: object
    states: dict[str, float]

    def step(
        self,
        known_signals: dict[str, float],
        parameters: dict[str, float] | None = None,
        Ts: float = 0.1,
    ):
        """Compute one residual output for the given signals and parameters.

        Parameters
        ----------
        known_signals : dict[str, float]
            Signal values keyed by model known signal names (`model.z`).
        parameters : dict[str, float], optional
            Parameter values keyed by model parameter names.
        Ts : float, optional
            Sampling time used by generated residual functions.
        """
        missing_signals = [
            signal for signal in self.model.z if signal not in known_signals
        ]
        if missing_signals:
            raise KeyError(
                f"Missing known signals for residual evaluation: {missing_signals}"
            )

        params = {} if parameters is None else parameters
        required_parameters = self.metadata.get("parameters", [])
        missing_parameters = [
            parameter for parameter in required_parameters if parameter not in params
        ]
        if missing_parameters:
            raise KeyError(
                f"Missing parameters for residual evaluation: {missing_parameters}"
            )

        # Generated residual functions index z by positions in model.z.
        z = np.array([known_signals[signal] for signal in self.model.z], dtype=float)

        residual_function = getattr(self.module, self.residual_name)
        state = self.states.get(self.residual_name, {})
        residual_value, updated_state = residual_function(z, state, params, Ts)
        self.states[self.residual_name] = updated_state
        return float(residual_value)


def build_runtime(
    model,
    Gamma,
    res_eq,
    residual_name="r",
    output_dir: str | Path = ".",
    diffres: str = "int",
    initial_state_by_residual: dict[str, float] | None = None,
):
    """Build a residual runtime object for FMU-like step evaluation.

    Parameters
    ----------
    model : DiagnosisModel
        Symbolic model used to generate residual code.
    Gamma : Gamma
        MSO selection for the model.
    res_eq : int
        Index of the equation to be used as residual in the generated code.
    residual_name : str, optional
        Name for the generated residual output. Defaults to "r1".
    output_dir : str | Path, optional
        Directory where generated code is written. Defaults to current directory.
    differs : str, optional
        Type of residual to generate, either "int" or "ext". Defaults to "int".
    initial_state_by_residual : dict[str, float], optional
        Initial state values for each residual. Keys are residual names, values
        are dicts of state variable names to initial values. Defaults to None (no states).
    """

    metadata = extract_fmu_variable_metadata(model, Gamma, res_eq)

    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)

    with _working_directory(destination):
        model.SeqResGen(
            Gamma,
            res_eq,
            residual_name,
            diffres=diffres,
            language="Python",
            batch=False,
        )

    module = _load_module(
        destination / f"{residual_name}.py", f"_fdt_extsim_residual_{residual_name}"
    )

    initial_states = (
        initial_state_by_residual.get(residual_name, {})
        if initial_state_by_residual
        else {}
    )
    states = {residual_name: initial_states}

    return FMURuntime(
        model=model,
        residual_name=residual_name,
        module=module,
        metadata=metadata,
        states=states,
    )
