"""
This file has some unusual formatting due to the way it is generated.
It is not intended to be manually edited.
"""

import importlib.util
import sys
from pathlib import Path


def generate_fmu_python_file(
    metadata, model_path, model_name, res_eq, output_path, residual_name="r"
):
    """Generate a Python file that defines an FMU class for the given model and metadata."""
    code = []
    code.append(write_imports())
    code.append(write_helper_functions(model_name, model_path))
    code.append(write_class_definition(model_name))
    code.append(write_init(metadata, residual_name, res_eq))
    code.append(write_do_step(metadata, residual_name))
    result = "".join(code)

    output_path = Path(output_path)
    # Allow callers to pass paths in directories that do not exist yet.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        f.write(result)

    return output_path


def write_imports():
    """Write import statements for the generated Python code."""
    code = """from pythonfmu import Fmi2Slave, Fmi2Causality, Fmi2Variability, Fmi2Initial, Real
import tempfile
import importlib
import importlib.util
import sys
from pathlib import Path

from faultdiagnosistoolbox.ExternalSimulation.fmu_runtime import (
    build_runtime
)

"""
    return code


def write_helper_functions(model_name, model_path):
    """Write helper functions for the generated Python code, including model import."""
    model_path_repr = repr(str(model_path))
    model_filename_repr = repr(Path(model_path).name)

    code = f"""def _import_module_from_path(module_name, path):
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Module file not found at: {{path}}")

    spec = importlib.util.spec_from_file_location(module_name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not create spec for module at: {{path}}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _load_model_module():
    module_name = "{model_name}"
    original_model_path = Path({model_path_repr})
    packaged_model_path = Path(__file__).resolve().with_name({model_filename_repr})

    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError:
        pass

    for candidate_path in (packaged_model_path, original_model_path):
        try:
            return _import_module_from_path(module_name, candidate_path)
        except FileNotFoundError:
            continue

    raise ModuleNotFoundError(
        f"Could not import model module '{{module_name}}'. "
        f"Tried module import, packaged file '{{packaged_model_path}}', "
        f"and source path '{{original_model_path}}'."
    )


_model_module = _load_model_module()
model = getattr(_model_module, "model", _model_module)

if not hasattr(model, "MTES") or not hasattr(model, "Matching"):
    raise AttributeError(
        "Imported model does not provide MTES/Matching. "
        "Expected a DiagnosisModel object, usually exposed as module variable 'model'."
    )


def _get_Gamma(res_eq):
    mtes = model.MTES()
    M0 = [e for e in mtes[0] if e != res_eq]
    return model.Matching(M0)

"""
    return code


def import_module_from_path(module_name, path):
    """Import a Python module from a file path."""
    path = Path(path)

    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)

    sys.modules[module_name] = module
    spec.loader.exec_module(module)

    return module


def write_class_definition(name):
    """Write the class definition for the generated Python code."""
    code = f"""class {name}(Fmi2Slave):
    author = "Fault Diagnosis Toolbox"
    description = "Residual FMU"

"""
    return code


def write_init(metadata, residual_name, res_eq):
    """Write the __init__ method for the generated Python code,
    including variable registration and runtime initialization."""
    code = f"""    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        {write_variable_initialization(metadata, residual_name)}

        {write_variable_registration(metadata, residual_name)}

        {write_runtime_initialization(res_eq, residual_name)}

"""
    return code


def write_variable_initialization(metadata, residual_name):
    """Write code to initialize variables in the __init__ method of the generated Python code."""
    code = []
    code.append(f"self.{residual_name} = 0.0")
    for signal in metadata["signals"]:
        code.append(f"self.{signal} = 0.0")
    for parameter in metadata["parameters"]:
        code.append(f"self.{parameter} = 0.0")
    if not code:
        return "pass"
    return "\n        ".join(code)


def write_variable_registration(metadata, residual_name):
    """Write code to register variables in the __init__ method of the generated Python code."""
    code = []
    code.append(f"""self.register_variable(Real("{residual_name}",
                causality=Fmi2Causality.output,
                variability=Fmi2Variability.continuous,
                initial=Fmi2Initial.exact))""")
    for signal in metadata["signals"]:
        code.append(f"""self.register_variable(Real("{signal}",
                causality=Fmi2Causality.input))""")
    for parameter in metadata["parameters"]:
        code.append(f"""self.register_variable(Real("{parameter}",
                causality=Fmi2Causality.parameter,
                variability=Fmi2Variability.tunable))""")
    if not code:
        return "pass"
    return "\n        ".join(code)


def write_runtime_initialization(res_eq, residual_name):
    """Write code to initialize the residual runtime in the __init__ method of the generated Python code."""
    code = [
        'tmp_dir = tempfile.mkdtemp(prefix="fdt_res_")',
        f"gamma = _get_Gamma({res_eq})",
        "self.runtime = build_runtime(",
        "    model=model,",
        "    Gamma=gamma,",
        f"    res_eq={res_eq},",
        f'    residual_name="{residual_name}",',
        "    output_dir=tmp_dir,",
        ")",
    ]
    return "\n        ".join(code)


def write_do_step(metadata, residual_name):
    """Write the do_step method for the generated Python code, which computes the residual output."""
    signal_lines = [
        f'"{signal}": self.{signal},' for signal in metadata.get("signals", [])
    ]
    parameter_lines = [
        f'"{parameter}": self.{parameter},'
        for parameter in metadata.get("parameters", [])
    ]

    signal_dict = "\n                ".join(signal_lines)
    parameter_dict = "\n                ".join(parameter_lines)

    return f"""
    def do_step(self, current_time: float, step_size: float) -> bool:
        known_signals = {{
            {signal_dict}
        }}
        parameters = {{
            {parameter_dict}
        }}

        out = self.runtime.step(known_signals, parameters=parameters, ts=step_size)
        self.{residual_name} = float(out)

        return True
"""
