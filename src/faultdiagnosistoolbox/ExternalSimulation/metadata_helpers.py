import numpy as np


def extract_fmu_variable_metadata(model, Gamma, res_eq):
    """Extract FMU-oriented variable metadata for residual generator.

    Parameters
    ----------
    model : DiagnosisModel
        Symbolic diagnosis model.
    Gamma : Gamma
        MSO selection for the model.
    res_eq : int
        Index of the equation to be used as residual in the generated FMU.

    Returns
    -------
    Object
        Metadata for the given residual generator, including lists of parameters
        and signals used by the residual equation.
    """
    if model.modelType != "Symbolic":
        raise ValueError("extract_fmu_variable_metadata requires a Symbolic model")

    # Local import avoids pulling CodeGeneration at module import time.
    from faultdiagnosistoolbox.CodeGeneration import UsedVars

    res_gen_eq = _residual_equation_set(Gamma, res_eq)

    return {
        "parameters": UsedVars(model.syme[res_gen_eq], model.parameters),
        "signals": UsedVars(model.syme[res_gen_eq], model.z),
    }


def _residual_equation_set(gamma, res_eq):
    """Build the equation set used by a sequential residual generator."""
    res_gen_eqs = np.array([], dtype=np.int64)
    for hall_component in gamma.matching:
        res_gen_eqs = np.concatenate((res_gen_eqs, hall_component.row))
    return np.concatenate((res_gen_eqs, [res_eq]))
