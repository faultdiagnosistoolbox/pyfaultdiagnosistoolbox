import numpy as np


def extract_fmu_variable_metadata(model, Gamma, res_eq, diffres="int"):
    """Extract FMU-oriented variable metadata for residual generator.

    Parameters
    ----------
    model : DiagnosisModel
        Symbolic diagnosis model.
    Gamma : Gamma
        MSO selection for the model.
    res_eq : int
        Index of the equation to be used as residual in the generated FMU.
    diffres : str, optional
        Residual treatment for differential residual equations, either "int" or
        "der". Defaults to "int".

    Returns
    -------
    Object
        Metadata for the given residual generator, including lists of parameters,
        signals, and state keys used by the residual equation.
    """
    if model.modelType != "Symbolic":
        raise ValueError("extract_fmu_variable_metadata requires a Symbolic model")

    # Local import avoids pulling CodeGeneration at module import time.
    from faultdiagnosistoolbox.CodeGeneration import DVar, IVar, UsedVars
    import faultdiagnosistoolbox as fdt

    res_gen_eq = _residual_equation_set(Gamma, res_eq)
    exact_i_state, exact_d_state = _matching_state_lists(
        model, Gamma, is_differential_constraint=fdt.IsDifferentialConstraint, ivar=IVar
    )
    residual_i_state, residual_d_state = _residual_state_lists(
        model,
        res_eq,
        diffres=diffres,
        is_differential_constraint=fdt.IsDifferentialConstraint,
        dvar=DVar,
        ivar=IVar,
    )
    state_keys = _ordered_unique(
        np.concatenate(
            (
                exact_i_state,
                residual_i_state,
                exact_d_state,
                residual_d_state,
            )
        )
    )

    return {
        "parameters": UsedVars(model.syme[res_gen_eq], model.parameters),
        "signals": UsedVars(model.syme[res_gen_eq], model.z),
        "state_keys": state_keys,
    }


def _residual_equation_set(gamma, res_eq):
    """Build the equation set used by a sequential residual generator."""
    res_gen_eqs = np.array([], dtype=np.int64)
    for hall_component in gamma.matching:
        res_gen_eqs = np.concatenate((res_gen_eqs, hall_component.row))
    return np.concatenate((res_gen_eqs, [res_eq]))


def _ordered_unique(values):
    """Return a deterministic list without duplicates while preserving order."""
    return list(dict.fromkeys(values.tolist()))


def _matching_state_lists(model, gamma, is_differential_constraint, ivar):
    """Collect state variables used by the exactly determined matching."""
    i_state = []
    d_state = []

    for hall_component in gamma.matching:
        if hall_component.matchType == "int":
            for expr, variable in zip(
                model.syme[hall_component.row], np.array(model.x)[hall_component.col]
            ):
                if is_differential_constraint(expr):
                    i_state.append(variable)
        elif hall_component.matchType == "der":
            d_state.append(ivar(model.syme[hall_component.row[0]]))
        elif hall_component.matchType == "mixed":
            for expr, variable in zip(
                model.syme[hall_component.row], np.array(model.x)[hall_component.col]
            ):
                if is_differential_constraint(expr):
                    if variable == ivar(expr):
                        i_state.append(variable)
                    else:
                        d_state.append(ivar(expr))

    return i_state, d_state


def _residual_state_lists(
    model, res_eq, diffres, is_differential_constraint, dvar, ivar
):
    """Collect state variables used by the residual equation."""
    expr = model.syme[res_eq]
    if not is_differential_constraint(expr):
        return [], []

    if diffres == "der":
        return [], [dvar(expr)]

    return [ivar(expr)], []
