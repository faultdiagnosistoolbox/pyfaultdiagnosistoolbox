"""This file contains some definitions for testing models for use by the
model_transforming_test.py file"""

import faultdiagnosistoolbox as fdt
import sympy as sym
import faultdiagnosistoolbox.ModelTransforming as mt

# Bad format test model
bad_subexpr_modelDef = {
    "type": "Symbolic",
    "x": ["x1", "x2"],
    "f": ["f1", "f2"],
    "z": ["y1", "y2"],
}
sym.var(bad_subexpr_modelDef["x"])
sym.var(bad_subexpr_modelDef["f"])
sym.var(bad_subexpr_modelDef["z"])

bad_subexpr_modelDef["rels"] = [
    -y1 + x1 + x2 + f1,
    -y2 + x1 + x2 + f2,
]

bad_subexpr_model = fdt.DiagnosisModel(bad_subexpr_modelDef, name="Non optimized model")

# Good format test model
good_subexpr_modelDef = {
    "type": "Symbolic",
    "x": ["x1", "x2", "x3"],
    "f": ["f1", "f2"],
    "z": ["y1", "y2"],
}
sym.var(good_subexpr_modelDef["x"])
sym.var(good_subexpr_modelDef["f"])
sym.var(good_subexpr_modelDef["z"])

good_subexpr_modelDef["rels"] = [
    -x3 + x1 + x2,
    -y1 + x3 + f1,
    -y2 + x3 + f2,
]

good_subexpr_model = fdt.DiagnosisModel(good_subexpr_modelDef, name="Optimized model")

# Model for identification test
subexpr_id_modelDef = {
    "type": "Symbolic",
    "x": ["x1", "x2", "x3"],
    "f": ["f1", "f2", "f3"],
    "z": ["y1", "y2", "y3"],
}
sym.var(subexpr_id_modelDef["x"])
sym.var(subexpr_id_modelDef["f"])
sym.var(subexpr_id_modelDef["z"])

subexpr_id_modelDef["rels"] = [
    -y1 + x1 * x2 + x3 + f1,
    -y2 + x1 * x2 + x3 + f2,
    -y3 + x1 * x2 + f3,
]

# Correct reduction and adjustments for identification test model
sym.var(["x0", "x4"])
subexpr_id_correct_reduction = [(x0, x1 * x2), (x4, x0 + x3)]

subexpr_id_reduced_modelDef = {
    "type": "Symbolic",
    "x": ["x0", "x1", "x2", "x3", "x4"],
    "f": ["f1", "f2", "f3"],
    "z": ["y1", "y2", "y3"],
}
subexpr_id_reduced_modelDef["rels"] = [
    -x0 + x1 * x2,
    -x4 + x0 + x3,
    -y1 + x4 + f1,
    -y2 + x4 + f2,
    -y3 + x0 + f3,
]

# Model with 12 equations for performance test
twelve_eq_modelDef = {
    "type": "Symbolic",
    "x": ["h1", "h2", "xc1", "xl2", "xl3", "dh1", "dh2", "x2"],
    "f": ["fh2", "ff1", "fc1", "fl2", "fl3", "fa"],
    "z": ["y1", "y2", "y3", "y4", "u"],
    "parameters": ["d1", "d2", "d3", "d4", "d5", "d6"],
}
sym.var(twelve_eq_modelDef["x"])
sym.var(twelve_eq_modelDef["f"])
sym.var(twelve_eq_modelDef["z"])
sym.var(twelve_eq_modelDef["parameters"])

twelve_eq_modelDef["rels"] = [
    -dh1 + d1 * u - d2 * xc1 * sym.sqrt(h1) + fa,
    -dh2 + d3 * xc1 * xl2 * sym.sqrt(h1) - d4 * sym.sqrt(h2),
    -y1 + h1,
    -y2 + h2 + fh2,
    -y3 + d5 * xc1 * sym.sqrt(h1) + ff1,
    -y4 + d6 * xl3 * sym.sqrt(h2),
    -xc1 + 1 - fc1,
    -xl2 + 1 - fl2,
    -xl3 + 1 - fl3,
    -x2 + d1 * u - h1,
    fdt.DiffConstraint("dh1", "h1"),
    fdt.DiffConstraint("dh2", "h2"),
]

twelve_eq_model = fdt.DiagnosisModel(
    twelve_eq_modelDef,
    name="Water Tank Model, reduced case, modified for modeltransforming test",
)
