import numpy as np
import faultdiagnosistoolbox as fdt
import sympy as sym

# %% Illustrative example from Section 3 in
# Mattias Krysander, Daniel Jung, "Structural Methods for Testable
# Signal Sets in Data-Driven Fault Diagnosis",
# IFAC World Congress, 2026

tss_illustrative_model = {
    "type": "Symbolic",
    "x": ["I_0", "I_1", "I_2", "V", "V_I"],
    "f": ["f_v", "f_R0", "f_R1", "f_R2"],
    "z": ["V_0", "y_1", "y_2"],
    "parameters": ["R0", "R1", "R2"],
}

I0, I1, I2, V, VI = sym.var(tss_illustrative_model["x"])
f_v, f_R0, f_R1, f_R2 = sym.var(tss_illustrative_model["f"])
V0, y1, y2 = sym.var(tss_illustrative_model["z"])
R0, R1, R2 = sym.var(tss_illustrative_model["parameters"])

tss_illustrative_model["rels"] = [
    V - V0 - f_v,
    V - VI - (R0 + f_R0) * I0,
    I0 - I1 - I2,
    VI - I1 * (R1 + f_R1),
    VI - I2 * (R2 + f_R2),
    y1 - I1,
    y2 - I2,
]

model = fdt.DiagnosisModel(tss_illustrative_model, name="TSS_Illustrative")


def test_tss_illustrative_example():
    tss = [e.tolist() for e in model.TSS()]
    sol = [[0, 1, 2, 3, 4, 5, 6], [3, 4, 5, 6], [0, 1, 2, 3, 4, 6], [0, 1, 2, 3, 4, 5]]
    assert len(tss) == len(sol)
    assert all([s_i in sol for s_i in tss])
