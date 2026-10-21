import sympy as sym
import faultdiagnosistoolbox as fdt

modelDef = {}
modelDef["type"] = "Symbolic"  # should always be 'Symbolic' for symbolic models
modelDef["x"] = ["x1", "x2", "x3"]  # Unknown variables
modelDef["f"] = ["f1", "f2"]  # Fault variables
modelDef["z"] = ["y1", "y2"]  # Measured variables (known)


# Make symbolic objects of all variables/parameters before writing down equations.
x1, x2, x3 = sym.symbols(modelDef["x"])
f1, f2 = sym.symbols(modelDef["f"])
y1, y2 = sym.symbols(modelDef["z"])

modelDef["rels"] = [
    -x3 + x1 + x2,
    -y1 + x3 + f1,
    -y2 + x3 + f2,
]

model = fdt.DiagnosisModel(modelDef, name="Three tank system")
