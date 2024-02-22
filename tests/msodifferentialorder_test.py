import faultdiagnosistoolbox as fdt
import sympy as sym
import numpy as np

x1, x2, x3, dx1, dx2, u, y, f = sym.symbols(["x1", "x2", "x3", "dx1", "dx2", "u", "y", "f"])
model_def = {"type": "Symbolic", "x": ["x1", "x2", "x3", "dx1", "dx2"], "f": ["f"], "z": ["y", "u"]}
model_def["rels"] = [
    -dx1 + -x1 + x2,
    -dx2 + -x2 + x3,
    -x3 + u,
    -y + x1 + f,
    fdt.DiffConstraint("dx1", "x1"),
    fdt.DiffConstraint("dx2", "x2"),
]

model_small = fdt.DiagnosisModel(model_def, name="Small model")
model_induction = fdt.models.induction_motor


def test_small_diff_order():
    mso = model_small.MSO()[0]
    eqorder, zidx, zorder, fidx, forder = model_small.MSOdifferentialOrder(mso)
    assert eqorder == [0, 1, 2, 0, 0, 1]
    assert (zidx == [0, 1]).all()
    assert (zorder == [2, 0]).all()
    assert (fidx == [0]).all()
    assert (forder == [0]).all()


def test_induction_diff_order():
    msos = model_induction.MSO()
    mso = np.sort(msos[0])
    eqorder, zidx, zorder, fidx, forder = model_induction.MSOdifferentialOrder(mso)
    assert eqorder == [2, 2, 1, 1, 1, 2, 1, 1, 1, 3, 2, 2]
    assert (zidx == [0, 1, 2, 3, 4]).all()
    assert (zorder == [2, 1, 3, 2, 2]).all()
    assert (fidx == [0, 1]).all()
    assert (forder == [0, 0]).all()

    mso = np.sort(msos[8])
    eqorder, zidx, zorder, fidx, forder = model_induction.MSOdifferentialOrder(mso)
    assert eqorder == [2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 3, 3]
    assert (zidx == [0, 1, 2, 3]).all()
    assert (zorder == [2, 2, 3, 3]).all()
    assert (fidx == [0, 1]).all()
    assert (forder == [0, 0]).all()
