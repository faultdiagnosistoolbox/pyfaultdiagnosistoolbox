import faultdiagnosistoolbox as fdt
import numpy as np
import pytest
import sympy as sym

# Model (6) in Frisk, Krysander, and Åslund, "Analysis and design of diagnosis
# systems based on the structural differential index", IFAC World Congress, 2017.
model_def = {
    "type": "VarStruc",
    "x": ["dx11", "dx12", "x11", "x12", "x2"],
    "f": ["f1", "f2", "f3"],
    "z": ["y1", "y2", "y3"],
}
model_def["rels"] = [
    ["dx11", "x2", "f1"],
    ["dx12", "x2"],
    ["y1", "x11", "f2"],
    ["y2", "x12"],
    ["y3", "x2", "f3"],
    fdt.DiffConstraint("dx11", "x11"),
    fdt.DiffConstraint("dx12", "x12"),
]
model = fdt.DiagnosisModel(model_def, name="Index example")

M1 = [1, 3, 4, 6]
M2 = [0, 2, 4, 5]
M3 = [0, 1, 2, 3, 5, 6]


def test_exactly_determined_index():
    assert model.DifferentialIndex([1, 3, 4]) == 1
    assert model.DifferentialIndex([3, 4, 6]) == 2


def test_mso_index():
    assert model.DifferentialIndex(M1) == 1
    assert model.DifferentialIndex(M2) == 1
    assert model.DifferentialIndex(M3) == 2


def test_residual_equation_index():
    assert (model.ResidualEquationIndex(M1) == [2, 1, 2, 1]).all()
    assert (model.ResidualEquationIndex(M3) == 2).all()


def test_invalid_sets():
    with pytest.raises(ValueError):
        model.DifferentialIndex([0])  # underdetermined
    with pytest.raises(ValueError):
        model.MinimalIndexResidualSets([0])
    with pytest.raises(ValueError):
        model.ResidualEquationIndex([1, 3, 4])
    with pytest.raises(ValueError):
        model.DifferentialIndex(M1 + [M1[0]])  # repeated equation


def test_equation_set_arguments():
    # An empty set means all equations, as in the other methods
    assert model.DifferentialIndex([]) == model.DifferentialIndex()
    assert model.DifferentialIndex(set(M1)) == model.DifferentialIndex(M1)
    assert model.DifferentialIndex(range(7)) == model.DifferentialIndex()


def test_pantelides_requires_square():
    assert model.Pantelides() is None  # 7 equations, 5 unknowns


def test_residuals():
    idx, R = model.DifferentialIndex(M1, residuals=True)
    assert idx == 1 and len(R) == 1
    assert model.DifferentialIndex(np.setdiff1d(M1, R)) == 1
    idx, R = model.DifferentialIndex([1, 3, 4], residuals=True)
    assert idx == 1 and len(R) == 0


def test_pso():
    # The whole model is a PSO set with redundancy 2
    idx, R = model.DifferentialIndex(residuals=True)
    idx_all, Rs = model.MinimalIndexResidualSets()
    assert idx == idx_all == 1
    assert any(set(R) == set(Ri) for Ri in Rs)
    assert all(model.DifferentialIndex(np.setdiff1d(np.arange(7), Ri)) == 1 for Ri in Rs)


def test_minimal_index_residual_sets_mso():
    idx, Rs = model.MinimalIndexResidualSets(M1)
    assert idx == 1
    assert sorted(int(Ri[0]) for Ri in Rs) == [3, 6]  # residual equations c22 and c32


def test_enumeration_fallback():
    # Equation 3 contains two derivatives, so all residual sets are enumerated
    m = fdt.DiagnosisModel({"type": "MatrixStruc", "X": [[0, 0, 1], [0, 1, 0], [1, 0, 0], [1, 2, 2]], "F": [], "Z": []})
    assert m.DifferentialIndex() == 1
    assert m.DifferentialIndex(residuals=True)[1].tolist() == [3]


@pytest.mark.parametrize("m", [fdt.models.electric_motor, fdt.models.three_tank, fdt.models.induction_motor])
def test_agrees_with_low_index(m):
    for mso in m.MSO():
        assert (m.DifferentialIndex(mso) <= 1) == m.IsLowIndex(mso)


@pytest.mark.parametrize("m", [fdt.models.electric_motor, fdt.models.three_tank, fdt.models.induction_motor])
def test_agrees_with_enumeration(m):
    sets = [np.sort(mso) for mso in m.MSO()] + [np.sort(m.GetDMParts().Mp.row)]
    for eq in sets:
        idx, R = m.DifferentialIndex(eq, residuals=True)
        idx_all, Rs = m.MinimalIndexResidualSets(eq)
        assert idx == idx_all
        assert any(set(R) == set(Ri) for Ri in Rs)


def test_random_models():
    rng = np.random.default_rng(0)
    checked = 0
    while checked < 200:
        n1, n2 = int(rng.integers(1, 4)), int(rng.integers(0, 3))
        nx = 2 * n1 + n2
        A = (rng.random((n1 + n2 + int(rng.integers(1, 4)), nx)) < 0.35).astype(int)
        D = np.zeros((n1, nx), dtype=int)
        for k in range(n1):
            D[k, k], D[k, n1 + k] = 3, 2
        X = np.vstack([A, D])
        X = X[np.any(X > 0, axis=1)][:, np.any(X > 0, axis=0)]
        m = fdt.DiagnosisModel({"type": "MatrixStruc", "X": X.tolist(), "F": [], "Z": []})
        if X.shape[1] == 0 or m.srank() < X.shape[1]:
            continue
        idx, R = m.DifferentialIndex(residuals=True)
        idx_all, Rs = m.MinimalIndexResidualSets()
        assert idx == idx_all
        assert any(set(R) == set(Ri) for Ri in Rs)
        res = m.ResidualSetAnalysis()
        assert res[0]["deg"] == idx - 1
        assert [r["deg"] for r in res if set(r["R"]) == set(R)] == [idx - 1]
        checked += 1


# Two storages with a common inflow u, level sensors y1, y2 and a flow sensor y3
dx1, dx2, x1, x2, y1, y2, y3, u = sym.symbols("dx1 dx2 x1 x2 y1 y2 y3 u")
storages = fdt.DiagnosisModel(
    {
        "type": "Symbolic",
        "x": ["dx1", "dx2", "x1", "x2"],
        "f": [],
        "z": ["y1", "y2", "y3", "u"],
        "rels": [
            -y2 + x2,
            -y1 + x1,
            dx1 + dx2 - u,
            fdt.DiffConstraint("dx1", "x1"),
            fdt.DiffConstraint("dx2", "x2"),
            -y3 + dx1,
        ],
    }
)


def test_residual_set_analysis_mso():
    res = storages.ResidualSetAnalysis(range(5))
    assert len(res) == 5
    assert all(r["index"] == 2 and r["deg"] == 1 for r in res)
    assert res[0]["R"].tolist() == [0]
    assert res[0]["zorder"] == {"y1": 1, "y2": 0, "u": 0}


def test_residual_set_analysis_pso():
    res = storages.ResidualSetAnalysis(range(6))
    assert len(res) == 11
    assert res[0]["R"].tolist() == [0, 1] and res[0]["deg"] == 0
    r45 = [r for r in res if r["R"].tolist() == [3, 4]][0]
    assert r45["index"] == 1 and r45["deg"] == 1
    assert r45["zorder"] == {"y1": 1, "y2": 1, "y3": 0, "u": 0}


def test_residual_set_analysis_bounds():
    assert [r["R"].tolist() for r in storages.ResidualSetAnalysis(range(6), max_zorder=0)] == [[0, 1]]
    res = storages.ResidualSetAnalysis(range(6), max_zorder={"y2": 0})
    assert all(r["zorder"]["y2"] == 0 for r in res)
    assert [0, 3] in [r["R"].tolist() for r in res]
    with pytest.raises(ValueError):
        storages.ResidualSetAnalysis(range(6), max_zorder={"y4": 0})


def test_residual_set_analysis_square():
    res = storages.ResidualSetAnalysis([1, 2, 3, 4])
    assert len(res) == 1 and res[0]["R"].tolist() == []
    assert res[0]["index"] == 2 and res[0]["deg"] == 1


def test_residual_set_analysis_mapping_bound():
    from types import MappingProxyType

    expected = [r["R"].tolist() for r in storages.ResidualSetAnalysis(range(6), max_zorder={"y2": 0})]
    res = storages.ResidualSetAnalysis(range(6), max_zorder=MappingProxyType({"y2": 0}))
    assert [r["R"].tolist() for r in res] == expected
