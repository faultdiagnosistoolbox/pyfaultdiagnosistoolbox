import numpy as np
import os
import re
import platform

# Initial model exploration
from watertank_model import model

msos = model.MSO()
mtes = model.MTES()

I_mixed = model.IsolabilityAnalysis()
I_der = model.IsolabilityAnalysis(causality="der")
I_int = model.IsolabilityAnalysis(causality="int")

ts = [0, 1, 5, 10, 13]
red = [1, 1, 5, 1, 0]
FSM = model.FSM([msos[ti] for ti in ts])


def test_define_model():
    assert model

    dm = model.GetDMParts()
    assert dm.M0 == []
    assert dm.Mm.row == []
    assert dm.Mm.col == []
    assert np.all(dm.Mp.row == np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]))
    assert np.all(dm.Mp.col == np.array([0, 1, 2, 3, 4, 5, 6]))


def test_isolability_analysis():
    assert np.all(np.abs((I_mixed - np.eye(6))) < 1e-10)
    assert np.all(np.abs((I_der - np.eye(6))) < 1e-10)
    assert np.all(np.abs((I_int - np.eye(6))) < 1e-10)


def test_mso():
    oi_mso = np.array([model.IsObservable(mi) for mi in msos])
    li_mso = np.array([model.IsLowIndex(mi) for mi in msos])
    oi_mtes = np.array([model.IsObservable(mi) for mi in mtes])
    li_mtes = np.array([model.IsLowIndex(mi) for mi in mtes])

    assert len(msos) == 17
    assert len(mtes) == 10
    assert np.sum(oi_mso) == 17
    assert np.sum(oi_mtes) == 10
    assert np.sum(li_mtes) == 8
    assert np.sum(li_mso) == 15


def test_fsm():
    assert np.all(
        np.array(FSM)
        == np.array(
            [[1, 0, 0, 0, 1, 0], [0, 1, 1, 0, 0, 0], [0, 1, 0, 1, 1, 0], [0, 1, 0, 0, 0, 1], [1, 0, 1, 1, 0, 1]]
        )
    )

    assert np.all(I_mixed == model.IsolabilityAnalysisArrs([msos[ti] for ti in ts]))


def test_codegen():
    for k, (test, redIdx) in enumerate(zip(ts, red)):
        mso = msos[test]
        red_k = mso[redIdx]
        m0 = [e for e in mso if e != red_k]
        resName = f"r{k + 1}"
        Gamma = model.Matching(m0)
        print(f"Generating code for {resName} ...", end="")
        try:
            model.SeqResGen(Gamma, red_k, resName, language="C", batch=True)
            assert True
        except:
            assert False

    for k, _ in enumerate(zip(ts, red)):
        resName = f"r{k + 1}"
        print(f"Compiling residual generator: {resName} ... ", end="")
        if platform.system() == "Windows":
            compile_cmd = f"python {resName}_setup.py build_ext --inplace"
        else:
            compile_cmd = f"python {resName}_setup.py build_ext --inplace > /dev/null"
        if os.system(compile_cmd) == 0:
            assert True
        else:
            assert False

    # Cleanup
    files_to_remove = [fi for fi in os.listdir() if re.match("r[1-5]*", fi)] + [
        "src/" + fi for fi in os.listdir("src") if re.match("r[1-5]*", fi)
    ]
    for fi in files_to_remove:
        os.remove(fi)
