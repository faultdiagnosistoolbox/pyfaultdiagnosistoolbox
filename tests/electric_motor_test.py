import numpy as np
import os

from electricmotor_model import model


msos = model.MSO()
mtes = model.MTES()


def test_msos():
    assert model
    assert len(msos) == 3
    assert len(mtes) == 3
    assert all([model.IsLowIndex(mi) for mi in msos])
    assert [len(mi) for mi in msos] == [4, 4, 6]


def test_psos():
    assert [model.IsPSO(mi) for mi in msos]
    assert [model.IsPSO(mi) for mi in mtes]


def test_integral():
    mso = np.sort(msos[2])
    red = mso[2]
    M0 = [ei for ei in mso if ei != red]
    Gamma_int = model.Matching(M0)
    try:
        model.SeqResGen(Gamma_int, red, "r_int")
        assert True
    except:
        assert False

    # Cleanup
    if "r_int.py" in os.listdir():
        os.remove("r_int.py")


def test_derivative():
    mso = np.sort(msos[2])
    red = mso[1]
    M0 = [ei for ei in mso if ei != red]
    Gamma_der = model.Matching(M0)
    try:
        model.SeqResGen(Gamma_der, red, "r_der")
        assert True
    except:
        assert False

    # Cleanup
    if "r_der.py" in os.listdir():
        os.remove("r_der.py")


def test_mixed():
    mso = np.sort(msos[2])
    red = mso[0]
    M0 = [ei for ei in mso if ei != red]
    Gamma_mixed = model.Matching(M0)

    try:
        model.SeqResGen(Gamma_mixed, red, "r_mixed")
        assert True
    except:
        assert False

    # Cleanup
    if "r_mixed.py" in os.listdir():
        os.remove("r_mixed.py")
