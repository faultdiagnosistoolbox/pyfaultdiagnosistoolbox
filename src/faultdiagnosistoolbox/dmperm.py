"""Dulmage-Mendelsohn an MSO functionality."""

import faultdiagnosistoolbox.dmpermlib as dmpermlib
import scipy.sparse as sp
import numpy as np
import copy
from dataclasses import dataclass, field


def CSCDict(A):
    """Compressed matrix format."""
    return {
        "nzmax": A.nnz,
        "m": A.shape[0],
        "n": A.shape[1],
        "p": A.indptr.astype(np.int64),
        "i": A.indices.astype(np.int64),
        "x": A.data.astype(np.float64),
        "nz": -1,
    }


def dmperm(A) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Dulmage-Mendelsohn decomposition of matrix.

    p, q, r, s, cc, rr, m = dmperm(A) finds the Dulmage-Mendelsohn decomposition
    of A.  p and q are row and column permutation vectors, respectively
    such that A[p][:, q] has block upper triangular form.  r and s are vectors
    indicating the block boundaries for the fine decomposition.  cc and rr
    are vectors of length five indicating the block boundaries of the
    coarse decomposition. m is a maximum matching such that m[j] = i if column j
    is matched to row i, or -1 if column j is unmatched.
    """

    if sp.issparse(A):
        return dmpermlib.dmperm(A)
    else:
        return dmpermlib.dmperm(sp.csc_matrix(A))


def srank(A) -> int:
    """Compute structural rank of matrix."""
    # _, _, _, _, rr, _ = dmperm(A)
    #    return rr[3]

    # Extract matching vector m(i) = j
    # if column j is matched to row i, or -1 if column j is unmatched.
    #
    # m = np.full(hod.shape[1], -1)
    # m[p[rr[0]:rr[1]]] = q[cc[1]:cc[2]]
    # m[p[rr[1]:rr[2]]] = q[cc[2]:cc[3]]
    # m[p[rr[2]:rr[3]]] = q[cc[3]:cc[4]]
    # return sum(m >= 0)

    # _, _, _, _, cc, _, _ = dmperm(A)
    # return sum(np.diff(cc[1:]))
    _, _, _, _, _, _, m = dmperm(A)
    return sum(m >= 0)


@dataclass
class EqBlock:
    """EqBlock class."""

    row: np.array = field(default_factory=lambda: np.array([], dtype=np.int64))
    col: np.array = field(default_factory=lambda: np.array([], dtype=np.int64))


@dataclass
class DMResult:
    """Dulmage-Mendelsohn decomposition base class."""

    Mm: EqBlock = field(default_factory=lambda: EqBlock([], []))
    M0: list = field(default_factory=lambda: [])
    Mp: EqBlock = field(default_factory=lambda: EqBlock([], []))
    rowp: list = field(default_factory=lambda: [])
    colp: list = field(default_factory=lambda: [])
    M0eqs: list = field(default_factory=lambda: [])
    M0vars: list = field(default_factory=lambda: [])


def MSO(X):
    """Compute set of MSO sets."""
    if sp.issparse(X):
        return dmpermlib.MSO(X)
    else:
        return dmpermlib.MSO(sp.csc_matrix(X))


def GetDMParts(X) -> DMResult:
    """Compute Dulmage-Mendelsohn decomposition."""
    if sp.issparse(X):
        p, q, r, s, _, _, _ = dmperm(X)
    else:
        p, q, r, s, _, _, _ = dmperm(sp.csc_matrix(X))

    m = X.shape[0]
    n = X.shape[1]

    res = DMResult()
    if (p.size == 0) or (q.size == 0):
        if m > n:
            res.Mm = EqBlock([], [])
            res.Mp = EqBlock(np.arange(0, m), np.arange(0, n))
        else:
            res.Mp = EqBlock([], [])
            res.Mm = EqBlock(np.arange(0, m), np.arange(0, n))

        res.rowp = np.arange(1, m)
        res.colp = np.arange(1, n)
        res.M0 = []
    else:
        idx = 0
        if s[1] > r[1]:  # Existance of M-
            res.Mm = EqBlock(np.sort(p[r[0] : r[1]]), np.sort(q[s[0] : s[1]]))
            idx = idx + 1
        else:
            res.Mm = EqBlock([], [])

        while (idx < r.size - 1) and (r[idx + 1] - r[idx]) == (s[idx + 1] - s[idx]):
            # M0 block exists
            res.M0.append(EqBlock(np.sort(p[r[idx] : r[idx + 1]]), np.sort(q[s[idx] : s[idx + 1]])))
            idx = idx + 1

        if idx < r.size - 1:  # M+ exists
            res.Mp = EqBlock(np.sort(p[r[idx] : r[idx + 1]]), np.sort(q[s[idx] : s[idx + 1]]))
        else:
            res.Mp = EqBlock([], [])
        res.rowp = p
        res.colp = q

        res.M0eqs = np.array([], dtype=np.int64)
        res.M0vars = np.array([], dtype=np.int64)

        for hc in res.M0:
            res.M0eqs = np.concatenate((res.M0eqs, hc.row))
            res.M0vars = np.concatenate((res.M0vars, hc.col))
    return res


def PSODecomposition(X) -> dict[list[EqBlock], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute PSO decomposition."""
    if not IsPSO(X):
        print("PSO Decomposition for PSO structures only, exiting...")
        return

    n = X.shape[0]
    m = X.shape[1]

    delrows = np.arange(0, n)
    eqclass = []
    trivclass = np.array([], dtype=np.int64)
    Mi = np.array([])

    while len(delrows) > 0:
        temprow = np.array([x for x in np.arange(0, n) if not x == delrows[0]])
        dm = GetDMParts(X[temprow, :])

        if len(dm.M0vars) > 0:
            ec_row = np.sort(np.append(temprow[dm.M0eqs], delrows[0]))
            ec = EqBlock(ec_row, np.sort(dm.M0vars))
            eqclass.append(ec)
            Mi = np.concatenate((Mi, ec.col)).astype(np.int64)
            delrows = [x for x in delrows if x not in ec.row]
        else:
            trivclass = np.concatenate((trivclass, [delrows[0]]))
            delrows = delrows[1 : len(delrows)]
    X0 = np.sort([x for x in np.arange(0, m) if not (x in Mi)])
    if len(X0) == 0:
        X0 = []

    res = {"eqclass": eqclass, "trivclass": trivclass, "X0": X0}

    p = np.array([], dtype=np.int64)
    q = np.array([], dtype=np.int64)

    for ec in eqclass:
        p = np.concatenate((p, ec.row))
        q = np.concatenate((q, ec.col))
    p = np.concatenate((p, trivclass))
    q = np.concatenate((q, X0))

    res["p"] = p
    res["q"] = q
    return res


def IsPSO(X, eq=None) -> bool:
    """Return True if PSO."""
    if eq is None:
        eqs = np.arange(0, X.shape[0])
    else:
        eqs = eq

    dm = GetDMParts(X[eqs, :])
    return (len(dm.Mm.row) == 0) and (len(dm.M0) == 0)


def IsObservable(Xin, eq=None) -> bool:
    """Return True if observable."""

    def DiffConstraints(X):
        diffEqs = np.where(np.any(X == 2, axis=1))[0]
        algEqs = np.array([ei for ei in np.arange(0, ne) if not (ei in diffEqs)])
        x1Idx = np.argwhere(np.any(X[diffEqs, :] == 2, axis=0)).flatten()
        dx1Idx = np.array(list(map(lambda x1i: np.argwhere(X[np.argwhere(X[:, x1i] == 2)[0][0], :] == 3)[0][0], x1Idx)))
        x2Idx = np.array([xi for xi in np.arange(0, nx) if (not (xi in x1Idx)) and (not (xi in dx1Idx))])
        return diffEqs, algEqs, x1Idx, dx1Idx, x2Idx

    if eq is None:
        X = Xin
    else:
        X = Xin[eq, :]
        # Remove zero columns
        zc = np.where(np.all(X == 0, axis=0))[0]
        X = np.delete(X, zc, axis=1)

    ne = X.shape[0]
    nx = X.shape[1]

    diffEqs, algEqs, x1Idx, dx1Idx, x2Idx = DiffConstraints(X)
    n1 = len(x1Idx)
    n2 = len(x2Idx)
    #    nalg = len(algEqs)

    # Model format:
    #  x1' = dx1
    #  0   = A11 x1 + A12 x2 + A13 dx1
    #
    #  AC-sF=0

    A = X[algEqs, :][:, np.hstack((x1Idx, x2Idx, dx1Idx)).astype(np.int64)]
    A11 = A[:, np.arange(0, n1)]
    A12 = A[:, np.arange(n1, n1 + n2)]
    A13 = A[:, np.arange(n1 + n2, 2 * n1 + n2)]
    AC = np.vstack((np.hstack((np.zeros((n1, n1 + n2)), np.eye(n1))), A)).astype(np.int64)
    F = np.vstack((np.hstack((np.eye(n1), np.zeros((n1, n1 + n2)))), np.zeros((ne - n1, n1 + n2 + n1)))).astype(
        np.int64
    )

    # Condition 1: srank(A-sF)=n1+n2
    obs = srank(np.hstack((np.logical_or(A11, A13), A12))) == n1 + n2

    if obs:
        # Condition 2: v(AC)=2*n1+n2
        obs = srank(AC) == 2 * n1 + n2
    if obs:
        # Condition 3: no s-argc in Hall components of AC-sF
        G3 = np.maximum(AC, 2 * F)
        dm = GetDMParts(G3)
        if len(dm.M0) > 0:
            hc_sarcs = list(map(lambda hc: np.any(G3[hc.row, :][:, hc.col] == 2), dm.M0))
            obs = len(hc_sarcs) == 0 or (not np.any(hc_sarcs))

    return obs


def IsHighIndex(X, eq=None) -> bool:
    """Return True if high differential index."""
    if eq is None:
        eq = np.arange(0, X.shape[0])
    X1 = X[eq, :]

    col_d1 = np.any(X1 == 3, axis=0)
    row_a = np.all(X1 <= 1, axis=1)
    row_d = np.any(X1 > 1, axis=1)
    col_2 = np.all(X1[row_d, :] == 0, axis=0)
    Xhod = X1[row_a, :][:, np.logical_or(col_d1, col_2)]
    nz = np.sum(np.all(Xhod == 0, axis=0))
    return srank(Xhod) < Xhod.shape[1] - nz


def IsLowIndex(X, eq=None) -> bool:
    """Return True if low differential index."""
    if eq is None:
        eq = np.arange(0, X.shape[0])
    return not IsHighIndex(X, eq)


def Mplus(X, causality="mixed"):
    """Compute over-determined part."""

    def Gp(Gi):
        dm = GetDMParts(Gi[2])
        if len(dm.Mp.row) == 0 or len(dm.Mp.col) == 0:
            return np.array([]), np.array([]), np.array([[]])
        G1 = copy.deepcopy(Gi)
        return (G1[0][dm.Mp.row], G1[1][dm.Mp.col], G1[2][dm.Mp.row, :][:, dm.Mp.col])

    def Gm(Gi):
        dm = GetDMParts(Gi[2])
        if len(dm.Mm.row) == 0 or len(dm.Mm.col) == 0:
            return np.array([]), np.array([]), np.array([[]])
        G1 = copy.deepcopy(Gi)
        return (G1[0][dm.Mm.row], G1[1][dm.Mm.col], G1[2][dm.Mm.row, :][:, dm.Mm.col])

    def G0(Gi):
        dm = GetDMParts(Gi[2])
        if len(dm.M0eqs) == 0 or len(dm.M0vars) == 0:
            return np.array([]), np.array([]), np.array([[]])
        G1 = copy.deepcopy(Gi)
        return (G1[0][dm.M0eqs], G1[1][dm.M0vars], G1[2][dm.M0eqs, :][:, dm.M0vars])

    # def Gadd(G1, G2):
    #     c1, x1, A1 = copy.deepcopy(G1)
    #     c2, x2, A2 = copy.deepcopy(G2)
    #
    #     c1 = np.unique(np.concatenate((c1, c2)))
    #     x1 = np.unique(np.concatenate((x1, x2)))
    #     A1[A2 == 1] = 1
    #     A1[A2 == 2] = 2
    #     return c1, x1, A1

    def CGX(Gi, X):
        if len(X) == 0:
            return []
        xIdx = [np.where(Gi[1] == ii)[0][0] for ii in X]
        idx = np.unique([e[0] for e in np.argwhere(Gi[2][:, xIdx] > 0)])
        if len(idx) == 0:
            return []

        return Gi[0][idx]

    def CGE(Gi, Ei):
        return Gi[0][np.where(np.any(Ei, axis=1))[0]]

    def GsubC(Gi, C):
        c, x, A = copy.deepcopy(Gi)
        Cidx = list(map(lambda ci: np.argwhere(Gi[0] == ci)[0][0], C))
        A = np.delete(A, Cidx, axis=0)
        c = np.array([ci for ci in c if not (ci in C)])
        return c, x, A

    def GsubX(Gi, X):
        c, x, A = copy.deepcopy(Gi)
        Xidx = list(map(lambda xi: np.argwhere(Gi[1] == xi)[0][0], X))
        A = np.delete(A, Xidx, axis=1)
        x = np.array([xi for xi in x if not (xi in X)])
        return c, x, A

    Xc = np.array(X.copy())
    # Represent graph as a tuple G=(constraints,variables, adjacency matrix)
    G = (np.arange(0, Xc.shape[0], dtype=np.int64), np.arange(0, Xc.shape[1], dtype=np.int64), Xc)

    if causality == "mixed":
        return Gp(G)[0]

    elif causality == "int":
        while True:
            # G := G+
            G = Gp(G)

            # G1 = G - Ed
            G1 = copy.deepcopy(G)
            G1[2][G1[2] == 3] = 0  # Zero the derivative edges

            # G := G - C(G,X(G1-))
            dm = GetDMParts(G1[2])

            if len(dm.Mm.col) > 0 and len(dm.Mm.row) > 0:
                G1m = Gm(G1)
                G = GsubC(G, CGX(G, G1m[1]))
            else:
                break
        return G[0]
    elif causality == "der":
        Xc = []
        while True:
            # Gnc := G - Xc
            Gnc = GsubX(G, Xc)

            # Gni := Gnc - C(Gnc,Ei)
            Ei = Gnc[2].copy()
            Ei[Ei != 2] = 0
            Gni = GsubC(Gnc, CGE(Gnc, Ei))

            # Xc:= Xc u X(Gni+ u Gni0)
            Gnip = Gp(Gni)
            Gni0 = G0(Gni)
            X1 = np.unique(np.concatenate((Gnip[1], Gni0[1]))).astype(np.int64)  # X1 = X(Gni+ u Gni0)
            Xc = np.unique(np.concatenate((Xc, X1))).astype(np.int64)

            if len(X1) == 0:
                break

        # G := (G-C(G,X\Xc))+
        X1 = np.array([xi for xi in G[1] if not (xi in Xc)])  # X1 = X\Xc
        Gres = Gp(GsubC(G, CGX(G, X1)))
        return Gres[0]
    else:
        return np.array([])


def MixedCausalityMatching(Gin):
    """Compute mixed causality matching."""

    def FindAdmissibleIntEdge(G):
        X = G[2]
        dm = GetDMParts(X)
        for hallComponent in dm.M0:
            Ei = np.argwhere(X[hallComponent.row, :][:, hallComponent.col] == 2)
            if len(Ei) > 0:
                Ei = Ei[0]  # Get the first one
                return np.array([G[0][hallComponent.row[Ei[0]]], G[1][hallComponent.col[Ei[1]]]])
        return []

    G = [Gin[0].copy(), Gin[1].copy(), Gin[2].copy()]
    Gamma = []
    e = FindAdmissibleIntEdge(G)
    while len(e) > 0:
        if len(Gamma) == 0:
            Gamma = [e]
        else:
            Gamma = np.concatenate((Gamma, [e]), axis=0)

        # G := G - C({e}) - X({e})
        rowIdx = np.argwhere(G[0] == e[0])[0, 0]
        colIdx = np.argwhere(G[1] == e[1])[0, 0]
        G[0] = np.delete(G[0], rowIdx)
        G[1] = np.delete(G[1], colIdx)
        G[2] = np.delete(G[2], rowIdx, axis=0)
        G[2] = np.delete(G[2], colIdx, axis=1)
        e = FindAdmissibleIntEdge(G)

    dm = GetDMParts(G[2])
    Gamma_p = np.stack((G[0][dm.rowp][:], G[1][dm.colp][:]), axis=1)
    return np.flipud(np.concatenate((Gamma, Gamma_p), axis=0))


def MTES(self):
    """Compute set of MTES sets."""
    S = {"eq": [], "f": [], "sr": []}
    m = MTES_initModel(self)  # overdetermined or empty
    if not (m is None) and m["sr"] > 0 and len(m["f"]) > 0:
        S = MTES_FindMTES(m, 0)
    #    return np.array(S['eq'], dtype=object)
    return S["eq"]


def MTES_storeFS(m):
    """internal."""
    #    eq = np.sort(np.hstack(m['e'])).tolist()
    #    f = np.sort(np.hstack(m['f'])).tolist()
    eq = np.sort(np.hstack(m["e"]))
    f = np.sort(np.hstack(m["f"]))
    return {"eq": [eq], "f": [f], "sr": [m["sr"]]}


def MTES_initModel(model):
    """internal."""
    dm = GetDMParts(model.X)
    row_over = dm.Mp.row
    col_over = dm.Mp.col

    ef = np.argwhere(np.any(model.F, axis=1)).flatten()
    ef = np.intersect1d(ef, row_over)
    idx = np.hstack((ef, np.setdiff1d(row_over, ef)))
    # m = {}
    # m['sr'] = len(row_over) - len(col_over)
    # m['X'] = model.X[idx, :][:, col_over]
    # m['f'] = list(map(lambda fi: np.argwhere(fi)[0], model.F[ef, :]))
    # m['delrow'] = 0
    # m['e'] = list(map(lambda ei: np.array([ei]), idx))
    if len(row_over) > 0:
        m = {
            "sr": len(row_over) - len(col_over),
            "X": model.X[idx, :][:, col_over],
            "f": list(map(lambda fi: np.argwhere(fi)[0], model.F[ef, :])),
            "delrow": 0,
            "e": list(map(lambda ei: np.array([ei]), idx)),
        }
    else:
        m = None

    return m


def MTES_GetPartialModel(m, rows):
    """internal."""
    n = {}
    variables = np.any(m["X"][rows, :], axis=0)
    n["X"] = m["X"][rows, :][:, variables]
    #    n['e'] = list(np.array(m['e'])[rows])
    n["e"] = [m["e"][k] for k in range(len(m["e"])) if k in rows]
    #    n['f'] = list(np.array(m['f'])[[ei for ei in rows if ei < len(m['f'])]])
    n["f"] = [m["f"][ei] for ei in rows if ei < len(m["f"])]
    n["sr"] = n["X"].shape[0] - n["X"].shape[1]
    n["delrow"] = np.sum(rows < m["delrow"])
    return n


def MTES_FindMTES(m, p):
    """internal."""
    m, _ = MTES_LumpExt(m, 0)
    if len(m["f"]) == 1:  # m is an MTES
        S = MTES_storeFS(m)
    else:  # recurse
        if p == 1:
            S = MTES_storeFS(m)
        else:
            S = {"eq": [], "f": [], "sr": []}
        row = m["delrow"]
        while len(m["f"]) > row:
            m, row = MTES_LumpExt(m, row)
        for delrow in np.arange(m["delrow"], len(m["f"])):
            # create the model where delrow has been removed
            m["delrow"] = delrow
            rows = np.delete(np.arange(0, m["X"].shape[0]), delrow)
            n = MTES_GetPartialModel(m, rows)
            Sn = MTES_FindMTES(n, p)  # make recursive call
            S = MTES_addResults(S, Sn)
    return S


def MTES_LumpExt(m, row):
    """MTES lumping."""
    n = {}

    no_rows = m["X"].shape[0]
    remRows = np.hstack((np.arange(0, row), np.arange(row + 1, no_rows)))
    remRowsf = np.hstack((np.arange(0, row), np.arange(row + 1, len(m["f"]))))

    dm = GetDMParts(m["X"][remRows, :])
    row_just = dm.M0eqs
    row_over = np.array(dm.Mp.row, dtype=np.int64)
    col_over = np.array(dm.Mp.col, dtype=np.int64)
    if len(row_just) > 0:
        eqcls = np.hstack((remRows[row_just], [row]))
        no_rows_before_row = np.sum(eqcls < row)
        row = row - no_rows_before_row
        no_rows_before = np.sum(eqcls < m["delrow"])
        n["delrow"] = m["delrow"] - no_rows_before
        eqclsf = eqcls[eqcls < len(m["f"])]
        row_overf = row_over[row_over < len(remRowsf)]
        if no_rows_before > 0:
            rowinsert = n["delrow"]
        else:
            rowinsert = row
        x1 = m["X"][remRows[row_over[0:rowinsert]], :][:, col_over]
        x3 = m["X"][remRows[row_over[rowinsert:]], :][:, col_over]
        x2 = np.any(m["X"][eqcls, :][:, col_over], axis=0).astype(np.int64)
        n["X"] = np.vstack((x1, x2, x3))

        foo = (
            remRows[row_over[0:rowinsert]].astype(np.int64)
            if len(row_over[0:rowinsert]) > 0
            else np.array([], dtype=np.int64)
        )
        #        e1 = list(np.array(m['e'], dtype=np.ndarray)[foo])
        e1 = [m["e"][idx] for idx in foo]
        #        e2 = list([np.hstack(np.array(m['e'], dtype=np.ndarray)[eqcls])])
        e2 = list([np.hstack([m["e"][idx] for idx in eqcls])])
        foo = (
            remRows[row_over[rowinsert:]].astype(np.int64)
            if len(row_over[rowinsert:]) > 0
            else np.array([], dtype=np.int64)
        )
        #        e3 = list(np.array(m['e'], dtype=np.ndarray)[foo])
        e3 = [m["e"][idx] for idx in foo]
        n["e"] = e1 + e2 + e3

        foo = (
            remRowsf[row_overf[0:rowinsert]].astype(np.int64)
            if len(row_overf[0:rowinsert]) > 0
            else np.array([], dtype=np.int64)
        )
        ef1 = list(np.array(m["f"], dtype=np.ndarray)[foo])
        ef2 = list([np.hstack(np.array(m["f"], dtype=np.ndarray)[eqclsf])])
        foo = (
            remRowsf[row_overf[rowinsert:]].astype(np.int64)
            if len(row_overf[rowinsert:]) > 0
            else np.array([], dtype=np.int64)
        )
        ef3 = list(np.array(m["f"], dtype=np.ndarray)[foo])
        n["f"] = ef1 + ef2 + ef3

        n["sr"] = m["sr"]

        if no_rows_before > 0:
            n["delrow"] = n["delrow"] + 1
    else:
        n = m

    row = row + 1
    return n, row


def MTES_addResults(S, Sn):
    """internal."""
    return {"eq": S["eq"] + Sn["eq"], "f": S["f"] + Sn["f"], "sr": S["sr"] + Sn["sr"]}
