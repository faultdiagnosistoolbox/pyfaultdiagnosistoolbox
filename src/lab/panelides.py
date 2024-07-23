"""
%load_ext autoreload
%autoreload 2
"""
import sympy as sym
import numpy as np
import faultdiagnosistoolbox as fdt
from faultdiagnosistoolbox.dmperm import GetDMParts


def GraphPantelides(Xi, V):
    def HighestOrderDerivatives(var_defs):
        hd = {}
        for dvi in var_defs:
            vi, der = dvi
            if vi in hd.keys():
                if hd[vi] < der:
                    hd[vi] = der
            else:
                hd[vi] = der
        return np.array([hd[v[0]] == v[1] for v in var_defs])

    def DifferentiateGraph(X, var_defs, der_eq):
        for der_i in der_eq:
            v_new = [(var_defs[k][0], var_defs[k][1] + 1) for k in np.argwhere(X[der_i] > 0).reshape(-1)]
            new_cols = 0
            for v_i in v_new:
                if not v_i in V:
                    var_defs.append(v_i)
                    new_cols += 1

            X = np.hstack((X, np.zeros((X.shape[0], new_cols))))
            X[der_i, :] = np.zeros(X.shape[1])
            X[der_i, [var_defs.index(vi) for vi in v_new]] = 1
        return V, X

    var_defs = V.copy()
    X = Xi.copy()
    nu = np.zeros(X.shape[0], dtype=np.int)
    while True:
        hod = HighestOrderDerivatives(var_defs)
        dm = GetDMParts(X[:, hod])
        if len(dm.Mp.row) > 0:
            der_eq = dm.Mp.row
            nu[der_eq] = nu[der_eq] + 1
            V, X = DifferentiateGraph(X, var_defs, der_eq)
        else:
            break
    s_idx = 1 * any([var_defs[k][1] == 0 for k in np.argwhere(hod).reshape(-1)]) + np.max(nu)
    return s_idx, nu

def Pantelides(model):
    if model.ne() == model.nx():
        der_col = np.argwhere(np.any(model.X == 2, axis=1)).reshape(-1)
        V = [(vi, 0) for vi in model.x] + [(model.x[k], 1) for k in der_col]
        X = model.X.copy()
        X[X == 3] = 1
        X = np.hstack((X, 1.0 * (X[:, der_col] == 2))).astype(np.int)
        X[X == 2] = 0
        return GraphPantelides(X, V)

    else:
        print('Pantelides only applicable to square systems')


# Define model

# x' = w
# y' = z
# w' = Tx
# z' = Ty-g
# 0  = x^2+y^2-L^2

model_def = {'type': 'Symbolic',
             'x': ['dx', 'dy', 'dw', 'dz', 'T', 'x', 'y', 'w', 'z'],
             'f': [], 'z': [],
             'parameters': ['g', 'L']}

sym.var(model_def['x'])
sym.var(model_def['parameters'])

model_def['rels'] = [
    -dx + w,
    -dy + z,
    -dw + T * x,
    -dz + T * y - g,
     x**2 + y**2 - L**2,
    fdt.DiffConstraint('dx','x'),
    fdt.DiffConstraint('dy', 'y'),
    fdt.DiffConstraint('dw', 'w'),
    fdt.DiffConstraint('dz', 'z')]

model = fdt.DiagnosisModel(model_def, name ='Pendulum')
s_idx, nu = model.Pantelides()
print(f"structural index = {s_idx}, nu = {nu}")




V2 = [('dx', 0), ('dy', 0), ('dw', 0), ('dz', 0), ('T', 0),
     ('x', 0), ('y', 0), ('w', 0), ('z', 0), ('x', 1), ('y', 1),
     ('w', 1), ('z', 1)]

#    d d d d           . . . .
#    x y w z T x y w z x y w z
X2 = np.array([[1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
              [0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
              [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
              [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
              [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
              [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]])



s_idx, nu = Pantelides(X, V)
print(f"Index = {s_idx}, nu = {nu}")

