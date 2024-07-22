# %%
import dmpermlib
from scipy.sparse import csc_matrix
import numpy as np

# import faultdiagnosistoolbox as fdt
from faultdiagnosistoolbox.dmperm import dmperm, MSO

# %%
A = np.array([[1, 0], [0, 4], [5, 5], [0, 7]])
Asp = csc_matrix(A)

# %%
X = np.array(
    [
        [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0],
        [1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0],
        [0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0],
        [1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1],
        [2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0],
        [0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0],
        [0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0],
        [0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0],
        [0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0],
        [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    ]
)
Xsp = csc_matrix(X)
Xsp.indptr = Xsp.indptr.astype(np.int64)
Xsp.indices = Xsp.indices.astype(np.int64)
Xsp.data = Xsp.data.astype(np.float64)

MSO(Xsp)  # fdt version
dmperm(Xsp)  # fdt version
# %%

dmpermlib.MSO(Xsp)

# %%
dmpermlib.dmperm(Xsp)
