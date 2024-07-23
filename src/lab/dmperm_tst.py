# %%
import faultdiagnosistoolbox as fdt
from scipy.sparse import csc_matrix
import faultdiagnosistoolbox.structuralanalysis as sa
import faultdiagnosistoolbox.dmpermlib as dm
import numpy as np


# %%
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


# %%
model = fdt.models.induction_motor
X = model.X
Xsp = csc_matrix(X)
# Xsp.indptr = Xsp.indptr.astype(np.int64)
# Xsp.indices = Xsp.indices.astype(np.int64)
# Xsp.data = Xsp.data.astype(np.float64)
Xdict = CSCDict(Xsp)

# %% Old code
sa.dmperm_internal(Xdict)
sa.findmso_internal(Xdict)

# %% New code
print(dm.dmperm(Xsp))
# print(dm.MSO(Xsp))
