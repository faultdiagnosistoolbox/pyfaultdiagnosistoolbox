import faultdiagnosistoolbox as fdt
import sympy as sym
import matplotlib.pyplot as plt
import numpy as np

pendulum_def = {'type': 'Symbolic',
                'x': ['dx', 'dy', 'dw', 'dz', 'T', 'x', 'y', 'w', 'z'],
                'f': [], 'z': [],
                'parameters': ['g', 'L']}

sym.var(pendulum_def['x'])
sym.var(pendulum_def['parameters'])

pendulum_def['rels'] = [
    -dx + w,
    -dy + z,
    -dw + T * x,
    -dz + T * y - g,
    x**2 + y**2 - L**2,
    fdt.DiffConstraint('dx', 'x'),
    fdt.DiffConstraint('dy', 'y'),
    fdt.DiffConstraint('dw', 'w'),
    fdt.DiffConstraint('dz', 'z')]

pendulum = fdt.DiagnosisModel(pendulum_def, name='Pendulum')

sidx, nu = pendulum.Pantelides()
print(sidx)

# _, ax = plt.subplots(num=10, clear=True)
# pendulum.PlotModel(verbose=True)

# %%
X = pendulum.X.copy()
ne, nx = X.shape
der = np.zeros((ne, 1), dtype=int)
X[X == 3] = 1
X = X - 1
Xder = X + der @ np.ones((1, nx))
Xder[X < 0] = -1
hd = np.max(Xder, axis=0)
hod = (Xder == (np.ones((ne, 1)) @ hd.reshape((1, -1)))).astype(int)

p, q, r, s, cc, rr, m2 = fdt.dmperm.dmperm(hod)
print(m2)

# %% Extract matching
m = np.full(hod.shape[1], -1)
m[p[rr[0]:rr[1]]] = q[cc[1]:cc[2]]
m[p[rr[1]:rr[2]]] = q[cc[2]:cc[3]]
m[p[rr[2]:rr[3]]] = q[cc[3]:cc[4]]
print(m)
# sum(m >= 0)
