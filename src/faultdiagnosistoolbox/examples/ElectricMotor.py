import matplotlib.pyplot as plt
import faultdiagnosistoolbox as fdt
import sympy as sym

# %% Define model object
model_def = {'type': 'Symbolic',
             'x': ['dIc', 'dw', 'dth', 'Ic', 'w', 'th', 'alpha', 'DT', 'Tm', 'Tl'],
             'f': ['fR', 'fi', 'fw', 'fD'],
             'z': ['V', 'yi', 'yw', 'yd'],
             'parameters': ['Ka', 'b', 'R', 'J', 'L']}

sym.var(model_def['x'])
sym.var(model_def['f'])
sym.var(model_def['z'])
sym.var(model_def['parameters'])

model_def['rels'] = [
  -V + Ic * (R + fR) + L * dIc + Ka * Ic * w,
  -Tm + Ka * Ic**2,
  -J * dw + DT - b * w,
  -DT + Tm - Tl,
  -dth + w,
  -dw + alpha,
  -yi + Ic + fi,
  -yw + w + fw,
  -yd + DT + fD,
  fdt.DiffConstraint('dIc','Ic'),
  fdt.DiffConstraint('dw','w'),
  fdt.DiffConstraint('dth','th')]

model = fdt.DiagnosisModel(model_def, name='Electric Motor')

model.Lint()

# Plot model structure
_, ax = plt.subplots(num=10)
model.PlotModel(ax=ax)

# Plot isolability properties and the Dulmage-Mendelsoh decomposition,
# including a canonical form of the over determined part
_, ax = plt.subplots(num=20)
model.IsolabilityAnalysis(ax=ax, causality='der')
ax.set_title('Isolability with derivative causality')

_, ax = plt.subplots(num=21)
model.IsolabilityAnalysis(ax=ax, causality='int')
ax.set_title('Isolability with integral causality')

_, ax = plt.subplots(num=22)
model.IsolabilityAnalysis(ax=ax)
ax.set_title('Isolability with mixed causality')

_, ax = plt.subplots(num=23)
model.PlotDM(fault=True, eqclass=True, ax=ax, verbose=True)

# Residual generator design
# Compute set of MSO sets
msos = model.MSO()
print(f"Found {len(msos)} MSO sets")

# Generate Python code for residual generators based on all MSOs
for i, mso_i in enumerate(msos):
    for k, red in enumerate(mso_i):
        print(f"EMSQResGen_{i + 1}_{k + 1}")
        Gamma = model.Matching([e for e in mso_i if e != red])
        model.SeqResGen(Gamma, red, f"EMSQResGen_{i + 1}_{k + 1}")
        print("")

plt.show()
