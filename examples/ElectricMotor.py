import numpy as np
import matplotlib.pyplot as plt
import sys
new_paths = ['../Misc/', '../src/','../src/faultdiagnosistoolbox/']
[sys.path.append(d) for d in new_paths if not d in sys.path]
from misc import *
import faultdiagnosistoolbox as fdt
import sympy as sym


# %% Define model object
modelDef={}
modelDef['type'] = 'Symbolic'
modelDef['x'] = ['dI','dw','dth', 'I','w','th','alpha','DT','Tm','Tl']
modelDef['f'] = ['fR','fi','fw','fD']
modelDef['z']= ['V','yi','yw','yd']
modelDef['parameters'] = ['Ka','b','R','J','L']

sym.var(modelDef['x'])
sym.var(modelDef['f'])
sym.var(modelDef['z'])
sym.var(modelDef['parameters'])

modelDef['rels'] = [
  -V + I*(R+fR) + L*dI + Ka*I*w,
  -Tm + Ka*I**2,
  -J*dw + DT-b*w,
  -DT + Tm-Tl,
  -dth + w,
  -dw + alpha,
  -yi + I + fi,
  -yw + w + fw,
  -yd + DT + fD,
  fdt.DiffConstraint('dI','I'),
  fdt.DiffConstraint('dw','w'),
  fdt.DiffConstraint('dth','th')]

model = fdt.DiagnosisModel( modelDef, name='Electric Motor' )
model.Lint()


# %% Plot model structure
plt.figure(10)
model.PlotModel()

# %% Plot isolability properties and the Dulmage-Mendelsoh decomposition, including a canonical form of the over determined part

plt.figure(20)
model.IsolabilityAnalysis(plot=True, causality='der')
plt.title('Isolability with derivative causality'); 

plt.figure(21)
model.IsolabilityAnalysis(plot=True, causality='int')
plt.title('Isolability with integral causality'); 

plt.figure(22)
model.IsolabilityAnalysis(plot=True)
plt.title('Isolability with mixed causality'); 

plt.figure(24)
model.PlotDM(fault=True, eqclass=True)


# %% Residual generator design
# Compute set of MSO sets

msos = model.MSO()
print("Found " + str(len(msos)) + " MSO sets")

# %% Generate Python code for residual generators based on MSO1

mIdx = 0
print(model.MSOCausalitySweep(msos[mIdx]))
for k,red in enumerate(msos[mIdx]):
    Gamma = model.Matching([e for e in msos[mIdx] if e != red])
    model.SeqResGen(Gamma, red, 'EMSQResGen' + str(mIdx+1) + str(k+1))
    print("")


# %%Generate Python code for residual generators based on MSO2
mIdx = 1
print(model.MSOCausalitySweep(msos[mIdx]))
for k,red in enumerate(msos[mIdx]):
    Gamma = model.Matching([e for e in msos[mIdx] if e != red])
    model.SeqResGen(Gamma, red, 'EMSQResGen' + str(mIdx+1) + str(k+1))
    print("")


# %% Generate Python code for residual generators based on MSO3
mIdx = 2
print(model.MSOCausalitySweep(msos[mIdx]))
for k,red in enumerate(msos[mIdx]):
    Gamma = model.Matching([e for e in msos[mIdx] if e != red])
    model.SeqResGen(Gamma, red, 'EMSQResGen' + str(mIdx+1) + str(k+1))
    print("")

