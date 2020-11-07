# %load_ext autoreload
# %autoreload 2
import faultdiagnosistoolbox as fdt
import sympy as sym
import numpy as np

modelDef = {}
modelDef['type'] = 'Symbolic'
modelDef['x'] = ['dI','dw','dth', 'I','w','th','alpha','T','Tm','Tl']
modelDef['f'] = ['fR','fi','fw','fT']
modelDef['z']= ['V','yi','yw','yT']
modelDef['parameters'] = ['Ka','b','R','J','L']

sym.var(modelDef['x'])
sym.var(modelDef['f'])
sym.var(modelDef['z'])
sym.var(modelDef['parameters'])

modelDef['rels'] = [
    -V + I*(R + fR) + L*dI + Ka*I*w,
    -Tm + Ka*I**2,
    -J*dw + T - b*w,
    -T + Tm - Tl,
    -yi + I + fi,
    -yw + w + fw,
    -yT + T + fT,
    fdt.DiffConstraint('dI','I'),
    fdt.DiffConstraint('dw','w')]

model = fdt.DiagnosisModel(modelDef, name = 'Electric motor')

msos = model.MSO()
mso = np.sort(msos[2])

# Integral
red = mso[2]
M0 = [ei for ei in mso if ei != red]
Gamma_int = model.Matching(M0)
model.SeqResGen(Gamma_int, red, 'r_int')

# Derivative
red = mso[1]
M0 = [ei for ei in mso if ei != red]
Gamma_der = model.Matching(M0)
model.SeqResGen(Gamma_der, red, 'r_der')

# Mixed
red = mso[0]
M0 = [ei for ei in mso if ei != red]
Gamma_mixed = model.Matching(M0)

model.SeqResGen(Gamma_mixed, red, 'r_mixed')
