import faultdiagnosistoolbox as fdt
import sympy as sym


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

model = fdt.DiagnosisModel(modelDef, name = 'Electric motor')