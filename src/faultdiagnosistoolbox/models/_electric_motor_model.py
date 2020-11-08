import faultdiagnosistoolbox as fdt
import sympy as sym

modelDef = {
    'type': 'Symbolic',
    'x': ['dI', 'dw', 'I', 'w', 'T', 'Tm', 'Tl'],
    'f': ['fR', 'fi', 'fw', 'fT'],
    'z': ['V', 'yi', 'yw', 'yT'],
    'parameters': ['Ka', 'b', 'R', 'J', 'L']}

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
    fdt.DiffConstraint('dI', 'I'),
    fdt.DiffConstraint('dw', 'w')]

model = fdt.DiagnosisModel(modelDef, name='Electric motor')
