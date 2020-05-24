import faultdiagnosistoolbox as fdt
import sympy as sym

modeldef = {}
modeldef['type'] = 'Symbolic'
modeldef['x'] = ['h1','h2','xc1','xl2','xl3','dh1','dh2']
modeldef['f'] = ['fh2','ff1','fc1','fl2','fl3','fa']
modeldef['z'] = ['y1','y2','y3','y4', 'u']
modeldef['parameters'] = ['d1','d2','d3','d4','d5','d6']

sym.var(modeldef['x'])
sym.var(modeldef['f'])
sym.var(modeldef['z'])
sym.var(modeldef['parameters'])

modeldef['rels'] = [
   -dh1 + d1 * u - d2 * xc1 * sym.sqrt(h1) + fa,
   -dh2 + d3 * xc1 * xl2 * sym.sqrt(h1) - d4 * sym.sqrt(h2),
   -y1 + h1,
   -y2 + h2 + fh2,
   -y3 + d5 * xc1 * sym.sqrt(h1) + ff1,
   -y4 + d6 * xl3 * sym.sqrt(h2),
   -xc1 + 1 - fc1,
   -xl2 + 1 - fl2,
   -xl3 + 1 - fl3,
   fdt.DiffConstraint('dh1', 'h1'),
   fdt.DiffConstraint('dh2', 'h2')]

model = fdt.DiagnosisModel(modeldef, name='Water Tank Model, reduced case')
