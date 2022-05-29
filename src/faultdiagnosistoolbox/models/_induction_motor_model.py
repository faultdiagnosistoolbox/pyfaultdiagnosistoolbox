import faultdiagnosistoolbox as fdt
import sympy as sym
modelDef = {}
modelDef['type'] = 'Symbolic'
modelDef['x'] = ['i_a', 'i_b', 'lambda_a', 'lambda_b', 'w',
                 'di_a', 'di_b', 'dlambda_a', 'dlambda_b', 'dw', 'q_a', 'q_b', 'Tl']
modelDef['f'] = ['f_a', 'f_b']
modelDef['z'] = ['u_a', 'u_b', 'y1', 'y2', 'y3']
modelDef['parameters'] = ['a', 'b', 'c', 'd', 'L_M', 'k', 'c_f', 'c_t']


(i_a, i_b, lambda_a, lambda_b, w,
 di_a, di_b, dlambda_a, dlambda_b, dw, q_a, q_b, Tl) = sym.symbols(modelDef['x'])
f_a, f_b = sym.symbols(modelDef["f"])
u_a, u_b, y1, y2, y3 = sym.symbols(modelDef["z"])
a, b, c, d, L_M, k, c_f, c_t = sym.symbols(modelDef["parameters"])

modelDef['rels'] = [
    -q_a + w * lambda_a,
    -q_b + w * lambda_b,
    -di_a + -a * i_a + b * c * lambda_a + b * q_b + d * u_a,
    -di_b + -a * i_b + b * c * lambda_b + b * q_a + d * u_b,
    -dlambda_a + L_M * c * i_a - c * lambda_a - q_b,
    -dlambda_b + L_M * c * i_b - c * lambda_b - q_a,
    -dw + -k * c_f * w + k * c_t * (i_a * lambda_b - i_b * lambda_a) - k * Tl,
    fdt.DiffConstraint('di_a', 'i_a'),
    fdt.DiffConstraint('di_b', 'i_b'),
    fdt.DiffConstraint('dlambda_a', 'lambda_a'),
    fdt.DiffConstraint('dlambda_b', 'lambda_b'),
    fdt.DiffConstraint('dw', 'w'),
    -y1 + i_a + f_a,
    -y2 + i_b + f_b,
    -y3 + w]

model = fdt.DiagnosisModel(modelDef, name='Induction motor')
