import faultdiagnosistoolbox as fdt
modelDef={}

modelDef['type'] = 'VarStruc'
modelDef['x']= ['p1','p2','p3','q0','q1','q2','q3','dp1','dp2','dp3']
modelDef['f'] = ['fV1','fV2','fV3','fT1','fT2','fT3']
modelDef['z'] = ['y1','y2','y3']
modelDef['parameters'] = ['Rv1', 'Rv2', 'Rv3', 'CT1', 'CT2', 'CT3']
modelDef['rels'] = [
    ['q1', 'p1', 'p2', 'fV1'],
    ['q2', 'p2', 'p3', 'fV2'],
    ['q3', 'p3', 'fV3'],
    ['dp1','q0','q1','fT1'],
    ['dp2','q1','q2','fT2'],
    ['dp3','q2','q3','fT3'],
    ['y1','p1'],
    ['y2','q2'],
    ['y3','q0'],
    fdt.DiffConstraint('dp1','p1'),
    fdt.DiffConstraint('dp2','p2'),
    fdt.DiffConstraint('dp3','p3')]

model = fdt.DiagnosisModel( modelDef, name='Three Tank System' )
