
# coding: utf-8

# In[1]:
import sys
import numpy as np
new_paths = ['../Misc/', '../src/']
[sys.path.append(d) for d in new_paths if not d in sys.path];
import faultdiagnosistoolbox as fdt

# In[3]:

modelDef = {}
modelDef['type'] = 'VarStruc'
modelDef['x'] = ['p0', 'p1', 'p2', 'p3', 'pf', 'q0','q1', 'q2', 'q3', 'qp','Rv0'];
modelDef['f'] = ['fC1','fC2','fC3','fv1','fv2','fv3','fpipe']
modelDef['z'] = ['y1','y2','y3','y4']
modelDef['rels'] = [
  ['qp','p0','fpipe'],
  ['qp','q0'],
  ['pf','p0','p1'],
  ['q0','pf', 'Rv0'],
  ['q1','p1','p2','fv1'],
  ['q2','p2','p3','fv2'],
  ['q3','p3','fv3'],
  ['p1','q0','q1','fC1'],
  ['p2','q1','q2','fC2'],
  ['p3','q2','q3','fC3'],
  ['y1','q0'],
  ['y2','p1'],
  ['y3','q3'],
  ['y4','p3'],
  ['Rv0','y2']]
model = fdt.DiagnosisModel(modelDef, name='MTES example')

# %%
modelDef={}
modelDef['type'] = 'VarStruc'
modelDef['x'] = ['x1','x2','x3']
modelDef['f'] = ['f1','f2','f3','f4','f5']
modelDef['z'] = ['y1','y2','y3', 'u']
modelDef['rels'] = [
  ['x1','u','f1'],
  ['x2','x1','x3','f2'],
  ['x3','x2'],
  ['y1','x2','f3'],
  ['y2','x2','f4'],
  ['y3','x3','f5']]
model = fdt.DiagnosisModel(modelDef, name='MTES small example')

# %%
m = MTES_initModel(model)
row = 0

n,nrow = MTES_LumpExt(m,row)



# In[ ]:
def MTES( self ):
    S = {'eq':[], 'f':[], 'sr':[]}    
    m = MTES_initModel(self) # overdetermined or empty
    if m['sr']>0 and len(m['f'])>0:
        S = MTES_FindMTES(m,0)
    return S
    
def MTES_initModel(model):
    dm = fdt.GetDMParts(model.X)
    row_over = dm.Mp.row
    col_over = dm.Mp.col

    ef = np.argwhere(np.any(model.F,axis=1)).flatten()
    ef = np.intersect1d(ef,row_over)
    idx = np.hstack((ef,np.setdiff1d(row_over,ef)))
    m = {}
    m['sr'] = len(row_over) - len(col_over)
    m['X'] = model.X[idx,:][:,col_over]
    m['f'] = list(map(lambda fi: np.argwhere(fi)[0], model.F[ef,:]))
    m['delrow']=0
    m['e'] = list(map(lambda ei: np.array([ei]), idx))
    return m
    
def MTES_FindMTES(m,p):
    return []

def MTES_LumpExt(m,row):
    n = {}

    no_rows = m['X'].shape[0]
    remRows = np.hstack((np.arange(0,row-1),np.arange(row+1,no_rows)))
    remRowsf = np.hstack((np.arange(0,row-1),np.arange(row+1,len(m['f']))))
  
    dm = fdt.GetDMParts(m['X'][remRows,:])
    row_just = dm.M0eqs
    row_over = dm.Mp.row
    col_over = dm.Mp.col
    if len(row_just)>0:
        eqcls = np.hstack((remRows[row_just],[row]))
        no_rows_before_row = np.sum(eqcls<row)
        row = row - no_rows_before_row
        no_rows_before = np.sum(eqcls<m['delrow'])
        n['delrow'] = m['delrow'] - no_rows_before
        eqclsf = eqcls[eqcls<len(m['f'])]
        row_overf = row_over[row_over<len(remRowsf)]
        if no_rows_before > 0:
            rowinsert = n['delrow']
        else:
            rowinsert = row

        n['X'] = np.vstack((m['X'][remRows[row_over[0:rowinsert]],:][:,col_over],
                        np.any(m['X'][eqcls,:][:,col_over],axis=0).astype(np.int64),
                        m['X'][remRows[row_over[rowinsert:]],:][:,col_over]))

        ef1 = list(np.array(m['f'])[remRowsf[row_overf[0:rowinsert]]])
        ef2 = list([np.array(m['f'])[eqclsf].flatten()])
        ef3 = list(np.array(m['f'])[remRowsf[row_overf[rowinsert:]]])
        n['f'] = ef1 + ef2 + ef3

        e1 = list(np.array(m['e'])[remRows[row_over[0:rowinsert]]])
        e2 = list([np.array(m['e'])[eqcls].flatten()])
        e3 = list(np.array(m['e'])[remRows[row_over[rowinsert:]]])
        n['e'] = e1 + e2 + e3
              
        n['sr'] = m['sr']

        if no_rows_before > 0: 
            n['delrow']= n['delrow']+1
    else:
        n = m

    row = row + 1
    return (n,row)

