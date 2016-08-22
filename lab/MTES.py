
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib qt5')
import matplotlib
import matplotlib.pyplot as plt
import sys
import numpy as np
import sympy as sym
import scipy.sparse as sp
new_paths = ['../Misc/', '../src/']
[sys.path.append(d) for d in new_paths if not d in sys.path];
from misc import *
import faultdiagnosistoolbox as fdt


# In[3]:

modelDef = {}
modelDef['type'] = 'VarStruc'
modelDef['x'] = ['p0', 'p1', 'p2', 'p3', 'pf', 'q0','q1', 'q2', 'q3', 'qp'] # 'Rv0'};
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


# In[4]:

len(model.MSO())


# In[10]:

map(lambda fi: np.where(model.F.toarray()[:,fi])[0][0], np.arange(0,model.F.shape[1]))


# In[ ]:
def MTES( self ):
    fidx = list(map(lambda fi: np.argwhere(fi)[0][0], np.transpose(self.F)))
    return MTES_TEsub(self.X, fidx,0)
    
def MTES_TEsub(X,f,p):
    S = {'eq':[], 'f':[], 'sr':[]}    
    m = MTES_initModel(X,f) # overdetermined or empty
    if m['sr']>0 and len(m['f'])>0:
        S = MTES_FindMTES(m,p)
    return S
    
def MTES_initModel(X,f):
    dm = fdt.GetDMParts(X)
    row_over = dm.Mp.row
    col_over = dm.Mp.col
    sr = len(row_over) - len(col_over)
    
    
#    m ={'sr':}

# %%
fidx = []
for ef in model.F:
    fi=np.argwhere(ef)
    if len(fi)>0:
        fidx.append(fi[0])
    else:
        fidx.append(np.array([],dtype=np.int64))
fidx = np.array(fidx)
                
X = model.X
dm = fdt.GetDMParts(X)
row_over = dm.Mp.row
col_over = dm.Mp.col
sr = len(row_over) - len(col_over)
f = fidx
f = [ef for ef in fidx if ef in row_over]

