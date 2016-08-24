
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
m,_ = MTES_LumpExt(m,row)
if len(m['f'])==1: # m is an MTES
    S = MTES_storeFS(m)
else: #recurse
    if p==0:
        S=MTES_storeFS(m)
    else:
        S = {'eq':[], 'f':[], 'sr':0}
    row=m['delrow']
    while len(m['f'])>=row:
        m,row=MTES_LumpExt(m,row)
    for delrow in np.arange(m['delrow'],len(m['f'])):
        # create the model where delrow has been removed
        m['delrow']=delrow
        rows = np.delete(np.arange(0,m['X'].shape[0]),delrow)
        n = MTES_GetPartialModel(m,rows)
        Sn = MTES_FindMTES(n,p) # make recursive call
        S = MTES_addResults(S,Sn)
        
# %%
def MTES_GetPartialModel(m,rows):
    n = {}
    vars = np.any(m['X'][rows,:],axis=0)
    n['X'] = m['X'][rows,:][:,vars]
    n['e'] = list(np.array(n['e'])[rows])
[ei for ei in rows if ei < len(m['f'])]
    return n    


    

# %%
function n = GetPartialModel(m,rows)
  n.sm = m.sm(rows,any(m.sm(rows,:),1));
  n.e =  m.e(rows);
  n.f = m.f(intersect(rows,1:length(m.f)));
  n.sr = size(n.sm,1)-size(n.sm,2);
  n.delrow = length(find(rows<m.delrow))+1;
end

function S = FindMTES(m,p)
  m = LumpExt(m,1);
  if length(m.f)==1 % if m is MTES
    S = storeFS(m); % then store m
  else %otherwise make recursive call
    if p == 1
      S = storeFS(m);
    else
      S.eq = {};
      S.f = {};
      S.sr = [];
    end
    row = m.delrow;
    while length(m.f)>=row % some rows are allowed to be removed
      [m,row] = LumpExt(m,row); % lump model w.r.t. row
    end
    for delrow = m.delrow:length(m.f)
      % create the model where delrow has been removed
      m.delrow = delrow;
      rows = [1:delrow-1 delrow+1:size(m.sm,1)];
      n = GetPartialModel(m,rows);

      Sn = FindMTES(n,p); % make recursive call
      S = addResults(S,Sn); % store results
    end   
  end
end


# %% Code
def MTES_storeFS(m):
    eq = list(np.sort(np.concatenate(m['e'])))
    f = list(np.sort(np.concatenate(m['f'])))
    return {'eq': eq, 'f': f, 'sr':m['sr']}
            
def MTES_FindMTES(m,p):
    pass
    
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

