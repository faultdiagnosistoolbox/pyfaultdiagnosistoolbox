import numpy as np
import dmperm as dmperm
import matplotlib.pyplot as plt

class DiagnosisModel(object):
  X = np.zeros((0,0),dtype='float64')
  F = np.zeros((0,0),dtype='float64')
  Z = np.zeros((0,0),dtype='float64')
  x = []
  f = []
  z = []
  ne = 0;
  nx = 0;

  def __init__(self, modeldef) : 
    self.ne=len(modeldef['rels'])
    self.nx=len(modeldef['x'])
    self.X = self._ModelStructure( modeldef['rels'], modeldef['x'])
    self.x = modeldef['x']
    self.F = self._ModelStructure( modeldef['rels'], modeldef['f'])
    self.f = modeldef['f']
    self.Z = self._ModelStructure( modeldef['rels'], modeldef['z'])
    self.z = modeldef['z']
    self.e = map(lambda x:"e"+np.str(x+1),np.arange(0,self.ne))
    
  def GetDMParts(self):
    return dmperm.GetDMParts(self.X)
    
  def _ModelStructure(self,rels,x) : 
    ne=len(rels)
    nx=len(x)

    X = np.matrix(np.zeros((ne,nx)),dtype='float64')
    for k,rel in enumerate(rels) :
      X[k, self._RelModelStructure(rel, x)]=1
    return X

  def _RelModelStructure(self,rel,x):
    return [xi for xi in range(0,len(x)) if x[xi] in rel]

  def PlotDM(self) :
    
    dm = dmperm.GetDMParts(self.X)

    plt.spy(self.X[dm.rowp,:][:,dm.colp],markersize=2, marker="o")

    if len(dm.Mm.row)>0:
      r = len(dm.Mm.row);
      c = len(dm.Mm.col);
      x1 = 0.5;
      x2 = x1+c;
      y1 = 0.5;
      y2 = y1+r;
      plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')    

    # Plot exactly determined part
    r = 1+len(dm.Mm.row);
    c = 1+len(dm.Mm.col);
    for hc in dm.M0:
      n = len(hc.row);
      x1 = c-0.5;
      x2 = x1+n;
      y1 = r-0.5;
      y2 = y1+n;
      plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')
      r = r+n;
      c = c+n;
      
    # Plot over determined part  
    if len(dm.Mp.row)>0:
      nr = len(dm.Mp.row);
      nc = len(dm.Mp.col);
      x1 = c-0.5;
      x2 = x1+nc;
      y1 = r-0.5;
      y2 = y1+nr;
      plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')    

    # Plot axis ticks
    plt.xticks(np.arange(0,self.X.shape[0]),[self.x[i] for i in dm.colp])
    plt.yticks(np.arange(0,self.X.shape[1]),[self.e[i] for i in dm.rowp])
    
