import numpy as np

class DiagnosisModel(object):
  X = np.zeros((0,0))
  F = np.zeros((0,0))
  Z = np.zeros((0,0))
  x = []
  f = []
  z = []
  ne = 0;
  nx = 0;

  def __init__(self, modeldef) : 
    self.ne=len(modeldef['rels'])
    self.nx=len(modeldef['x'])
    self.X = ModelStructure( modeldef['rels'], modeldef['x'])
    self.x = modeldef['x']
    self.F = ModelStructure( modeldef['rels'], modeldef['f'])
    self.f = modeldef['f']
    self.Z = ModelStructure( modeldef['rels'], modeldef['z'])
    self.z = modeldef['z']

#  def PlotModel(self) :
#    return plt.spy(self.X, markersize=2, precision=0)

def ModelStructure(rels,x) : 
    ne=len(rels)
    nx=len(x)

    X = np.matrix(np.zeros((ne,nx)))
    for k,rel in enumerate(rels) :
      X[k, RelModelStructure(rel, x)]=1
    return X

def RelModelStructure(rel,x):
  return [xi for xi in range(0,len(x)) if x[xi] in rel]

