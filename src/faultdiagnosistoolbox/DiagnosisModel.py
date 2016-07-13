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
    self.X = self.ModelStructure( modeldef['rels'], modeldef['x'])
    self.x = modeldef['x']
    self.F = self.ModelStructure( modeldef['rels'], modeldef['f'])
    self.f = modeldef['f']
    self.Z = self.ModelStructure( modeldef['rels'], modeldef['z'])
    self.z = modeldef['z']

  def ModelStructure(self,rels,x) : 
    ne=len(rels)
    nx=len(x)

    X = np.matrix(np.zeros((ne,nx)))
    for k,rel in enumerate(rels) :
      X[k, self.RelModelStructure(rel, x)]=1
    return X

  def RelModelStructure(self,rel,x):
    return [xi for xi in range(0,len(x)) if x[xi] in rel]

