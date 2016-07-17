import numpy as np
import dmperm as dmperm
import matplotlib.pyplot as plt
import scipy.sparse as sp

class DiagnosisModel(object):
  X = np.zeros((0,0),dtype='float64')
  F = np.zeros((0,0),dtype='float64')
  Z = np.zeros((0,0),dtype='float64')
  x = []
  f = []
  z = []
  name = ""

  def __init__(self, modeldef) : 
    self.X = self._ModelStructure( modeldef['rels'], modeldef['x'])
    self.x = modeldef['x']
    self.F = self._ModelStructure( modeldef['rels'], modeldef['f'])
    self.f = modeldef['f']
    self.Z = self._ModelStructure( modeldef['rels'], modeldef['z'])
    self.z = modeldef['z']
    self.e = map(lambda x:"e"+np.str(x+1),np.arange(0,self.ne()))

  def ne(self):
    return self.X.shape[0]

  def nx(self):
    return self.X.shape[1]

  def nf(self):
    return self.F.shape[1]

  def nz(self):
    return self.Z.shape[1]
    
  def GetDMParts(self):
    return dmperm.GetDMParts(self.X)
    
  def _ModelStructure(self,rels,x) : 
    ne=len(rels)
    nx=len(x)
    
    X = np.zeros((ne,nx),dtype='int64')
    for k,rel in enumerate(rels):
      if self._IsDifferentialConstraint(rel):
        if (rel[0] in x) and (rel[1] in x):
          dvIdx = x.index(rel[0])
          ivIdx = x.index(rel[1])
          X[k,dvIdx] = 3;
          X[k,ivIdx] = 2;
      else:        
        X[k, self._RelModelStructure(rel, x)]=1
    return sp.csc_matrix(X)

  def _RelModelStructure(self,rel,x):
    return [xi for xi in range(0,len(x)) if x[xi] in rel]

  def _IsDifferentialConstraint( self, rel ):
    return (len(rel)==3 and rel[2]=="diff");

  def MSO(self):
    return dmperm.MSO(self.X)
  
  def IsPSO( self, *args ):
    if len(args)>0:
        eq=args[0]
    else:
        eq=np.arange(0,self.X.shape[0])
        
    dm = dmperm.GetDMParts(self.X[eq,:])    
    return (len(dm.Mm.row)==0) and (len(dm.M0)==0)
  
  def PlotDM(self, **options) :
    labelVars = False;
    if options.has_key('verbose'):
      labelVars = options['verbose'];
    elif self.nx()+self.nf()+self.nz()<20:
      labelVars = True;
    
    dm = dmperm.GetDMParts(self.X)
    X=self.X[dm.rowp,:][:,dm.colp].todense()
    
    plt.spy(X==1,markersize=4, marker="o")
    
    for idx,val in enumerate(np.argwhere(X==3)):
      plt.text(val[1]-0.06,val[0]+0.15, 'I', color='b')
      
    for idx,val in enumerate(np.argwhere(X==2)):
      plt.text(val[1]-0.06,val[0]+0.15, 'D', color='b')

    if len(dm.Mm.row)>0:
      r = len(dm.Mm.row);
      c = len(dm.Mm.col);
      x1 = -0.5;
      x2 = x1+c;
      y1 = -0.5;
      y2 = y1+r;
      plt.plot( [x1, x1, x2, x2, x1],[y1, y2, y2, y1, y1],'b')    

    # Plot exactly determined part
    r = len(dm.Mm.row);
    c = len(dm.Mm.col);
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
    if labelVars:
      plt.xticks(np.arange(0,self.X.shape[1]),[self.x[i] for i in dm.colp])
      plt.yticks(np.arange(0,self.X.shape[0]),[self.e[i] for i in dm.rowp])

    # Change plot range
    plt.axis([-0.7,self.X.shape[1]-0.3,-0.7,self.X.shape[0]-0.3])

    plt.xlabel('Variables')
    plt.ylabel('Equations')

    
  def PlotModel(self, **options):

    labelVars = False;
    if options.has_key('verbose'):
      labelVars = options['verbose'];
    elif self.nx()+self.nf()+self.nz()<20:
      labelVars = True;
    
    X = self.X.todense()
    F = self.F.todense()
    Z = self.Z.todense()
    
    plt.spy(np.concatenate((X==1,np.zeros(F.shape),np.zeros(Z.shape)),axis=1),markersize=4,marker="o", color="b")
    plt.spy(np.concatenate((np.zeros(X.shape),F,np.zeros(Z.shape)),axis=1),markersize=4,marker="o",color="r")
    plt.spy(np.concatenate((np.zeros(X.shape),np.zeros(F.shape),Z),axis=1),markersize=4,marker="o",color="k")
    
    for idx,val in enumerate(np.argwhere(X==3)):
      plt.text(val[1]-0.06,val[0]+0.15, 'I',color="b")
      
    for idx,val in enumerate(np.argwhere(X==2)):
      plt.text(val[1]-0.06,val[0]+0.15, 'D',color="b")

    # Plot axis ticks
    if labelVars:
      plt.xticks(np.arange(0,self.nx()+self.nf()+self.nz()),self.x + self.f + self.z)
      plt.yticks(np.arange(0,self.X.shape[0]),self.e)

    # Plot variable set divisors
    plt.plot( [self.nx()-1+0.5, self.nx()-1+.5], [0, self.ne()-1 ], color="k", linestyle="dashed")
    plt.plot( [self.nx()+self.nf()-1+0.5, self.nx()+self.nf()-1+.5], [0, self.ne()-1 ], color="k", linestyle="dashed")
    
    # Change plot range
    #plt.axis([-0.7,self.X.shape[1]-0.3,self.X.shape[0]-0.3, -0.7])

    if len(self.name)>0:
      plt.title(self.name)

    plt.xlabel('Variables')
    plt.ylabel('Equations')
      

def DiffConstraint(dvar,ivar):
  return [dvar, ivar, "diff"];

