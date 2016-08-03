import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import sympy as sym
import copy
import dmperm as dmperm
import Matching as match
import StructurePlotting as smplot
import CodeGeneration as codegen
import SensorPlacement as sensplace
from VarIdGen import VarIdGen
import sys

class DiagnosisModel(object):
    def __init__(self, modeldef, name='') : 
        self.X = np.zeros((0,0),dtype='float64')
        self.F = np.zeros((0,0),dtype='float64')
        self.Z = np.zeros((0,0),dtype='float64')
        self.x = []
        self.f = []
        self.z = []
        self.rels = []
        self.name = name
        self.parameters = []

        self.vGen = VarIdGen()
        
        if modeldef['type'] is 'VarStruc' or modeldef['type'] is 'Symbolic':
            self.X = _ModelStructure( modeldef['rels'], modeldef['x'])
            self.x = modeldef['x']
            self.F = _ModelStructure( modeldef['rels'], modeldef['f'])
            self.f = modeldef['f']
            self.Z = _ModelStructure( modeldef['rels'], modeldef['z'])
            self.z = modeldef['z']
            self.e = map(lambda x:self.vGen.NewE(),np.arange(0,self.ne()))
            self.modelType = modeldef['type']

            if 'parameters' in modeldef:
                self.parameters = modeldef['parameters']
        elif modeldef['type'] is 'MatrixStruc':
            self.X = sp.csc_matrix(modeldef['X'])
            ne = self.X.shape[0]
            if len(modeldef['F'])>0:
                self.F = sp.csc_matrix(modeldef['F'])
            else:
                self.F = sp.csc_matrix(np.zeros((ne,0),dtype=np.int64))
                
            if len(modeldef['Z'])>0:
                self.Z = sp.csc_matrix(modeldef['Z'])
            else:
                self.Z = sp.csc_matrix(np.zeros((ne,0),dtype=np.int64))

            if modeldef.has_key('x'):
                self.x = modeldef['x']
            else:
                self.x = map(lambda x:"x"+np.str(x+1),np.arange(0,self.X.shape[1]))
            if modeldef.has_key('f'):
                self.f = modeldef['f']
            else:
                self.f = map(lambda x:"f"+np.str(x+1),np.arange(0,self.F.shape[1]))
            if modeldef.has_key('z'):
                self.z = modeldef['z']
            else:
                self.z = map(lambda x:"z"+np.str(x+1),np.arange(0,self.Z.shape[1]))

            self.e = map(lambda x:"e"+np.str(x+1),np.arange(0,self.ne()))
            self.modelType = modeldef['type']
        else:
            print 'Model definition type ' + modeldef['type'] + ' is not supported (yet)'

        if modeldef['type'] is 'Symbolic':
            self.syme = np.array(_ToEquations(modeldef['rels']))

        self.P = np.arange(0,len(self.x))
        self.Pfault = []

    def copy(self):
        return copy.deepcopy(self)
    
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
        
    def MSO(self):
        return dmperm.MSO(self.X)

    def sprank(self):
        return dmperm.sprank(self.X)
  
    def IsPSO( self, *args ):
        if len(args)>0:
            eq=args[0]
        else:
            eq=np.arange(0,self.X.shape[0])
        
        dm = dmperm.GetDMParts(self.X[eq,:])    
        return (len(dm.Mm.row)==0) and (len(dm.M0)==0)
  
    def IsHighIndex(self, **opts):
        if opts.has_key('eq'):
            eq = opts['eq']
        else:
            eq = np.arange(0,self.X.shape[0])
        return dmperm.IsHighIndex(self.X,eq=eq)

    def IsLowIndex(self, **opts):
        if opts.has_key('eq'):
            eq = opts['eq']
        else:
            eq = np.arange(0,self.X.shape[0])
        return dmperm.IsLowIndex(self.X,eq=eq)

    def IsUnderdetermined(self):
        dm = dmperm.GetDMParts(self.X)
        return len(dm.Mm.row)>0

    def FSM(self, eqs_sets, plot=False ):
        r=np.zeros((len(eqs_sets),self.F.shape[1]),dtype=np.int64)
        for idx,eqs in enumerate(eqs_sets):
            r[idx,:] = np.any(self.F[eqs,:].todense(),axis=0)

        if plot:
            plt.spy(r, markersize=10, marker='o')
            plt.xticks(np.arange(0,self.nf()),self.f)
            plt.yticks(np.arange(0,len(eqs_sets)), ["eq. set "+str(k+1) for k in np.arange(0,len(eqs_sets))])
            plt.gca().xaxis.tick_bottom()

        return r

    def Matching(self,eqs):
        return match.Matching(self.X.todense(),np.array(eqs))

    def MSOCausalitySweep(self, mso, diffres='int', causality=''):
        def IsDifferentialStructure( X, eq ):
            return np.any(X[eq,:]==3)

        if (causality is 'der') and (not diffres is 'der'):
            diffres = 'der'
        if (causality is 'int') and (not diffres is 'int'):
            diffres = 'int'
        res = []
        X = self.X.todense()
        for red in mso:
            m0 = np.sort([e for e in mso if e != red])
            Gamma = self.Matching(m0)
        
            if not IsDifferentialStructure(X,red):
                res.append(Gamma.matchType)
            elif diffres is 'int':
                if Gamma.matchType is 'der':
                    res.append('mixed')
                elif Gamma.matchType is 'int' or Gamma.matchType is 'mixed':
                    res.append(Gamma.matchType)
                elif Gamma.matchType is 'algebraic':
                    res.append('int')
            elif diffres is 'der':
                if Gamma.matchType is 'int':
                    res.append('mixed')
                elif Gamma.matchType is 'der' or Gamma.matchType is 'mixed':
                    res.append(Gamma.matchType)
                elif Gamma.matchType is 'algebraic':
                    res.append('der')
        if len(causality)>0:
            res = np.array(map(lambda c: c is causality or c is 'algebraic', res))
        return res
    
    def IsDynamic(self):
        return np.any(self.X.todense()==2)

    def IsStatic(self):
        return not self.IsDynamic()

    def Redundancy(self, eqs=''):
        if len(eqs)==0:
            eqs = np.arange(0,self.ne())
        
        dm = dmperm.GetDMParts( self.X[eqs,:] )
        return len(dm.Mp.row)-len(dm.Mp.col)

    def MTESRedundancy( self ):
        eqs = np.argwhere(np.any(self.F.todense(),axis=1))[:,0]
        return self.Redundancy(eqs) + 1
    
    def PlotDM(self, **options) :
        labelVars = False
        if options.has_key('verbose'):
            labelVars = options['verbose']
        elif self.X.shape[0]<40:
            labelVars = True
        if options.has_key('eqclass'):
            eqclass=options['eqclass']
        else:
            eqclass=False
        if options.has_key('fault'):
            fault=options['fault']
        else:
            fault=False
        smplot.PlotDM(self, verbose=labelVars, eqclass=eqclass, fault=fault)
            
    def PlotModel(self, **options):
        labelVars = False;
        if options.has_key('verbose'):
            labelVars = options['verbose'];
        elif self.nx()+self.nf()+self.nz()<40:
            labelVars = True;

        smplot.PlotModel(self, verbose=labelVars)
        
    def PlotMatching( self, Gamma, **options):
        q = np.array([],dtype=np.int64)
        for g in Gamma.matching:
            q = np.concatenate((q,g.col))
        
        # Determine if axis should be labeled
        labelVars = False;
        if options.has_key('verbose'):
            labelVars = options['verbose'];
        elif len(q)<40:
            labelVars = True;
        smplot.PlotMatching(self, Gamma, verbose=labelVars)
        
    def PossibleSensorLocations(self, x=-1):

        if x==-1:
            self.P = np.arange(0,len(self.x)) # Assume all possible sensor locations
        else:
            if issubclass(type(x[0]),str): # list of strings
                self.P = np.array([self.x.index(xi) for xi in x if xi in self.x])
            else:
                self.P = x.copy()

    def SensorLocationsWithFaults(self, x=[]):
        if len(x)>0:
            if issubclass(type(x[0]),str): # list of strings
                self.Pfault = np.array([self.x.index(xi) for xi in x if xi in self.x])
            else:
                self.Pfault = x.copy()
        else:
            self.Pfault = []

    def SensorPlacementIsolability(self):
        return sensplace.SensorPlacementIsolability(self)

    def SensorPlacementDetectability(self):
        return sensplace.SensorPlacementDetectability(self)

    def AddSensors(self,sensors,name=[],fault=[]):
        if issubclass(type(sensors[0]),str): # list of strings, convert to indices into self.x
            s = np.array([self.x.index(xi) for xi in sensors if xi in self.x])
        else:
            s = sensors
            
        ns = len(s)
        nx = self.X.shape[1]
        nz = self.Z.shape[1]
        nf = self.F.shape[1]
        ne = self.X.shape[0]
        Xs = np.zeros((ns,nx),np.int64)
        Fs = np.zeros((ns,nf),np.int64)
        Zs = np.zeros((ns,nz+ns),np.int64)
    
        fs = np.zeros(ns).astype(np.bool)
        for sIdx, si in enumerate(s):
            Xs[sIdx,si] = 1
            Zs[sIdx,sIdx+nz] = 1
            if si in self.Pfault:
                nf = nf+1
                Fs = np.hstack((Fs,np.zeros((ns,1),dtype=np.int64))) # Add column for new fault
                Fs[sIdx,-1] = 1 # Add fault
                fs[sIdx] = True
            else:
                Fs[sIdx,:] = np.zeros((1,nf))
    
        self.X = sp.csc_matrix(np.vstack((self.X.toarray(),Xs)))
        self.Z = np.hstack((self.Z.toarray(),np.zeros((ne,ns),dtype=np.int64)))
        self.Z = sp.csc_matrix(np.vstack((self.Z,Zs)))
        self.F = self.F.toarray()
        if np.sum(fs)>0:
            self.F = np.hstack((self.F,np.zeros((ne,np.sum(fs)),dtype=np.int64)))
        self.F = sp.csc_matrix(np.vstack((self.F,Fs)))
    
        self.e = self.e + map(lambda x:self.vGen.NewE(),s)    
        
        for idx,zi in enumerate(s):
            if len(name)==0:
                znum = np.sum(np.array(s)[0:idx] == s[idx])+1
                if znum>1:
                    zName = "z" + str(znum) + self.x[zi]
                else:
                    zName = "z" + self.x[zi]
            else:
                zName = name[idx]
            self.z = self.z + [zName]
                
            if fs[idx]:
                if len(fault)==0:
                    fName = 'f' + zName
                else:
                    fName = fault[idx]
                self.f = self.f + [fName]
                
            if self.modelType is 'Symbolic':
                if fs[idx]:
                    rel = sym.Eq(sym.symbols(zName),sym.symbols(self.x[zi])+sym.symbols(fName))
                else:
                    rel = sym.Eq(sym.symbols(zName),sym.symbols(self.x[zi]))
                self.syme = np.concatenate((self.syme,[rel]))

    def DetectabilityAnalysis(self):
        dm = dmperm.GetDMParts(self.X)
        df = [self.f[fidx] for fidx in np.arange(0,self.F.shape[1]) if np.argwhere(self.F[:,fidx])[:,0] in dm.Mp.row]
        ndf = [self.f[fidx] for fidx in np.arange(0,self.F.shape[1]) if np.argwhere(self.F[:,fidx])[:,0] not in dm.Mp.row]
        return (df,ndf)
    
    def IsolabilityAnalysis( self, plot=False, permute=True, causality='mixed' ):

        MplusCausal = lambda X: dmperm.Mplus(X,causality=causality)
        ne = self.X.shape[0]
        nf = len(self.f)

        plusRow = MplusCausal( self.X )

        # Determine equations for each fault
        feq = map(lambda fi: np.argwhere(self.F[:,fi])[0][0],np.arange(0,nf))

        # Determine non-detectable faults
        ndrows = [x for x in np.arange(0,ne) if x not in plusRow]
        ndf = [self.f[fi] for fi in np.arange(0,len(self.f)) if feq[fi] in ndrows]        
        df = [self.f[fi] for fi in np.arange(0,len(self.f)) if not feq[fi] in ndrows]

        im = np.ones((nf,nf),dtype=np.int64)
        for fi in np.arange(0,nf):
            # Decouple fi
            fieqs = [x[0] for x in np.argwhere(self.F[:,fi].todense()==0)]
            X = self.X[fieqs,:]
            plusRow = MplusCausal(X)
            fisolrows = [fieqs[ei] for ei in plusRow]
            idx = [fj for fj in np.arange(0,nf) if np.any(self.F[fisolrows,:].todense(),axis=0)[(0,fj)]]
            im[idx,fi]=0;

        if plot:
            if permute:
                p,q,_,_,_,_ = dmperm.dmperm(im)
            else:
                p = np.arange(0,nf)
                q = p
            plt.spy(im[p,:][:,q],markersize=8,marker="o", color="b")
            plt.xticks(np.arange(0,nf),[self.f[fi] for fi in p])
            plt.yticks(np.arange(0,nf),[self.f[fi] for fi in p])
            if len(self.name)>0:
                titleString = "Isolability matrix for '%s'" % self.name
            else:
                titleString = "Isolability matrix"
            if causality is 'der':
                titleString = "%s (derivative causality)" % titleString
            elif causality is 'int':
                titleString = "%s (integral causality)" % titleString
            plt.title( titleString )
            plt.gca().xaxis.tick_bottom()
        return im        

    def IsolabilityAnalysisArrs(self,arrs,permute=True, plot=False):
        FSM = self.FSM(arrs)
        return self.IsolabilityAnalysisFSM(FSM,permute=permute,plot=plot)

    def IsolabilityAnalysisFSM(self,FSM,permute=True, plot=False):
        nf = FSM.shape[1]
        nr = FSM.shape[0]
        im = np.ones((nf,nf),dtype=np.int64)
        for f in FSM:
            zIdx = np.array([[x0,y0] for x0 in np.where(f>0)[0] for y0 in np.where(f==0)[0]])
            if len(zIdx)>0:
                im[zIdx[:,0],zIdx[:,1]]=0

        if plot:
            if permute:
                p,q,_,_,_,_ = dmperm.dmperm(im);
            else:
                p = np.arange(0,nf)
                q = p
            plt.spy(im[p,:][:,q],markersize=10,marker='o')
            plt.xticks(np.arange(0,self.nf()),self.f)
            plt.yticks(np.arange(0,self.nf()),self.f)
            plt.gca().xaxis.tick_bottom()
            if len(self.name)>0:
                plt.title("Isolability matrix for a given FSM in '" + self.name + "'")
            else:
                plt.title('Isolability matrix for a given FSM')
        return im 

    def SeqResGen(self, Gamma, resEq, name, diffres='int', language='Python', batch=False, external=[], api='Python'):
        codegen.SeqResGen(self, Gamma, resEq, name, diffres=diffres, language=language, batch=batch, external=external)

    def Lint(self):
        war = False
        err = False

        dm = dmperm.GetDMParts(self.X)

        if len(self.name)>0:
            print 'Model: ' + self.name
        else:
            print "Model information"

        sys.stdout.write("\n  Type:" + self.modelType)

        nd = np.sum(self.X==3)
        if nd>0:
            sys.stdout.write(", dynamic\n")
        else:
            sys.stdout.write(", static\n")

        print '\n  Variables and equations'
        print '    ' + str(self.nx()) + ' unknown variables'
        print '    ' + str(self.nz()) + ' known variables'
        print '    ' + str(self.nf()) + ' fault variables'
        print '    ' + str(self.ne()) + ' equations, including ' + str(nd) + ' differential constraints'
        print '\n  Degree of redundancy: ' + str(self.Redundancy())

        if self.Redundancy()>0 and len(self.f)>0:
            print '  Degree of redundancy of MTES set: ' + str(self.MTESRedundancy())
        print '\n'


        if self.ne() != self.F.shape[0] or self.ne() != self.Z.shape[0]:
            print 'Error! Inconsistent numnber of rows in incidence matrices'
            err = err+1

        if self.nx() != len(self.x):
            print 'Error! Inconsistent number of unknown variables'
            err = err+1

        if self.nz() != len(self.z):
            print 'Error! Inconsistent number of known variables'
            err = err+1

        if self.nf() != len(self.f):
            print 'Error! Inconsistent number of fault variables'
            err = err+1

        if self.ne()!= len(self.e):
            print 'Error! Inconsistent number of equations'
            err = err+1

        if len([v for v in self.P if not v in np.arange(0,self.nx())])>0:
            print 'Error! Possible sensor locations outside set of unknown variables'
            err = err+1

        if len([v for v in self.Pfault if not v in np.arange(0,self.nx())])>0:
            print 'Error! Possible sensor locations with faults outside set of unknown variables'
            err = err+1

        if np.any(np.sum(self.F>0,axis=0)>1):
            print 'Error! Fault variables can only appear in 1 equation, rewrite model with intermediate variables'
            err = err+1;

        xIdx = np.where(np.all(self.X.toarray()==0,axis=0))[0]
        for ii in xIdx:
            print 'Warning! Variable ' + self.x[ii] + ' is not included in model'
            war = war + 1;

        zIdx = np.where(np.all(self.Z.toarray()==0,axis=0))[0]
        for ii in zIdx:
            print 'Warning! Variable ' + self.z[ii] + ' is not included in model'
            war = war + 1;

        fIdx = np.where(np.all(self.F.toarray()==0,axis=0))[0]
        for ii in fIdx:
            print 'Warning! Variable ' + self.f[ii] + ' is not included in model'
            war = war + 1;
        if self.IsUnderdetermined():
            print 'Warning! Model is underdetermined'
            war = war + 1;

        print '  Model validation finished with %d errors and %d warnings.' % (err, war)

        
def DiffConstraint(dvar,ivar):
    return [dvar, ivar, "diff"];

def _ModelStructure(rels,x) : 
    ne=len(rels)
    nx=len(x)

    X = np.zeros((ne,nx),dtype='int64')
    for k,rel in enumerate(rels):
        if IsDifferentialConstraint(rel):
            if (rel[0] in x) and (rel[1] in x):
                dvIdx = x.index(rel[0])
                ivIdx = x.index(rel[1])
                X[k,dvIdx] = 3;
                X[k,ivIdx] = 2;
        else:        
            X[k, _RelModelStructure(rel, x)]=1
    return sp.csc_matrix(X)

def _RelModelStructure(rel,x):
    if IsSymbolic(rel):
        relVars = [str(e) for e in rel.atoms(sym.Symbol)]
    else:
        relVars = rel;
    return [xi for xi in range(0,len(x)) if x[xi] in relVars]

def IsDifferentialConstraint( rel ):
    return isinstance(rel,list) and len(rel)==3 and rel[2] is "diff";

def IsSymbolic(v):
    return isinstance(v, tuple(sym.core.all_classes))

def _ToEquations(rels):
    def _ToEquation(rel):
        if IsSymbolic(rel) and not isinstance(rel,sym.Equality):
            return sym.Eq(rel)
        else:
            return rel
    return map(lambda r: _ToEquation(r), rels)
