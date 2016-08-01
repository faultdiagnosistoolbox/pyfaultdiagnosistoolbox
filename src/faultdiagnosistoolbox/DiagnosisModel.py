import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import sympy as sym

import dmperm as dmperm
import Matching as match
import StructurePlotting as smplot
import CodeGeneration as codegen

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
        
        if modeldef['type'] is 'VarStruc' or modeldef['type'] is 'Symbolic':
            self.X = _ModelStructure( modeldef['rels'], modeldef['x'])
            self.x = modeldef['x']
            self.F = _ModelStructure( modeldef['rels'], modeldef['f'])
            self.f = modeldef['f']
            self.Z = _ModelStructure( modeldef['rels'], modeldef['z'])
            self.z = modeldef['z']
            self.e = map(lambda x:"e"+np.str(x+1),np.arange(0,self.ne()))
            self.type = modeldef['type']

            if 'parameters' in modeldef:
                self.parameters = modeldef['parameters']
        elif modeldef['type'] is 'MatrixStruc':
            self.X = sp.csc_matrix(modeldef['X'])
            self.F = sp.csc_matrix(modeldef['F'])
            self.Z = sp.csc_matrix(modeldef['Z'])

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
            self.type = modeldef['type']
        else:
            print 'Model definition type ' + modeldef['type'] + ' is not supported (yet)'

        if modeldef['type'] is 'Symbolic':
            self.syme = np.array(_ToEquations(modeldef['rels']))

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
